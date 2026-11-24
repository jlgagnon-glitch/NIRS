#!/usr/bin/env python3
"""NIRS mVO2 recovery analysis.

Reads OxySoft .xlsx exports, applies the Ryan blood-volume correction,
measures mVO2 in each 5-s cuff occlusion, fits the recovery curve and
writes figures, CSV tables and a PDF report.

Usage:
    python3 nirs_mvo2.py --data DATA_FOLDER --protocol protocol.csv --out OUTPUT_FOLDER
"""

# ---------------------------------------------------------------- 1. Imports
import argparse
import csv
import re
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from openpyxl import load_workbook

# ---------------------------------------------------------------- 2. Settings
O2HB_COL = 6           # OxySoft column holding Rx1-Tx3 O2Hb (1-based)
HHB_COL = 7            # OxySoft column holding Rx1-Tx3 HHb
EVENT_COL = 16         # column holding event codes (A1, A2, ...)
TEXT_COL = 17          # column holding event text
N_CUFFS = 6            # 5-s occlusions per series
CYCLE = 10.0           # seconds from one cuff to the next (5 s off + 5 s on)
CUFF_DELAY = 5.0       # first cuff goes on 5 s after the series marker
REF_WINDOW = 0.2       # seconds before cuff-on used as the zero for the correction
SLOPE_START = 1.0      # start each slope this many seconds after cuff-on
SLOPE_TRIM = 0.3       # end each slope this many seconds before cuff-off
REST_START, REST_END = 1.0, 30.5   # resting-slope window after the first long cuff goes on
SHIFT_RANGE = (-2.0, 3.0)          # allowed offset (s) between marker and real cuff timing
K_BOUNDS = (0.5, 2000.0)           # limits on the time constant k (s)
USABLE_K = (5.0, 1990.0)           # k must fall in this range to be usable
USABLE_R2 = 0.80                   # and R^2 must be at least this
FITS = ["1-param", "2-param", "3-param", "anchored"]
SERIES = ["Series 1 (pre-exercise)", "Series 2 (post-exercise)"]


# ---------------------------------------------------------------- 3. Reading the data
def read_oxysoft(path):
    """Return time (s), O2Hb, HHb and the event list from one OxySoft export."""
    sheet = load_workbook(path, read_only=True, data_only=True)["Export 1"]
    rate, rows, events, in_data = 10.0, [], [], False
    for row in sheet.iter_rows(values_only=True):
        if not row or row[0] is None:
            continue
        if row[0] == "Export sample rate":
            rate = float(row[1])
        elif not in_data and list(row[:3]) == [1, 2, 3]:
            in_data = True                   # the row 1, 2, 3, ... marks the start of data
        elif in_data and isinstance(row[0], (int, float)) and row[O2HB_COL - 1] is not None:
            t = float(row[0]) / rate
            rows.append((t, float(row[O2HB_COL - 1]), float(row[HHB_COL - 1])))
            code = row[EVENT_COL - 1] if len(row) >= EVENT_COL else None
            if code is not None and t > 0.5:
                text = row[TEXT_COL - 1] if len(row) >= TEXT_COL else ""
                events.append((t, str(code).strip(), str(text or "").strip()))
    data = np.array(rows)
    return data[:, 0], data[:, 1], data[:, 2], events


def find_file(folder, subject):
    """Find the one .xlsx file in folder whose name contains 'subject' + the ID."""
    hits = [p for p in Path(folder).glob("*.xlsx")
            if re.search(r"subject" + subject, p.name, re.IGNORECASE)]
    if len(hits) != 1:
        raise SystemExit(f"Expected one file for subject {subject}, found {len(hits)}")
    return hits[0]


def read_protocol(path):
    """Read protocol.csv: one row per subject with its four marker times."""
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


# ---------------------------------------------------------------- 4. Signal processing
def ryan_correct(o2, hhb, ref):
    """Ryan correction with absolute-value beta; ref is a boolean mask for the zero."""
    d_o2 = o2 - o2[ref].mean()
    d_hhb = hhb - hhb[ref].mean()
    d_thb = d_o2 + d_hhb
    total = np.abs(d_o2) + np.abs(d_hhb)
    beta = np.divide(np.abs(d_o2), total, out=np.full_like(total, 0.5), where=total > 0)
    return d_o2 - d_thb * (1 - beta), d_hhb - d_thb * beta


def step_size(thb):
    """Size of the change in tHb over three samples (big at cuff on/off)."""
    step = np.zeros_like(thb)
    step[1:-2] = np.abs(thb[3:] - thb[:-3])
    return step


def nearest_step(t, step, t0, tol):
    """Time and size of the largest tHb step within +/- tol seconds of t0."""
    near = (t >= t0 - tol) & (t <= t0 + tol)
    i = np.argmax(np.where(near, step, -1))
    return t[i] + 0.1, step[i]


def locate_cuffs(t, step, start):
    """Fit the 5-s on / 5-s off grid to the tHb steps after a series marker."""
    def grid_score(shift):
        times = [start + shift + CUFF_DELAY + CYCLE * i + extra
                 for i in range(N_CUFFS) for extra in (0.0, CUFF_DELAY)]
        return sum(nearest_step(t, step, x, 0.4)[1] for x in times)
    shifts = np.arange(SHIFT_RANGE[0], SHIFT_RANGE[1] + 0.01, 0.1)
    best = max(shifts, key=grid_score)
    cuffs = []
    for i in range(N_CUFFS):
        on = nearest_step(t, step, start + best + CUFF_DELAY + CYCLE * i, 0.6)[0]
        off = nearest_step(t, step, start + best + CYCLE * (i + 1), 0.6)[0]
        cuffs.append((on, off))
    return best, cuffs


def slope(t, y, t0, t1):
    """Least-squares slope of y between t0 and t1, and its r^2."""
    w = (t >= t0) & (t <= t1)
    m, _ = np.polyfit(t[w], y[w], 1)
    r = np.corrcoef(t[w], y[w])[0, 1]
    return m, r ** 2


def corrected_slope(t, o2, hhb, on, t0, t1):
    """mVO2: slope of Ryan-corrected HHb, zeroed just before cuff-on."""
    ref = (t >= on - REF_WINDOW) & (t <= on)
    _, hhb_c = ryan_correct(o2, hhb, ref)
    return slope(t, hhb_c, t0, t1)


# ---------------------------------------------------------------- 5. Exponential fits
def model(t, end, delta, k):
    return end - delta * np.exp(-t / k)


def fit_recovery(t, y, rest, anchor=None):
    """1-, 2- and 3-parameter fits of y(t) = End - Delta*exp(-t/k), plus the anchored fit
    (End = rest and the curve passes through the end-of-exercise mVO2 at time t_ex)."""
    lo, hi = K_BOUNDS
    delta0 = rest - y[0]
    tries = {
        "1-param": (lambda t, k: model(t, rest, delta0, k), [20], (lo, hi)),
        "2-param": (lambda t, d, k: model(t, rest, d, k), [delta0, 20], ([-50, lo], [50, hi])),
        "3-param": (model, [y[-1], y[-1] - y[0], 15], ([-50, -50, lo], [50, 50, hi])),
    }
    if anchor is not None:
        t_ex, ex = anchor
        tries["anchored"] = (lambda t, k: rest - (rest - ex) * np.exp(-(t - t_ex) / k), [20], (lo, hi))
    fits = {}
    for name, (func, p0, bounds) in tries.items():
        try:
            p, _ = curve_fit(func, t, y, p0=p0, bounds=bounds, maxfev=50000)
        except (RuntimeError, ValueError):
            continue
        if name == "1-param":
            end, delta, k = rest, delta0, p[0]
        elif name == "2-param":
            end, delta, k = rest, p[0], p[1]
        elif name == "anchored":
            end, delta, k = rest, (rest - ex) * np.exp(t_ex / p[0]), p[0]
        else:
            end, delta, k = p
        resid = y - model(t, end, delta, k)
        r2 = 1 - np.sum(resid ** 2) / np.sum((y - y.mean()) ** 2)
        fits[name] = dict(End=end, Delta=delta, k=k, rate=60 / k, R2=r2)
    return fits


def is_usable(fit, slopes):
    return (USABLE_K[0] <= fit["k"] < USABLE_K[1] and fit["R2"] >= USABLE_R2
            and min(slopes) >= 0)


# ---------------------------------------------------------------- 6. One subject
def analyze(t, o2, hhb, marks):
    """Resting mVO2, the twelve 5-s slopes and the fits for one subject."""
    thb = o2 + hhb
    step = step_size(thb)
    cuff1 = nearest_step(t, step, marks["occl1_start"], 2.0)[0]
    rest, _ = corrected_slope(t, o2, hhb, cuff1, cuff1 + REST_START, cuff1 + REST_END)
    cuff2 = nearest_step(t, step, marks["occl2_start"], 2.0)[0]
    end_ex, _ = corrected_slope(t, o2, hhb, cuff2, cuff2 + SLOPE_START, cuff2 + CUFF_DELAY - SLOPE_TRIM)
    result = dict(rest=rest, cuff1=cuff1, cuff2=cuff2, end_ex=end_ex, series={})
    for name, key in zip(SERIES, ["series1_start", "series2_start"]):
        shift, cuffs = locate_cuffs(t, step, marks[key])
        occ = []
        for on, off in cuffs:
            m, r2 = corrected_slope(t, o2, hhb, on, on + SLOPE_START, off - SLOPE_TRIM)
            occ.append(dict(on=on, off=off, slope=m, r2=r2))
        mids = np.array([(c["on"] + c["off"]) / 2 for c in occ])
        tt = mids - mids[0]
        y = np.array([c["slope"] for c in occ])
        anchor = None
        if name == SERIES[1]:   # release of the 2nd long cuff, on the same clock as tt
            anchor = (occ[0]["on"] - CUFF_DELAY - mids[0], end_ex)
        result["series"][name] = dict(shift=shift, occ=occ, t=tt, y=y, anchor=anchor,
                                      fits=fit_recovery(tt, y, rest, anchor))
    return result


# ---------------------------------------------------------------- 7. Figures
def plot_fits(sid, res, path):
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)
    for ax, name in zip(axes, SERIES):
        s = res["series"][name]
        tt = np.linspace(-10, 52, 300)
        ax.axhline(res["rest"], color="grey", ls=":", label="resting mVO2")
        for fname, color in zip(FITS, ["tab:blue", "tab:orange", "tab:green", "tab:purple"]):
            if fname in s["fits"]:
                f = s["fits"][fname]
                ax.plot(tt, model(tt, f["End"], f["Delta"], f["k"]), color=color,
                        label=f"{fname}: k = {f['k']:.1f} s, R² = {f['R2']:.2f}")
        ax.plot(s["t"], s["y"], "ko", label="mVO2 (slope of each 5-s cuff)")
        if s["anchor"] is not None:
            ax.plot(*s["anchor"], "D", color="tab:purple", label="end-of-exercise mVO2")
        ax.set_title(name)
        ax.set_xlabel("Time since first 5-s cuff (s)")
        ax.legend(fontsize=8)
    ys = np.concatenate([res["series"][n]["y"] for n in SERIES] + [[res["rest"], res["end_ex"], 0]])
    pad = 0.2 * (ys.max() - ys.min())
    axes[0].set_ylim(ys.min() - pad, ys.max() + 2 * pad)
    axes[0].set_ylabel("mVO2 (µM/s)")
    fig.suptitle(f"Subject {sid}: exponential fits")
    fig.savefig(path, dpi=110, bbox_inches="tight")
    plt.close(fig)


def plot_windows(sid, t, o2, hhb, res, path):
    fig, axes = plt.subplots(2, N_CUFFS, figsize=(16, 6.5))
    for row, name in enumerate(SERIES):
        for i, c in enumerate(res["series"][name]["occ"]):
            ax = axes[row, i]
            on, off = c["on"], c["off"]
            _, hhb_c = ryan_correct(o2, hhb, (t >= on - REF_WINDOW) & (t <= on))
            near = (t >= on - 1.5) & (t <= off + 1.5)
            w = (t >= on + SLOPE_START) & (t <= off - SLOPE_TRIM)
            ax.axvspan(0, off - on, color="grey", alpha=0.15)
            ax.plot(t[near] - on, hhb_c[near], ".-", ms=2)
            line = np.polyval(np.polyfit(t[w], hhb_c[w], 1), t[w])
            ax.plot(t[w] - on, line, color="tab:orange", lw=2)
            ax.set_title(f"#{i + 1} {c['slope']:.3f} µM/s, r² {c['r2']:.2f}", fontsize=8)
        axes[row, 0].set_ylabel(name, fontsize=8)
    fig.suptitle(f"Subject {sid}: corrected HHb in each 5-s cuff (grey = cuff on)")
    fig.tight_layout()
    fig.savefig(path, dpi=100)
    plt.close(fig)


def plot_signals(sid, t, o2, hhb, res, events, path):
    c1 = res["cuff1"]
    o2_c, hhb_c = ryan_correct(o2, hhb, (t >= c1 - 10.5) & (t <= c1 - 0.5))
    fig, axes = plt.subplots(3, 1, figsize=(14, 11))
    axes[0].plot(t, o2_c, color="tab:red", lw=0.8, label="O2Hb (corrected)")
    axes[0].plot(t, hhb_c, color="tab:blue", lw=0.8, label="HHb (corrected)")
    for time, code, text in events:
        axes[0].axvline(time, color="grey", ls="--", lw=0.7)
        axes[0].text(time, axes[0].get_ylim()[1], f" {code} {text}", rotation=90,
                     va="top", fontsize=7)
    axes[0].set_title("Whole recording (zero = 10 s before the first long cuff)")
    axes[0].legend(loc="lower left")
    for ax, name in zip(axes[1:], SERIES):
        occ = res["series"][name]["occ"]
        near = (t >= occ[0]["on"] - 8) & (t <= occ[-1]["off"] + 4)
        for c in occ:
            ax.axvspan(c["on"], c["off"], color="grey", alpha=0.15)
        ax.plot(t[near], o2_c[near], color="tab:red", lw=1)
        ax.plot(t[near], hhb_c[near], color="tab:blue", lw=1)
        ax.set_title(name + " (grey = cuff on)")
    axes[-1].set_xlabel("Time (s)")
    fig.suptitle(f"Subject {sid}: Ryan-corrected O2Hb and HHb, Rx1-Tx3")
    fig.tight_layout()
    fig.savefig(path, dpi=100)
    plt.close(fig)


# ---------------------------------------------------------------- 8. Tables and report
def write_tables(results, out):
    with open(out / "subjects.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["subject", "excluded", "rest_mVO2", "end_exercise_mVO2", "note"])
        for sid, res in results.items():
            w.writerow([sid, res["excluded"], f"{res['rest']:.5f}", f"{res['end_ex']:.5f}", res["note"]])
    with open(out / "slopes.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["subject", "series", "cuff", "on_s", "off_s", "mVO2_uM_per_s", "r2"])
        for sid, res in results.items():
            for name in SERIES:
                for i, c in enumerate(res["series"][name]["occ"]):
                    w.writerow([sid, name, i + 1, f"{c['on']:.1f}", f"{c['off']:.1f}",
                                f"{c['slope']:.5f}", f"{c['r2']:.3f}"])
    with open(out / "fits.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["subject", "excluded", "series", "fit", "End", "Delta", "k_s",
                    "rate_per_min", "R2", "usable"])
        for sid, res in results.items():
            for name in SERIES:
                s = res["series"][name]
                for fname, fit in s["fits"].items():
                    w.writerow([sid, res["excluded"], name, fname, f"{fit['End']:.5f}",
                                f"{fit['Delta']:.5f}", f"{fit['k']:.2f}", f"{fit['rate']:.3f}",
                                f"{fit['R2']:.3f}", is_usable(fit, s["y"]) and not res["excluded"]])


def write_report(results, out):
    from reportlab.lib.pagesizes import letter, landscape
    from reportlab.lib.styles import ParagraphStyle
    from reportlab.lib.units import inch
    from reportlab.pdfbase import pdfmetrics
    from reportlab.pdfbase.ttfonts import TTFont
    from reportlab.platypus import Image, PageBreak, Paragraph, SimpleDocTemplate, Table
    font = Path(matplotlib.get_data_path()) / "fonts" / "ttf" / "DejaVuSans.ttf"
    pdfmetrics.registerFont(TTFont("DejaVu", str(font)))
    body = ParagraphStyle("body", fontName="DejaVu", fontSize=9, leading=12)
    head = ParagraphStyle("head", fontName="DejaVu", fontSize=14, leading=18, spaceAfter=6)
    story = [Paragraph("NIRS mVO2 recovery analysis", head)]
    rows = [["Subject", "Rest (µM/s)", "End-ex. (µM/s)"] + [f"S2 {f} rate (R²)" for f in FITS] + ["Note"]]
    for sid, res in results.items():
        s2 = res["series"][SERIES[1]]
        cells = []
        for fname in FITS:
            fit = s2["fits"].get(fname)
            if fit is None:
                cells.append("-")
            else:
                ok = is_usable(fit, s2["y"]) and not res["excluded"]
                cells.append(f"{fit['rate']:.2f} ({fit['R2']:.2f})" + ("" if ok else " *"))
        rows.append([sid, f"{res['rest']:.4f}", f"{res['end_ex']:.4f}"] + cells + [res["note"]])
    table = Table([[Paragraph(str(c), body) for c in r] for r in rows],
                  colWidths=[0.6 * inch, 0.8 * inch, 0.8 * inch] + [1.15 * inch] * 4 + [2.2 * inch])
    story += [table, Paragraph("* = not usable (k outside 5-1990 s, R² below 0.80, a negative slope, "
                               "or subject excluded). Rates are in 1/min.", body)]
    for fname in FITS:
        rates = []
        for res in results.values():
            s2 = res["series"][SERIES[1]]
            fit = s2["fits"].get(fname)
            if fit and not res["excluded"] and is_usable(fit, s2["y"]):
                rates.append(fit["rate"])
        if len(rates) > 1:
            story.append(Paragraph(f"Series 2, {fname}: n = {len(rates)}, mean {np.mean(rates):.2f} "
                                   f"± {np.std(rates, ddof=1):.2f}, median {np.median(rates):.2f} /min", body))
    for sid in results:
        for kind in ["signals", "fits", "windows"]:
            story += [PageBreak(), Image(str(out / "figures" / f"{sid}_{kind}.png"),
                                         width=9 * inch, height=6.2 * inch, kind="proportional")]
    SimpleDocTemplate(str(out / "report.pdf"), pagesize=landscape(letter)).build(story)


# ---------------------------------------------------------------- 9. Main program
def main():
    parser = argparse.ArgumentParser(description="NIRS mVO2 recovery analysis")
    parser.add_argument("--data", required=True, help="folder with the OxySoft .xlsx files")
    parser.add_argument("--protocol", required=True, help="protocol.csv with marker times")
    parser.add_argument("--out", default="results", help="folder for the output")
    args = parser.parse_args()
    out = Path(args.out)
    (out / "figures").mkdir(parents=True, exist_ok=True)
    results = {}
    for row in read_protocol(args.protocol):
        sid = row["subject"].zfill(4)
        path = find_file(args.data, sid)
        print(f"Subject {sid}: {path.name}")
        t, o2, hhb, events = read_oxysoft(path)
        marks = {k: float(row[k]) for k in
                 ["occl1_start", "series1_start", "occl2_start", "series2_start"]}
        res = analyze(t, o2, hhb, marks)
        res["excluded"] = row.get("exclude", "").strip().lower() == "yes"
        res["note"] = row.get("note", "")
        results[sid] = res
        plot_signals(sid, t, o2, hhb, res, events, out / "figures" / f"{sid}_signals.png")
        plot_fits(sid, res, out / "figures" / f"{sid}_fits.png")
        plot_windows(sid, t, o2, hhb, res, out / "figures" / f"{sid}_windows.png")
    write_tables(results, out)
    write_report(results, out)
    print(f"Done. Results are in {out.resolve()}")


if __name__ == "__main__":
    main()
