# NIRS mVO₂ Recovery Analyzer

A browser tool, and a matching Python script, for measuring skeletal-muscle oxygen consumption (mVO₂) recovery from repeated short arterial occlusions recorded with an Artinis PortaLite and exported from OxySoft.

**Open the tool:** https://jlgagnon-glitch.github.io/NIRS/ (once GitHub Pages is enabled)

Everything runs in your own browser. The spreadsheets you load are never uploaded anywhere.

## Using the browser tool

1. Open the link above in Chrome, Safari, Edge or Firefox. You can also download `index.html` and open it from your computer; it works offline.
2. Drag your OxySoft `.xlsx` exports onto the page, or click **Add .xlsx files**. Each file name must contain `subject` followed by the subject ID, for example `subject0103.xlsx`.
3. Optionally click **Load protocol CSV** to set the protocol marker times yourself (see `protocol_template.csv`). Without it, the tool works out the times from the OxySoft event markers.
4. Check each subject. Every value starts at the default, and you only adjust it if the data call for it:
   - drag the protocol marker lines on the whole-recording graph;
   - drag the orange handles to change a cuff's slope window, or the grey dashed handles to move its on/off edges;
   - **Exclude** a bad cuff, or exclude a whole subject with a note.

   Edited values are highlighted, and **Reset** returns them to the defaults.
5. **Export CSVs** gives `slopes.csv`, `fits.csv`, `group_summary.csv` and `subjects.csv` (resting and end-of-exercise mVO₂). **Print / PDF report** prints a summary plus two pages per subject; choose *Save as PDF* in the print dialog. **Save session** stores your edits in a file so the analysis can be reproduced or shared.

## Expected protocol

For each series: a 5-minute arterial occlusion, then six 5-s occlusions (5 s off, 5 s on, repeated) starting at the cuff release. The protocol has one series before and one after a hand-grip exercise. The default channel is Rx1–Tx3 (OxySoft columns 6 and 7), with event codes and text in columns 16 and 17.

## Method (defaults)

- **Blood-volume correction:** Ryan et al. (2012), with an absolute-value β computed at every sample: β = |ΔO₂Hb| / (|ΔO₂Hb| + |ΔHHb|) and HHb꜀ = ΔHHb − ΔtHb·β. Δ values are measured from the 0.2 s before each cuff inflation.
- **Finding the cuffs:** a 5-s on / 5-s off grid is aligned to the sharp tHb steps within −2 to +3 s of the series marker. Each edge is then refined to the nearest step (±0.6 s).
- **mVO₂:** the least-squares slope of corrected HHb from 1.0 s after cuff-on to 0.3 s before cuff-off (µM/s).
- **Resting mVO₂:** the same slope from 1.0 s to 30.5 s after the first 5-minute cuff goes on.
- **End-of-exercise mVO₂:** the same slope from 1.0 s to 4.7 s after the second 5-minute cuff goes on (as the hand grip stops), i.e. the window of a 5-s cuff.
- **Recovery model:** y(t) = End − Delta·e^(−t/k), with t = 0 at the first 5-s cuff. The rate constant is 60/k (min⁻¹). Three fits are made:
  - 1-parameter: End = resting mVO₂, Delta = End − the first point;
  - 2-parameter: End = resting mVO₂;
  - 3-parameter: all three parameters free;
  - anchored (post-exercise series only): End = resting mVO₂, and the curve passes through the end-of-exercise mVO₂ at the release of the second 5-minute cuff (recovery is assumed to pause while the cuff is on), so only k is fitted.

  k is bounded to 0.5–2000 s.
- **Usable fit:** k between 5 and 1990 s, R² ≥ 0.80, no negative mVO₂ slope among the included cuffs, and the subject not excluded.

All of these can be changed under **Settings**.

## Python version

`nirs_mvo2.py` runs the same analysis from the command line. On default settings it gives the same results as the browser tool.

```
python3 -m pip install -r requirements.txt
python3 nirs_mvo2.py --data DATA_FOLDER --protocol protocol.csv --out results
```

`protocol.csv` needs one row per subject with the columns `subject, occl1_start, series1_start, occl2_start, series2_start, exclude, note`, as in `protocol_template.csv`.

## Reference

Ryan TE, Erickson ML, Brizendine JT, Young H-J, McCully KK. Noninvasive evaluation of skeletal muscle mitochondrial capacity with near-infrared spectroscopy: correcting for blood volume changes. *J Appl Physiol* 113: 175–183, 2012.
