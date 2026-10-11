# The constructions of Theorem B

- `recipes.json` — the two constructions as the viewer (`index.html`) holds them, `ARC_REC2` (pattern 22) and
  `ARC_REC3` (pattern 36): placed corners, steps, branch hints and default corners.
- `construct.py` — evaluates a construction in floating point or in interval arithmetic and lists the incidences
  each step realises. At the viewer's default corners it reproduces the viewer's nodes to 3e-14.
- `iv.py` — interval arithmetic at 60 significant digits with outward rounding; cos and sin of a rational angle by
  the Taylor series with its remainder.
- `certify.py` — the computation in the proof of Theorem B, at the corners of Table 4 read as exact decimals:
  every meeting exists and its branch is decisive; the node-on-side incidences are exactly those realised by the
  steps (106 and 110), each once, every other node/side pair being separated; and no leg meets a leg of another
  triangle of its family (50 pairs per figure). Writes `certify.log`, `incidences_22.csv`, `incidences_36.csv`.
- `seedtable.py` — writes `../../paper/tex/seeds.tex`, Table 4: the placed corners and their values.
- `steptable.py` — writes `../../paper/tex/steps_arcs.tex`, the step tables of Appendix B.

Run `python3 certify.py` (standard library only; about a second), `python3 seedtable.py` and `python3 steptable.py`.
