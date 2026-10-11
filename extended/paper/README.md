# The paper

*Extended circular-arc Śrī Yantras*; the compiled paper is `../../sri-yantra-extended.pdf`.

* `tex/` — the LaTeX source. Build in `tex/` with `pdflatex paper.tex`, twice. Three of its inputs are generated:
  `tables13.tex` (Tables 5 and 6) by `../check/gentables.py`, `seeds.tex` (Table 4) by `../check/arcs/seedtable.py`
  and `steps_arcs.tex` (Tables 7 and 8) by `../check/arcs/steptable.py`.
* `figs/` — the figures. `kites` (Figure 1) is drawn by `../check/drawkites.py` and `straight_nearmiss` (Figure 2) by
  `../check/drawstraight.py`; `../check/fig1/README.md` describes how their points were chosen. `arcx2` and `arcx3`
  (Figure 3) are exported from the viewer `../../index.html` in the paper's tones.
* `tools/` — the export of Figure 3: `add_paper_palette.py` adds the paper's tones to a copy of the viewer,
  `export_paper.py` exports the two figures as SVG with headless Chromium (Playwright), and `svg2.py` converts them to
  PDF and PNG; `export_viewer.py` exports any figure of the viewer in one of its own palettes.

      python3 tools/add_paper_palette.py ../../index.html index_paper.html
      python3 tools/export_paper.py index_paper.html figs arcx2,arcx3
      python3 tools/svg2.py figs/arcx2.svg figs/arcx3.svg
