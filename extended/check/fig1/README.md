# Choosing the points of Figures 1 and 2

Figure 2 shows a member of the four-dimensional family of straight figures that are laminar, have both triples
of outer corners collinear, and satisfy every concurrence except (t2, t6, t6). The point was chosen for even rings
and wide clearances. The search, in the order it was run:

1. `maximin.py` — from the laminar points of the certifier (`../../certifier/cones/pts_t.json`), maximise the smallest
   relative clearance of the 25 laminarity inequalities. Writes `maximin.json`.
2. `front.py` — from the best of those, trace the largest smallest clearance as the miss is bounded by
   0.30 … 0.04 (tip clearances kept in [0.02, 0.12]). Writes `front.json`. This is the source of the numbers in the
   introduction (about 16% at a miss of 10%, about 4% at 4%).
3. `RINGOBJ=range python3 nm.py 0.06 0.17,0.14,0.12,0.10 nmD.json` — minimise the spread of cell areas in the five
   rings (sum over rings of log(largest/smallest)), with clearances at least 0.06 and the miss bounded.
4. `RINGOBJ=range python3 nmtip.py` — the same at miss <= 0.14 with tip clearances in [0.05, 0.10].
5. `RINGOBJ=range python3 nmtip2.py` — the same with the apexes of U0 and D0 tied to the end rungs,
   y0 - zD = y2 - y1 and zU - y10 = y10 - y9 (two more equations; the family is then two-dimensional).
6. `RINGOBJ=range python3 tip2b.py` — refine at misses 0.14, 0.15, 0.16; the figure is the one at 0.15.
7. `pick.py` — polish it onto the family and write `../straight_pt.json`.
8. The yantra of the nested-kite figure (`../classical_pt.json`) is Chiodo's straight construction at the
   viewer parameters of the preset `Fig1paper2` (exactp = 1.31067, 0.87705, 0.74617, 0.27759), in which the outer
   corners of U1, U2, U3 lie on one line (`fig1paper2_preset.json`). `viewer_fig1.py` reads its nine triangles from
   the viewer (`T9.json`) and `viewer_point.py` converts them to the ladder model and checks them. (`classical.py`
   projects the point of Figure 2 onto the classical figures instead; it is not used for the figure.)

Then `../straightfig.py` checks the point, `../drawstraight.py` draws `../../paper/figs/straight_nearmiss.{pdf,png}`
(Figure 2; panels (b), (c) carry triangles out with `walls.py`, in the positions of the wall patterns 22 and 36)
and `../drawkites.py` draws `../../paper/figs/kites.{pdf,png}` (Figure 1). Figure 3 is exported from a copy of the viewer
with the same tones: `../../paper/tools/add_paper_palette.py index.html index_paper.html`, then
`../../paper/tools/export_paper.py index_paper.html ../../paper/figs arcx2,arcx3` and `../../paper/tools/svg2.py` on the two SVGs.

`ringpaint.py` sorts the cells into rings (the cells at odd depth from the outside: in the extended figure the new
ring is depth 1, the tips of U0 and D0 at their apexes and the four triangles at the corners of their bases) and the
cells between rings (a pale tint of the innermost ring they touch), and builds the kites; `palette.py` holds the ring tones (OKLCH,
rainbow order from the centre, lightness rising outward so that the rings stay apart in greyscale).

`model.py` holds the ladder model (half-widths, residuals of E1–E4, C_D, C_U, the 25 laminarity inequalities, the
extension) and `arr.py` the planar arrangement (cells, their areas and their depth from the outside); the evenness
of the rings is measured on the cells at depths 2, 3, 5, 7, 9 of the eleven-triangle arrangement (the fourteen cells
between the yantra and the fifth kite, and the rings of 14, 10, 10, 8).
The json files are the outputs of steps 1–6, so any step can be rerun from the previous one.
