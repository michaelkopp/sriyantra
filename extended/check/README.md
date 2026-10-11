# Checks, tables and figures

Python 3. `straight.py`, `fan.py`, `gentables.py` and `arcs/` need the standard library only; the figure scripts
need NumPy and Matplotlib.

| file | what it does | paper |
|---|---|---|
| `straight.py` | verifies the half-widths and the four quartics of the ladder model in exact rational arithmetic | Lemma 3.1 |
| `fan.py` | checks the identities for Φ numerically in the fan chart, from `../certifier/cones/tri.json` | Proposition 6.1 |
| `gentables.py` | writes `../paper/tex/tables13.tex`, the conditions and clauses in the thirteen points, from `tri.json`, the data the certifiers read | Tables 5, 6 |
| `straightfig.py` | checks the point of Figure 2 (`straight_pt.json`): the equations that hold, the one that fails, every laminarity inequality with its margin, the tip clearances and the spread of the cell areas in each ring | Figure 2 |
| `drawkites.py` | draws Figure 1 from `classical_pt.json` into `../paper/figs/kites.{pdf,png}` | Figure 1 |
| `drawstraight.py` | draws Figure 2 from `straight_pt.json` into `../paper/figs/straight_nearmiss.{pdf,png}` | Figure 2 |
| `fig1/` | the choice of the points of Figures 1 and 2, the ladder model, the planar arrangement, the rings and the palette; see `fig1/README.md` | Figures 1–3 |
| `arcs/` | the computation in the proof of Theorem B, and the tables of the two constructions; see `arcs/README.md` | Section 9, Tables 4, 7, 8 |

    python3 straight.py
    python3 fan.py
    python3 gentables.py
    python3 straightfig.py
    python3 drawkites.py
    python3 drawstraight.py
    python3 arcs/certify.py
    python3 arcs/seedtable.py
    python3 arcs/steptable.py
