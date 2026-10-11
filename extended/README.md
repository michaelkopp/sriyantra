# Extended circular-arc Śrī Yantras: programs and data

Everything the paper *Extended circular-arc Śrī Yantras* (`../sri-yantra-extended.pdf`) uses, as described in its
Section 10. The two constructed figures of Theorem B are also in the viewer (`../index.html`, *Arcs extended 2* and
*Arcs extended 3*) and in the GeoGebra files `../sri-yantra-arcs-extended-2.ggb` and
`../sri-yantra-arcs-extended-3.ggb`.

| folder | contents | paper |
|---|---|---|
| `certifier/cones/` | the certificate of Theorem A: branch and bound over the 58,786 cones, the soundness tests, and the records of the run | Sections 5, 7; Table 2 |
| `certifier/fan/` | the certificate of the regions B1, B2, B3 and D of the fan chart around the vertex T* | Section 6; Proposition 6.2 |
| `certifier/cells/` | the second certification of Theorem A by 32,805 cells, with its records | remark in Section 7 |
| `check/` | exact checks of the algebra, the generated tables, the figures of Part I, and `check/arcs/`, the computation in the proof of Theorem B | Lemma 3.1, Proposition 6.1, Section 9, Tables 4–8, Figures 1–2 |
| `paper/` | the LaTeX source and the figures | |

Each folder has a README with its files, how to run them and what they print.

## Requirements

* Python 3; NumPy and SciPy (HiGHS) for the two certifiers' linear programs; NumPy and Matplotlib for the figures.
  The checks in `check/` and `check/arcs/` need the standard library only.
* Node.js for `certifier/fan/` and the branch and bound of `certifier/cells/`.

## Quick checks

    python3 check/arcs/certify.py                                  # Theorem B, about a second: ALL CERTIFIED
    (cd certifier/cones && python3 test_cone.py)                   # one of the eight soundness tests: PASS
    (cd certifier/cones && python3 summarize.py roots_all_surv.json all_*.jsonl)   # the records of Table 2: CLOSED
    (cd certifier/cells && python3 summarize.py)                   # the second certification: CLOSED
    (cd certifier/fan && node fanRBB.js 0.1 3e6 '[0,[2,inf],0,0,[-2,2],[-2,2],[-2,2],[-2,2],[-2,2],[-2,2]]')   # B1: "open":0

The full runs are given in the READMEs: the cone certifier took 71 hours of processor time on ten cores, the region D
and the box B3 of the fan chart several minutes each.
