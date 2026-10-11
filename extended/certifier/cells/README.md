# The second certification: cells

An independent certificate of Theorem A, by a decomposition into cells (the remark at the end of Section 7 of the
paper). Outside six cells it shares no decomposition with `../cones/`; the six cells are closed there.
Node.js for the branch and bound, Python 3 with NumPy and SciPy (HiGHS) for the LP certificates.

## Cells

Gauge y7 = ∞. On the line the twelve finite points are, in order,

    y8 < y9 < y10 < zU < zD < y0 < y1 < y2 < y3 < y4 < y5 < y6.

Six of them, y8, y0, y2, y3, y4, y5, are known through the five gaps between them. In chart j (j = 0, …, 4) the j-th
of these gaps is the largest and is the unit, and the logarithms of the four others lie in [−∞, 0]. Each of the four
is cut at −8, −6, −4, −3, −2, −1.5, −1, −0.5 into nine bins, so the cells are 5 × 9^4 = 32,805 boxes that cover the
configuration space. A cell is named `chart.i.j.k.l` by its bins.

In a box of a cell, `g8.js` solves for the other six points in turn, each from a condition affine in its gap: y9
(Q43), zU (apexU), y1 (DOWN), y6 (Q21), y10 (Q59), zD (apexD). The remaining equations (the target N9 = P9 and UP)
and the 25 clauses are then tested. Every quantity is the logarithm of a positive number, enclosed in an interval
with outward rounding (`lg.js`); differences are taken only where their sign is known.

## Closing a box

A box is closed when

1. the interval tests of `g8.js` refute it (a solve step has no root in the range, or the target, UP or a clause
   has the wrong sign; `order …`, `target`, `UP`, `nc …`, `corner …` in the records), or
2. an LP certificate refutes it (`bracketlp.py`; `LP`), or
3. an LP certificate shows that every point of the box satisfying the relaxed conditions lies in one of the regions
   B1, B2, B3 of Proposition 6.2, certified in `../fan/` (`handoff`; one LP for each half-space outside the region).

**LP certificates.** Variables: L_ab = log|p_a − p_b| for the 66 pairs of finite points. Every condition and clause
is linear in L. For three points i < j < k on the line, L_ik = L_ij + sp(L_jk − L_ij), sp(u) = log(1 + e^u),
relaxed by a tangent line and the chord of sp over the enclosure of L_jk − L_ij. The trinomials e^m1 + e^m2 = 1 and
the corner clauses e^ma + e^mb ≤ 1 are relaxed by tangent lines of fixed dyadic slope. Bounds on L are the
enclosures chart A computes for the box. When the enclosure of L_jk − L_ij is bounded on one side only, the upper
bounds sp(u) ≤ sp(u1) (u ≤ u1) and sp(u) ≤ u + sp(−u0) (u ≥ u0) replace the chord. Infeasibility is certified by
Farkas multipliers from HiGHS, rounded to nearby rationals and checked exactly with Python Fractions:

    w = Aᵀy + Eᵀz,   Σ_j min_{x_j ∈ [l_j, u_j]} w_j x_j  >  y·b + z·c,   y ≥ 0.

**Splitting.** Besides the four gap ratios of the box, the branch and bound splits the range of a solved gap when
the solve step cannot decide whether its root lies in the range, and when the step is decided but its root
enclosure is still wider than max(0.05, 30 × box width). The pieces cover the range, so the split is exhaustive.

## Files

* `runcellsLP.js` — runs a set of cells on N workers, appends one record per cell to the output file, and skips
  the cells already there when resumed.
* `bbZLP.js` — branch and bound in one cell; `g8.js` — the solve steps and interval tests; `lg.js` — interval
  arithmetic on logarithms.
* `runcells.js`, `bbZ.js` — the same with the interval tests alone (no LP); they wrote `results2.jsonl`.
* `test_g8.js` — the interval tests at sample points (`v8_l.json`, points of V8; `lam_l.json`, laminar points;
  `etaconv.js` converts them to cell coordinates): the points of V8 are refuted by the clauses and by nothing
  else, and the laminar points by the target.
* `bracketlp.py` — the LP relaxation and the exact check of its certificates; `lpserverA.py` serves it to the
  workers.
* `test_lp.py` — soundness of the LP certificates on sample points: boxes around points of V8 are never refuted
  (clauses left out), boxes around laminar points are never refuted (target left out), the handoff succeeds around
  the laminar points with t3 < 0.19, and boxes around points of V8 with some enclosures made one-sided are never
  refuted. It must print PASS.
* `tri.json` (conditions and clauses), `u11.json`, `v8_t.json` and `pts_t.json` (sample points for the test),
  `budget44.txt`, `open4.txt`, `cells8.txt` (the cells of runs 4, 5 and 8).
* `results2.jsonl` … `results8.jsonl` — the records of the runs: one line per cell with its status, the number of
  boxes and how they were closed.
* `summarize.py` — checks that the records close every cell (below).

## Runs

    node test_g8.js
    python3 test_lp.py                                                                    # PASS
    node runcells.js --jobs 10 --budget 5e5                                               # run 2: all cells
    node runcellsLP.js --jobs 10 --budget 2e5 --c0 0 --from results2.jsonl --out results3.jsonl
    node runcellsLP.js --jobs 10 --budget 1e6 --c0 0 --ids $(cat budget44.txt) --out results4.jsonl
    node runcellsLP.js --jobs 10 --budget 4e6 --c0 0 --ids $(cat open4.txt) --out results5.jsonl
    node runcellsLP.js --jobs 10 --budget 2e6 --c0 0 --ids 0.5.5.5.7 --sub 2 --out results6.jsonl
    node runcellsLP.js --jobs 10 --budget 3e6 --c0 0 --box '0:[[-2,-1.75],[-2,-1.75],[-2,-1.75],[-1,-0.75]]' --boxid 0.5.5.5.7/0.0.0.0 --sub 2 --out results7.jsonl --redo 1
    node runcellsLP.js --jobs 3 --budget 6e6 --minw 1e-4 --c0 0 --box '0:[[-2,-1.75],[-2,-1.75],[-2,-1.75],[-1,-0.75]]' --boxid 0.5.5.5.7/0.0.0.0 --sub 2 --out results7.jsonl --redo 1
    node runcellsLP.js --jobs 10 --budget 3e6 --c0 0 --ids $(cat cells8.txt) --out results8.jsonl

Run 2 delegates the boxes with u(F2) < 0.05 to the neighbourhood of T* (the default `--c0 0.05`); with `--c0 0`
nothing is delegated, and a record counts only when nothing was delegated. `--from FILE` runs the cells that FILE
records as open and the closed cells with delegated boxes; `--ids` runs the listed cells (`budget44.txt`,
`open4.txt` and `cells8.txt` are the cells left open by the run before); `--sub 2` splits each selected cell into
16 sub-cells, named `CELL/i.j.k.l`; `--redo 1` runs again the cells that the output file records as open; `--minw`
is the smallest box width. A run can be stopped and restarted; finished cells are kept. Runs 6 and 7 split the cell
0.5.5.5.7 into 16 sub-cells and its sub-cell 0.0.0.0 into 16 again. The six cells of `cells8.txt` that have no
record in `results8.jsonl` are those closed in `../cones/`.

## Result

    python3 summarize.py

    cells 32805; closed by their own records 32798, through sub-cells 1 ['0.5.5.5.7']; not closed here 6 ['0.0.0.0.0', '0.0.0.0.1', '0.0.0.0.2', '0.0.0.1.2', '0.0.1.2.3', '0.1.2.3.4']
    six cells: 32610 cones, 262314 boxes, 0 not closed
    CLOSED

`summarize.py` checks that the cells form the full grid and that in every chart the bins tile [−∞, 0]; that a
cell counts as closed only by a closed record with nothing delegated, or through sub-cells that tile it and are
closed; and that the six remaining cells, in which the four gap ratios are all extreme, are closed by the cone
certifier restricted to them (`../cones/six_*.jsonl`, see `../cones/README.md`, with the region D of `../fan/`).
