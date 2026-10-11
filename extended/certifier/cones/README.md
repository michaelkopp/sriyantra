# The cone certifier

The certificate of Theorem A (Sections 5 and 7 of the paper, Table 2): the variety V8 of the extension
conditions has no point in the laminar set L*, anywhere in the configuration space of the thirteen points.
Python 3 with NumPy and SciPy (HiGHS).

## Coordinates

Gauge y7 = ∞ (as in `bracketlp.py`). The twelve finite points lie on the line in the order

    P0..P11 = y8, y9, y10, zU, zD, y0, y1, y2, y3, y4, y5, y6

and gap k (k = 0..10) is the distance P_k → P_{k+1}, with logarithm l_k. The scale is fixed by l_root = 0.
Every configuration of the 13 points (in their cyclic order, modulo projective maps) is one sequence l, and every
degeneration (points merging, scales separating) is a direction to infinity.

**Cones.** The Cartesian tree of l has the largest gap at its root and the Cartesian trees of the gaps to its
left and right as subtrees. For a fixed binary tree T, the sequences with Cartesian tree T form the cone

    s_k = l_k − l_parent(k) ≤ 0   for the ten non-root nodes k.

The 58,786 binary trees on eleven gaps give cones that cover the whole space. Inside a cone the ten coordinates
s_k are independent, so every degeneration is s → −∞ in some coordinates: a half-infinite box.

**Log-distances (Lemma 5.1).** For points P_a, P_b the gaps a..b−1 form a range whose largest gap is its lowest
common ancestor m in T, and

    L(P_a, P_b) = l_m + x_ab,     x_ab = log(1 + Σ exp(l_k − l_m))  (k in the range, k ≠ m),

where l_k − l_m is the sum of s along the path from m to k. So every log-distance is a linear form in s plus a
correction x_ab whose enclosure is exact as s → −∞. Conditions, clauses and the fan coordinates are linear in the
log-distances, hence linear in s plus a sum of corrections.

## Closing a box

A box (cone T, intervals for s) is closed when

1. an interval test refutes a condition or a clause (linear part exact, corrections as intervals): test (a);
2. an LP certificate refutes it: test (b);
3. the mean-value relaxation of `mvf.py` refutes it, possibly after certified contraction: tests (c) and (d)
   (in the whole-space mode, `CELLS=all`);
4. its fan coordinates lie in the region D of Proposition 6.2 (certified in `../fan/`): test (e).

In the run restricted to the six cells of `../cells/` (environment variable `CELLS`, the default of `cone.py`),
a box is also closed when its enclosure of the gap ratios misses the six cells, and the LP carries the rows of
the cells it meets (one LP over the smallest box of gap ratios containing them, or one LP per cell when that fails
and the box meets at most two).

The LP has the variables (s, x): ten cone coordinates with the box as bounds, and the 66 corrections with their
enclosures as bounds. Every log-distance is substituted as L = l_m + x, exactly: the s-coefficients are computed in
integer or rational arithmetic, and the tangent and chord rows of the three-point relations
L_ik = L_ij + sp(L_jk − L_ij) carry exactly complementary coefficients (q and 1 − q, q rational), so every row is a
valid inequality on the whole box, half-infinite boxes included. All coefficients are dyadic rationals (slopes
rounded to multiples of 2^-30, with the constants recomputed for the rounded slope: the minimum of sp(u) − s u for
tangents, a line above both end points for chords), so the substitution is exact in floating point.
Infeasibility is shown by Farkas multipliers that are checked in exact rational arithmetic.

## Mean-value relaxation and certified contraction (`mvf.py`)

A box that the LP does not close is given a second relaxation, formed directly in the cone coordinates. Every
condition and clause f is a smooth function of s. A coordinate whose interval is half-infinite or wider than 1 is
replaced inside the corrections by u_i = exp(s_i) ∈ [e^a, e^b] (or [0, e^b]), because the corrections depend on s_i
only through products of such exponentials; then f = k·s_U + h(z), z = (s_F, u), with k the integer coefficients of
those coordinates in the linear part (kept exact) and h smooth on the box, up to u = 0. With the interval gradient G
of h over the box, its centre z_c, half-widths r and any dyadic vector g0 (the gradient at z_c, rounded to a
multiple of 2^-30),

    h(z) − g0·z ∈ [h(z_c) − g0·z_c − R, h(z_c) − g0·z_c + R],   R = Σ_i max|G_i − g0_i| r_i

(h(z_c) and R enlarged by a relative 1e-12). A condition gives two rows, a clause one. Tangent lines of exp
(u ≥ α s + α(1 − log α), α dyadic, constant rounded down) and, on a finite interval, the chord (constant rounded up)
link u_i with s_i. The box holds no point of V8 ∩ L* when this system is infeasible, which is shown by Farkas
multipliers checked in exact rational arithmetic.

**Contraction.** When the system is feasible, each coordinate (s_i, and u_i) is minimised and maximised over it by LP.
A slab {x_j ≤ β} with β below the minimum by a margin is cut off the box when the system together with x_j ≤ β has an
exact Farkas certificate, and likewise at the upper end; bounds on u_i give bounds on s_i through s = log u, rounded
outward. Every point of V8 ∩ L* in the box therefore lies in the contracted box, and the branch and bound continues
with it. Up to four rounds are made, each with the relaxation of the current box.

**Trinomials whose terms differ on the u-coordinates.** For log(e^ma + e^mb) with different coefficients of the
u-coordinates in ma and mb, two variables are added: δ = mb − ma (a mean-value row, with the enclosure of mb − ma as
bounds) and w = sp(δ), sp(t) = log(1 + e^t), bounded below by tangents of sp and above by its chord (both with dyadic
slopes, constants rounded outward); the trinomial is then ma + w, with a mean-value row for ma.

**Gradient enclosures.** For a correction x = log(1 + Σ_j e^{w_j}), ∂x/∂s_i = E_J/(1 + E_J + E_N) with J the paths
through i and E the sums of the exponentials; it increases with E_J and decreases with E_N, so the enclosures of E_J and
E_N bound it. Likewise ∂x/∂u_i = F_J/(1 + u_i F_J + E_N), F_J = Σ_J e^{w_j − s_i}.

## Branching

A box that is not closed is split in one coordinate. After a failed mean-value relaxation it is the coordinate with
the largest share Σ max|G_i − g0_i| r_i of the remainders. Otherwise the coordinate is chosen where the chord of exp is
loosest: for every correction term exp(w), w = sum of s along a path, the largest gap between exp and its chord over
the enclosure of w is shared among the coordinates of the path in proportion to their widths (a half-infinite
coordinate takes it whole), and the coordinate with the largest share is split. A finite interval is split where the
range of exp halves; a half-infinite interval (−∞, a] at a − max(2, |a|).

## Files

* `cone.py` — trees, enclosures, interval tests, the LP and its exact check, the region D.
* `bracketlp.py` — the conditions and clauses as linear forms in the 66 log-distances (from `tri.json`), the
  relaxation rows and the exact verification of Farkas certificates.
* `mvf.py` — the mean-value relaxation and the certified contraction.
* `bbcone.py` — branch and bound inside one cone.
* `roots_all.py` — tests all 58,786 cones on their root box (s ≤ 0) and writes the 49,627 that remain to
  `roots_all_surv.json`. `roots.py` does the same for the run restricted to the six cells (`roots_surv.json`).
* `runcones.py` — runs the cones of a root file (resumable, any number of workers; `REDO=1` runs again the cones
  whose last record has open boxes, with a larger budget).
* `run_all.sh` — the whole configuration space (`CELLS=all`) on N workers: runs five of the tests below, then the
  workers; with `REDO=1` the cones still open are shared out through a claim queue, and `NOREGIONS=0` lets boxes
  close in D.
* `summarize.py` — checks that every cone of a root file has a closed record.
* `test_cone.py`, `test_cone2.py`, `test_mvf.py`, `test_mvf_tail.py`, `test_contract.py`, `test_split.py`,
  `test_grad.py`, `test_comb.py` — soundness checks; each must print PASS. They build boxes (finite, half-infinite
  and wide coordinates) around Newton-polished points of V8 (equations, no clauses; `v8_t*.json`) and of the
  laminar part of V7 (clauses, no target; `pts_t*.json`): none may be refuted, and the contracted box must still
  contain the point. `test_grad.py` checks the gradient enclosures against central differences and closed forms at
  random points of random boxes.
* `tri.json` — the conditions and clauses in the thirteen points (Tables 5 and 6 of the paper are generated from it
  by `../../check/gentables.py`).

## Running

    python3 test_cone.py              # each of the eight tests prints PASS
    python3 test_cone2.py 3 40
    python3 test_mvf.py 1 20
    python3 test_mvf_tail.py
    python3 test_contract.py 1 6
    python3 test_split.py
    python3 test_grad.py 1 20
    python3 test_comb.py

    python3 roots_all.py              # 58,786 cones: 9,159 closed at the root, 49,627 to run
    ./run_all.sh 10                                   # 5,000 boxes per cone
    NOREGIONS=0 BUDGET=300000 REDO=1 ./run_all.sh 10  # the cones left open, with a larger budget
    python3 summarize.py roots_all_surv.json all_*.jsonl

`run_all.sh` is resumable: run it again after an interruption and the finished cones are kept.

## Records

`all_0.jsonl` … `all_9.jsonl` hold one record for every run of a cone (the ten workers of `run_all.sh`, first with
the budget of 5,000 boxes per cone, then the cones left open with `REDO=1`), and `all_*.log` their logs.
`cloud_all_0.jsonl`, `cloud_all_1.jsonl` are the records of a second run started from the other end of the list on
another machine; `run_all.sh` skips the cones closed there, and every cone also has a closed record in `all_*.jsonl`.

    python3 summarize.py roots_all_surv.json all_*.jsonl

    cones 49627 of 49627 done, 0 missing, 0 with open boxes; 4522045 boxes, 253924 s CPU
    closing reasons: {'LP': 1996263, 'iv': 154707, 'MVF': 132058, 'D': 2808}
    LP calls 4364969, MVF calls 912937, slabs cut off by certified contraction 2259157
    cones closed with the help of region D: 16 [22585, 22586, 22599, 22600, 22717, 22718, 22731, 22732, 23014, 23015, 23028, 23029, 23146, 23147, 23160, 23161]
    CLOSED

These are the numbers of Table 2. With the 9,159 cones closed at the root, every one of the 58,786 cones is closed,
so V8 ∩ L* is empty.

## The six cells

The second certification (`../cells/`) leaves six cells to this certifier, run with the default `CELLS`:

    python3 roots.py                  # 21,877 cones outside the six cells, 4,299 closed at the root, 32,610 to run
    python3 runcones.py roots_surv.json W N six_W.jsonl 5000     # W = 0..N-1
    python3 summarize.py roots_surv.json six_0.jsonl six_1.jsonl

    cones 32610 of 32610 done, 0 missing, 0 with open boxes; 262314 boxes, 22241 s CPU
    closing reasons: {'LP': 120586, 'iv': 23697, 'cell': 1220, 'D': 1959}
    LP calls 235438, MVF calls 0, slabs cut off by certified contraction 0
    cones closed with the help of region D: 16 [22585, 22586, 22599, 22600, 22717, 22718, 22731, 22732, 23014, 23015, 23028, 23029, 23146, 23147, 23160, 23161]
    CLOSED

(`cell`: the box misses the six cells. The mean-value relaxation is used only in the whole-space mode.)
