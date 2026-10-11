# The fan-chart certifier

The certificate of Proposition 6.2 of the paper: V8 ∩ L* has no point in the region R = B1 ∪ B2 ∪ B3 ∪ D of the
fan chart. Plain `node`, no dependencies. The cone certifier (`../cones/`) closes a box when its fan coordinates
lie in D, and the cell certifier (`../cells/`) when they lie in B1, B2 or B3.

## The fan chart

zD = 0, zU = ∞, y10 = 1, with heights p_k = t_k t_{k+1} ⋯ t_10, t_k ∈ (0, 1), k = 1, …, 10
(p_0 = zD, p_1 = y0, …, p_11 = y10). Every bracket factors exactly as

    p_b − p_a = (t_b ⋯ t_10) · g(a, b),     g(a, b) = 1 − t_a ⋯ t_{b−1}   (a ≥ 1),

so every condition and every clause is  f = Σ κ_k log t_k + Σ P_ab log g(a, b)  (`funcs.json`, `clauses.json`).
`python3 gapform.py` derives these factorisations from `tri.json` and prints them, with the closed form of Φ.

**Φ** (Proposition 6.1) = L_target − L_Q21 + L_Q59 + log(m1 + m2)_DOWN has no monomial part:

    Φ = log(1 + R) + Σ n_ab log g(a, b),     R = m2/m1 (DOWN),

    n = +1 on (1,8) (2,4) (2,5) (3,7) (4,10) (5,8) (7,9) (8,11) (9,10),
    n = −1 on (1,4) (2,8) (2,9) (3,5) (4,8) (5,7) (7,10) (8,10) (9,11).

At every point of V8, Φ = 0.

## Adapted coordinates (`fanR.js`)

    e = log t3,  d = log(t6/t3),  log t1,  log t10,
    r2 = log(t2/t3),  r4 = log(t4/t6),  r5 = log(t5 t6/t3),  r7 = log(t7 t6/t3),  r8 = log(t8/t3),  r9 = log(t9/t3).

The linear parts become Q21: r5, Q43: r7, Q59: r4 + r5 − r8, apexU: r8 + r9 − r4 − r5, apexD: r2 − r7 − r8,
DOWN.m1: r8, target: −r4, so the scale e cancels. Each of the six conditions Q21, Q43, Q59, apexU, apexD, DOWN is
solved for its coordinate by interval propagation. Φ/t3 is evaluated as a signed sum of the terms −n (P/t3) h(P)
and (R/t3) h̃(R), with h(P) = −log(1 − P)/P and h̃(R) = log(1 + R)/R, which stays accurate as t3 → 0. A box is
refuted when propagation empties it, or when the target, UP, a clause, or Φ/t3 has a sign on it. An upper tail of
a relative coordinate is continued in the absolute coordinate log t_k, and a coordinate with t_k close to 1 in
log(−log t_k) (cluster mode). Caps t_k ≤ θ, as in D, enter the propagation.

**Arithmetic** (`lg.js`): every quantity is the logarithm of a positive number, enclosed in an interval with
outward padding; sums are evaluated as log-sum-exp, and a difference is formed only where its sign is known, so
that clusters of points at any scale are evaluated without loss.

## Files

* `fanRBB.js` — branch and bound over one box of the adapted coordinates.
* `fanD.js` — the region D: covers the part of D outside B1 ∪ B2 by 25 pieces and closes each with `fanRBB.js`.
* `fanR.js` — the adapted coordinates, the propagation and the tests; `fanB.js` — the forms of the conditions,
  clauses and Φ in the fan chart; `lg.js` — the arithmetic.
* `funcs.json`, `clauses.json` — the conditions and the 25 clauses in the fan chart; `tri.json` — the same in the
  thirteen points, from which `gapform.py` derives the factorisations.
* `ntest.js` — enclosure test: at random configurations (some t_k within 1e-9 to 1e-2 of 1), every form enclosed
  over a box around the point must contain its value there.
* `v6test.js` — soundness sample for D: points of V6 (the six solved conditions) found by Newton; propagation with
  the caps of D must not refute small boxes around them.
* `certificates.log` — the output of the four runs below.

## The certificates

The first argument is the bound c on t3; the JSON list overrides the root box in the order
e, d, log t1, log t10, r2, r4, r5, r7, r8, r9 (0 keeps the default; t1 and t10 are arbitrary throughout).

    node fanRBB.js 0.1 3e6 '[0,[2,inf],0,0,[-2,2],[-2,2],[-2,2],[-2,2],[-2,2],[-2,2]]'                    # B1
    node fanRBB.js 0.1 3e6 '[0,[0.5,2.2],0,0,[-0.4,0.4],[-0.4,0.4],[-0.4,0.4],[-0.4,0.4],[-0.4,0.4],[-0.4,0.4]]'   # B2
    node fanRBB.js 0.2 2e7 '[[-2.3026,-1.6094],[1.0,inf],0,0,[-0.25,0.25],[-0.25,0.25],[-0.25,0.25],[-0.25,0.25],[-0.25,0.25],[-0.25,0.25]]'  # B3
    node fanD.js 0.1 0.95 1e6                                                                              # D

Each run of `fanRBB.js` must print `"open":0`, and `fanD.js` `"closed":true`:

    B1: {"boxes":3567,"refuted":1784,"open":0,...}
    B2: {"boxes":81685,"refuted":40843,"open":0,...}
    B3: {"boxes":1485235,"refuted":742618,"open":0,...}
    D:  25 pieces, 2,256,565 boxes in all, each with "open":0; {"c":0.1,"theta":0.95,"K":[2,5,7,8,9],"closed":true,...}

These are the counts in the proof of Proposition 6.2. On one core B1 and B2 take seconds, B3 and D a few
minutes each; `certificates.log` has the full output.

    node ntest.js
    node v6test.js 0.1 0.95

    {"checks":144000,"bad":0}
    {"found":1093,"valid":896,"inD":887,"inPieces":92,"refuted":0}

`ntest.js` must report no bad enclosure, and `v6test.js` must refute none of the boxes.
