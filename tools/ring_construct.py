#!/usr/bin/env python3
"""The Chiodo-type construction of a ring-held circular-arc Sri Yantra.

The order from fits.py is combinatorial. This executes it as geometry, and
solves the single arranged coincidence it leaves over. What comes out is a
construction: choose the seed parameters and the ring's radius, solve one
equation for one of them, and every incidence of the figure - and the
concyclicity of the whole grown ring - follows from the steps.

    ring 1   twelve seed parameters, E1's radius, one equation
    ring 2   eleven seed parameters, E2's radius, one equation
"""
import json, math, sys
import numpy as np

import os as _os
# recipe.json sits beside this file in the repository, and one directory up
# under paper/ in the working tree the scripts were developed in
for _p in (_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), 'recipe.json'),
           '../paper/recipe.json'):
    if _os.path.exists(_p):
        REC = json.load(open(_p)); break
else:
    raise SystemExit('recipe.json not found beside ring_construct.py')
SIDES = [s['nodes'] for s in REC['sides']]
MIR = REC['recipe']['mir']
INC114 = REC['inc']
NS, NN = len(SIDES), len(MIR)
_key = {frozenset(s): i for i, s in enumerate(SIDES)}
SMIR = [_key[frozenset(MIR[n] for n in SIDES[s])] for s in range(NS)]
NREP = [min(n, MIR[n]) for n in range(NN)]
SREP = [min(s, SMIR[s]) for s in range(NS)]
NORB, SORB = {}, {}
for n in range(NN):
    NORB.setdefault(NREP[n], []).append(n)
for s in range(NS):
    SORB.setdefault(SREP[s], []).append(s)
SIDES_OF = [[s for s in range(NS) if n in SIDES[s]] for n in range(NN)]
AXIS = [MIR[n] == n for n in range(NN)]
RIM = [n in set(REC['outer']) for n in range(NN)]


# ---------------------------------------------------------------- geometry
def unit(g):
    g = np.asarray(g, float)
    n = np.linalg.norm(g)
    g = g / n if n else g
    return g if g[np.argmax(abs(g))] > 0 else -g


def through3(p1, p2, p3):
    """closed form instead of an SVD - the scoring runs this a few million times"""
    (x1, y1), (x2, y2), (x3, y3) = p1, p2, p3
    r1, r2, r3 = x1*x1 + y1*y1, x2*x2 + y2*y2, x3*x3 + y3*y3
    d3 = lambda a1,b1,c1,a2,b2,c2,a3,b3,c3: (a1*(b2*c3-b3*c2) - b1*(a2*c3-a3*c2)
                                             + c1*(a2*b3-a3*b2))
    a =  d3(x1,y1,1, x2,y2,1, x3,y3,1)
    b = -d3(r1,y1,1, r2,y2,1, r3,y3,1)
    c =  d3(r1,x1,1, r2,x2,1, r3,x3,1)
    d = -d3(r1,x1,y1, r2,x2,y2, r3,x3,y3)
    return unit([a, b, c, d])


def through_axis2(p1, p2):
    """centred on the axis, so b = 0 and two points suffice"""
    (x1, y1), (x2, y2) = p1, p2
    r1, r2 = x1*x1 + y1*y1, x2*x2 + y2*y2
    return unit([y1 - y2, 0.0, r2 - r1, r1*y2 - r2*y1])


def on(g, p):
    x, y = p
    return g[0] * (x * x + y * y) + g[1] * x + g[2] * y + g[3]


def _roots(a, b, c):
    if abs(a) < 1e-15:
        return [] if b == 0 else [-c / b]
    d = b * b - 4 * a * c
    if d < 0:
        return []
    s = math.sqrt(d)
    return [(-b + s) / (2 * a), (-b - s) / (2 * a)]


def meet_line(g, L, hint):
    A, B, C = L
    a, b, c, d = g
    if abs(B) > abs(A):
        p, q = -A / B, -C / B
        pts = [(x, p * x + q) for x in
               _roots(a * (1 + p * p), b + c * p + 2 * a * p * q, a * q * q + c * q + d)]
    else:
        p, q = -B / A, -C / A
        pts = [(p * y + q, y) for y in
               _roots(a * (1 + p * p), c + b * p + 2 * a * p * q, a * q * q + b * q + d)]
    if not pts:
        return None
    return min(pts, key=lambda P: (P[0] - hint[0]) ** 2 + (P[1] - hint[1]) ** 2)


def meet(g1, g2, hint):
    L = (g1[1] * g2[0] - g2[1] * g1[0], g1[2] * g2[0] - g2[2] * g1[0],
         g1[3] * g2[0] - g2[3] * g1[0])
    return None if max(abs(x) for x in L) < 1e-14 else meet_line(g1, L, hint)


def meet_axis(g, hint):
    return meet_line(g, (1.0, 0.0, 0.0), hint)


def meet_r(g, r, hint):
    return meet(g, (1.0, 0.0, 0.0, -r * r), hint)


# --------------------------------------------------------------- the figures
def setup(which):
    """which = 'ring1' or 'ring2'; returns the reference data for that family"""
    key = 'arcs' if which == 'ring1' else 'arcs2'
    P = REC['defaults'][key]['P']
    G = [unit(g) for g in REC['defaults'][key]['G']]
    rings = REC['rings'][key]
    tips = rings[0 if which == 'ring1' else 1]
    R = float(np.mean([math.hypot(*t['tip']) for t in tips]))
    inst = {}
    for t in tips:
        p = tuple(t['tip'])
        ss = tuple(s for s in range(NS) if abs(on(G[s], p)) < 1e-7)
        inst.setdefault(tuple(sorted({SREP[s] for s in ss})), []).append((p, ss))
    return dict(which=which, P=P, G=G, R=R, TIP=inst)


TIPNEED = {}


def compile_order(order, tips_meta):
    tipkey = {t['i']: tuple(sorted({t['a'], t['b']})) for t in tips_meta}
    # the vertex where the ring circle crosses the axis lies on ONE circle
    # orbit and is determined outright: need 0, not 1. Without this the fit
    # at that vertex - "this circle must pass through the top of the ring" -
    # was invisible to fit_step, and every order pinned there failed to compile
    TIPNEED.clear()
    TIPNEED.update({i: (0 if len(ab) == 1 else 1) for i, ab in tipkey.items()})
    pos = {(r[0], r[1]): k for k, r in enumerate(order)}
    steps = []
    for k, r in enumerate(order):
        kind, key = r[0], r[1]
        if kind == 's':
            inc = ([('n', NREP[n]) for s in SORB[key] for n in SIDES[s]] +
                   [('t', i) for i, ab in tipkey.items() if key in ab])
        elif kind == 'n':
            inc = [('s', SREP[s]) for s in SIDES_OF[key]]
        else:
            inc = [('s', x) for x in tipkey[key]]
        steps.append(dict(obj=(kind, key), at=k,
                          before=[q for q in sorted(set(inc)) if pos[q] < k]))
    return steps, tipkey


def need_of(kind, key):
    if kind == 's':
        return 2 if SMIR[SORB[key][0]] == SORB[key][0] else 3
    if kind == 't':
        return TIPNEED.get(key, 1)
    return 0 if (AXIS[key] and RIM[key]) else 1 if (AXIS[key] or RIM[key]) else 2


def fit_step(steps):
    f = [s for s in steps if s['before'] and len(s['before']) == need_of(*s['obj']) + 1]
    assert len(f) == 1, f
    return f[0]


def seed_spec(steps):
    out = []
    for st in steps:
        k, key = st['obj']
        if k == 'n' and not st['before']:
            out.append((key, 'pole' if (AXIS[key] and RIM[key]) else
                        'axis' if AXIS[key] else 'rim' if RIM[key] else 'free'))
    return out


def base_seeds(F, steps):
    v = []
    for key, kind in seed_spec(steps):
        x, y = F['P'][key]
        if kind == 'axis':
            v.append(y)
        elif kind == 'rim':
            v.append(math.atan2(y, x))
        elif kind == 'free':
            v += [x, y]
    return v


def run(F, steps, tipkey, vals, R, fit_as_extra=True):
    """execute; the arranged coincidence is left unimposed and reported"""
    P, Gc, T = {}, {}, {}
    fit = None
    it = iter(vals)
    for st in steps:
        kind, key = st['obj']
        bef = st['before']
        nd = need_of(kind, key)
        if kind == 'n':
            if not bef:
                if AXIS[key] and RIM[key]:
                    P[key] = tuple(F['P'][key])
                elif AXIS[key]:
                    P[key] = (0.0, next(it))
                elif RIM[key]:
                    th = next(it)
                    P[key] = (math.cos(th), math.sin(th))
                else:
                    P[key] = (next(it), next(it))
            else:
                for n in NORB[key]:
                    cs = sorted({s for o in bef if o[0] == 's' for s in SORB[o[1]]
                                 if s in SIDES_OF[n] and s in Gc})
                    if len(cs) < max(nd, 1):
                        return None
                    if AXIS[key] and RIM[key]:
                        P[n] = tuple(F['P'][n])
                    elif AXIS[key]:
                        P[n] = meet_axis(Gc[cs[0]], F['P'][n])
                    elif RIM[key]:
                        P[n] = meet_r(Gc[cs[0]], 1.0, F['P'][n])
                    else:
                        P[n] = meet(Gc[cs[0]], Gc[cs[1]], F['P'][n])
                    if P[n] is None:
                        return None
                    if st is FITSTEP[0] and n == NORB[key][0]:
                        fit = on(Gc[cs[nd]], P[n])
            for n in NORB[key]:
                if n not in P:
                    P[n] = (-P[MIR[n]][0], P[MIR[n]][1])
        elif kind == 't':
            for p0, ss in F['TIP'][tipkey[key]]:
                src = [s for s in ss if s in Gc]
                q = (meet_axis((1.0, 0.0, 0.0, -R * R), p0) if not src
                     else meet_r(Gc[src[0]], R, p0))
                if q is None:
                    return None
                T[(key, p0)] = q
        else:
            for s in SORB[key]:
                pts = []
                for o in bef:
                    if o[0] == 'n':
                        c = [n for n in NORB[o[1]] if n in SIDES[s]]
                        if c:
                            pts.append(P[c[0]])
                    else:
                        for (ti, p0), q in T.items():
                            if ti == o[1] and s in [x for x in F['TIP'][tipkey[ti]]
                                                    if x[0] == p0][0][1]:
                                pts.append(q)
                if len(pts) < nd:
                    return None
                Gc[s] = through_axis2(*pts[:nd]) if nd == 2 else through3(*pts[:nd])
                if st is FITSTEP[0] and s == SORB[key][0] and len(pts) > nd:
                    fit = on(Gc[s], pts[nd])
    return P, Gc, T, fit


FITSTEP = [None]


def solve_fit(F, steps, tipkey, vals, R, idx, span=1e-4, iters=200):
    """bisect the arranged coincidence in seed parameter `idx`"""
    def f(t):
        r = run(F, steps, tipkey, [*vals[:idx], t, *vals[idx + 1:]], R)
        return None if r is None else r[3]
    t0 = vals[idx]
    a, b = t0 - span, t0 + span
    fa, fb = f(a), f(b)
    while (fa is None or fb is None or fa * fb > 0) and span > 1e-14:
        span /= 2
        a, b = t0 - span, t0 + span
        fa, fb = f(a), f(b)
    if fa is None or fb is None or fa * fb > 0:
        return None
    for _ in range(iters):
        m = (a + b) / 2
        fm = f(m)
        if fm is None:
            return None
        if fa * fm <= 0:
            b, fb = m, fm
        else:
            a, fa = m, fm
    return (a + b) / 2


def report(which, orderfile, cubesfile, idx=0, dR=0.0):
    import re
    F = setup(which)
    order = json.load(open(orderfile))
    tips = json.loads(re.search(r'DATA = json.loads\(r"""(.*?)"""\)',
                                open(cubesfile).read(), re.S).group(1))['tips']
    steps, tipkey = compile_order(order, tips)
    FITSTEP[0] = fit_step(steps)
    v = base_seeds(F, steps)
    R = F['R'] + dR
    print('%s: %d objects, %d seed parameters, ring radius %.12f'
          % (which, len(steps), len(v), R))
    k, key = FITSTEP[0]['obj']
    print('  the one arranged coincidence: step %d, %s %d with %s'
          % (FITSTEP[0]['at'], {'s': 'circle', 'n': 'node'}[k], key, FITSTEP[0]['before']))
    r0 = run(F, steps, tipkey, v, R)
    print('  residual before solving: %.3e' % r0[3])
    t = solve_fit(F, steps, tipkey, v, R, idx)
    if t is None:
        print('  could not bracket a root'); return None
    v[idx] = t
    P, Gc, T, fit = run(F, steps, tipkey, v, R)
    rad = [math.hypot(*q) for q in T.values()]
    worst = max(abs(on(unit(Gc[s]), P[n])) for n, s in INC114)
    print('  solved parameter %d -> %.15f   residual %.2e' % (idx, t, fit))
    print('  ring: %d tips, radius %.14f, spread %.2e' %
          (len(rad), sum(rad) / len(rad), (max(rad) - min(rad)) / (sum(rad) / len(rad))))
    print('  worst of the 114 incidences: %.2e' % worst)
    return dict(seeds=v, R=R, P=P, fit=fit)


if __name__ == '__main__':
    report('ring1', '../order_ring1_f1.json', '../ring1_cubes.py', idx=0)
    print()
    report('ring2', '../order_ring2_f1.json', '../ring2_cubes.py', idx=0)
