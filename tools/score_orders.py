#!/usr/bin/env python3
"""Score the collected orders: they are the same mathematics, so what is left
to choose between them is how the arithmetic behaves.

Stage 1, on all of them: amplification. Nudge each free parameter and see how
far the figure moves. A chain of sixty circle-through-three-points multiplies a
change in a seed by a factor that varies wildly between orders, and that factor
is what decides whether a slider does anything before the construction leaves
the region where its steps are defined.

Stage 2, on the best few: the actual reach of every slider, measured by pushing
outward until the arranged coincidence can no longer be solved.
"""
import json, math, sys, time
import ring_construct as C

TIPS_META = None
BUILD = None


def compile_one(which, order, tips_meta):
    """the order -> an executable recipe, in the format the viewer runs"""
    import emit_ringbuild as E
    F = C.setup(which)
    steps, tipkey = C.compile_order(order, tips_meta)
    C.FITSTEP[0] = C.fit_step(steps)
    TIPS, tipidx = [], {}
    for orb, inst in sorted(F['TIP'].items()):
        for p, ss in inst:
            tipidx[(orb, p)] = len(TIPS)
            TIPS.append(dict(orb=list(orb), sides=list(ss), h=[p[0], p[1]]))
    out, seeds, fit = [], [], None
    P, Gc, T = {}, {}, {}
    for st in steps:
        kind, key = st['obj']
        bef, nd = st['before'], C.need_of(*st['obj'])
        if kind == 'n' and not bef:
            x, y = F['P'][key]
            k = ('pole' if (C.AXIS[key] and C.RIM[key]) else 'axis' if C.AXIS[key]
                 else 'rim' if C.RIM[key] else 'free')
            seeds.append(dict(i=key, kind=k, x=x, y=y,
                              a=(math.atan2(y, x) if k == 'rim' else None)))
            P[key] = (x, y)
            if C.MIR[key] != key:
                P[C.MIR[key]] = (-x, y)
            continue
        if kind == 'n':
            for n in C.NORB[key]:
                cs = sorted({s for o in bef if o[0] == 's' for s in C.SORB[o[1]]
                             if s in C.SIDES_OF[n] and s in Gc})
                hx, hy = F['P'][n]
                if C.AXIS[key] and C.RIM[key]:
                    continue
                if C.AXIS[key]:
                    out.append(['na', n, cs[0], hy]); P[n] = C.meet_axis(Gc[cs[0]], (hx, hy))
                elif C.RIM[key]:
                    out.append(['nr', n, cs[0], hx, hy]); P[n] = C.meet_r(Gc[cs[0]], 1.0, (hx, hy))
                else:
                    out.append(['nx', n, cs[0], cs[1], hx, hy])
                    P[n] = C.meet(Gc[cs[0]], Gc[cs[1]], (hx, hy))
                if P[n] is None:
                    return None
                if st is C.FITSTEP[0] and n == C.NORB[key][0]:
                    fit = dict(kind='n', node=n, circle=cs[nd])
                break
            for n in C.NORB[key]:
                if n not in P:
                    P[n] = (-P[C.MIR[n]][0], P[C.MIR[n]][1])
        elif kind == 't':
            for p0, ss in F['TIP'][tipkey[key]]:
                j = tipidx[(tipkey[key], p0)]
                src = [s for s in ss if s in Gc]
                if src and nd > 0:
                    out.append(['t', j, src[0], p0[0], p0[1]])
                    T[j] = C.meet_r(Gc[src[0]], F['R'], p0)
                else:
                    # the axis vertex is (0, R) outright, whether or not its
                    # circle is drawn yet; producing it as circle-meets-E made
                    # the fit there a tautology, on(circle, circle meets E) = 0
                    out.append(['ta', j, p0[1]])
                    T[j] = C.meet_axis((1.0, 0.0, 0.0, -F['R'] ** 2), p0)
                if T[j] is None:
                    return None
                # a ring vertex can carry the coincidence too - it wants one
                # circle and has two - and this branch used to drop it, so any
                # order with the fit pinned to a tip compiled to nothing
                if st is C.FITSTEP[0] and fit is None and len(src) > nd:
                    fit = dict(kind='t', tip=j, circle=src[nd])
        else:
            for s in C.SORB[key]:
                pts, refs = [], []
                for o in bef:
                    if o[0] == 'n':
                        c = [n for n in C.NORB[o[1]] if n in C.SIDES[s]]
                        if c:
                            pts.append(P[c[0]]); refs.append(c[0])
                    else:
                        for (orb, p0), j in tipidx.items():
                            if orb == tipkey[o[1]] and s in TIPS[j]['sides']:
                                pts.append(T[j]); refs.append([j])
                if len(refs) < nd:
                    return None
                if nd == 2:
                    out.append(['s2', s, refs[0], refs[1]])
                    Gc[s] = C.through_axis2(pts[0], pts[1])
                else:
                    out.append(['s3', s, refs[0], refs[1], refs[2]])
                    Gc[s] = C.through3(pts[0], pts[1], pts[2])
                if st is C.FITSTEP[0] and s == C.SORB[key][0]:
                    fit = dict(kind='s', circle=s, point=refs[nd])
                break
            for s in C.SORB[key]:
                if s not in Gc:
                    g = Gc[C.SMIR[s]]
                    Gc[s] = [g[0], -g[1], g[2], g[3]]
    if fit is None:
        return None
    return dict(R=F['R'], seeds=seeds, tips=TIPS, steps=out, fit=fit,
                mir=C.MIR, smir=C.SMIR)


def run(rb, vals, R):
    P = [None] * 69
    G = [None] * 27
    T = [None] * len(rb['tips'])
    mir, smir = rb['mir'], rb['smir']

    def put(i, q):
        P[i] = q
        if mir[i] != i:
            P[mir[i]] = (-q[0], q[1])
    k = 0
    for s in rb['seeds']:
        if s['kind'] == 'pole':
            put(s['i'], (s['x'], s['y']))
        elif s['kind'] == 'axis':
            put(s['i'], (0.0, vals[k])); k += 1
        elif s['kind'] == 'rim':
            put(s['i'], (math.cos(vals[k]), math.sin(vals[k]))); k += 1
        else:
            put(s['i'], (vals[k], vals[k + 1])); k += 2
    pt = lambda r: (T[r[0]] if isinstance(r, list) else P[r])
    for st in rb['steps']:
        op = st[0]
        if op in ('s2', 's3'):
            g = (C.through_axis2(pt(st[2]), pt(st[3])) if op == 's2'
                 else C.through3(pt(st[2]), pt(st[3]), pt(st[4])))
            G[st[1]] = g
            m = smir[st[1]]
            if G[m] is None:
                G[m] = g if m == st[1] else [g[0], -g[1], g[2], g[3]]
        else:
            if op == 'na':
                q = C.meet_axis(G[st[2]], (0.0, st[3]))
            elif op == 'nr':
                q = C.meet_r(G[st[2]], 1.0, (st[3], st[4]))
            elif op == 'nx':
                q = C.meet(G[st[2]], G[st[3]], (st[4], st[5]))
            elif op == 't':
                q = C.meet_r(G[st[2]], R, (st[3], st[4]))
            else:
                q = (0.0, R if st[2] > 0 else -R)
            if q is None:
                return None
            if op in ('t', 'ta'):
                T[st[1]] = q
            else:
                put(st[1], q)
    f = rb['fit']
    q = (T[f['tip']] if f['kind'] == 't' else
         pt(f['point']) if f['kind'] == 's' else P[f['node']])
    g = G[f['circle']]
    if g is None or q is None:
        return None
    return P, G, T, C.on(g, q)


def seeds_of(rb):
    v = []
    for s in rb['seeds']:
        if s['kind'] == 'free':
            v += [s['x'], s['y']]
        elif s['kind'] == 'rim':
            v.append(s['a'])
        elif s['kind'] == 'axis':
            v.append(s['y'])
    return v


def solve(rb, vals, R, idx, span=0.5):
    def f(t):
        w = list(vals); w[idx] = t
        r = run(rb, w, R)
        return None if r is None else r[3]
    t0 = vals[idx]
    xs, ys = [], []
    d = 1e-9
    while d <= span:
        for sgn in (1, -1):
            y = f(t0 + sgn * d)
            if y is not None:
                xs.append(t0 + sgn * d); ys.append(y)
        d *= 1.7
    y0 = f(t0)
    if y0 is not None:
        xs.append(t0); ys.append(y0)
    o = sorted(range(len(xs)), key=lambda i: xs[i])
    xs = [xs[i] for i in o]; ys = [ys[i] for i in o]
    best = None
    for i in range(1, len(xs)):
        if ys[i - 1] * ys[i] < 0 and (best is None or abs(xs[i] - t0) < abs(best[1] - t0)):
            best = (xs[i - 1], xs[i], ys[i - 1], ys[i])
    if best is None:
        return None
    a, b, fa, fb = best
    for _ in range(120):
        m = (a + b) / 2
        fm = f(m)
        if fm is None:
            return None
        if fa * fm <= 0:
            b, fb = m, fm
        else:
            a, fa = m, fm
    return (a + b) / 2


def pick_solve_index(rb, vals, R):
    """the parameter the coincidence is most sensitive to"""
    base = run(rb, vals, R)
    if base is None:
        return None
    best, bi = 0.0, None
    for i in range(len(vals)):
        w = list(vals); w[i] += 1e-7
        r = run(rb, w, R)
        if r is None:
            continue
        s = abs(r[3] - base[3])
        if s > best:
            best, bi = s, i
    return bi


def amplification(rb, vals, R, idx, d=1e-7):
    """how far the figure moves per unit of parameter, with the fit re-solved"""
    t = solve(rb, vals, R, idx)
    if t is None:
        return None
    v = list(vals); v[idx] = t
    A = run(rb, v, R)
    if A is None:
        return None
    out = []
    for i in list(range(len(vals))) + ['R']:
        if i == idx:
            continue
        w = list(v); RR = R
        if i == 'R':
            RR = R + d
        else:
            w[i] += d
        tt = solve(rb, w, RR, idx)
        if tt is None:
            out.append(None); continue
        w[idx] = tt
        B = run(rb, w, RR)
        if B is None:
            out.append(None); continue
        out.append(max(math.hypot(A[0][n][0] - B[0][n][0], A[0][n][1] - B[0][n][1])
                       for n in range(69)) / d)
    return v, out
