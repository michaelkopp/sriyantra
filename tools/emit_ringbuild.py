#!/usr/bin/env python3
"""Compile the two ring orders into explicit recipes the viewer can run.

Opcodes, in the spirit of ARCREC.steps:

    ["s2", s, p, q]           circle s, centred on the axis, through p and q
    ["s3", s, p, q, r]        circle s through three points
    ["nx", n, s1, s2, hx, hy] node n where circles s1 and s2 cross
    ["na", n, s, hy]          node n where circle s meets the axis
    ["nr", n, s, hx, hy]      node n where circle s meets the rim
    ["t",  k, s, hx, hy]      ring tip k where circle s meets E
    ["ta", k, hy]             ring tip k where E meets the axis

A point argument is a number for a node, or [k] for ring tip k.
"""
import json, math, re, sys
import ring_construct as C


def emit(which, orderfile, cubesfile, solve_idx):
    F = C.setup(which)
    order = json.load(open(orderfile))
    tips = json.loads(re.search(r'DATA = json.loads\(r"""(.*?)"""\)',
                                open(cubesfile).read(), re.S).group(1))['tips']
    steps, tipkey = C.compile_order(order, tips)
    C.FITSTEP[0] = C.fit_step(steps)

    # flatten the tip instances: one entry per actual ring vertex
    TIPS, tipidx = [], {}
    for orb, inst in sorted(F['TIP'].items()):
        for p, ss in inst:
            tipidx[(orb, p)] = len(TIPS)
            TIPS.append(dict(orb=list(orb), sides=list(ss), h=[p[0], p[1]]))

    out, seeds = [], []
    P, Gc, T = {}, {}, {}
    vals = iter(C.base_seeds(F, steps))
    fit = None
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
                    out.append(['na', n, cs[0], hy])
                    P[n] = C.meet_axis(Gc[cs[0]], (hx, hy))
                elif C.RIM[key]:
                    out.append(['nr', n, cs[0], hx, hy])
                    P[n] = C.meet_r(Gc[cs[0]], 1.0, (hx, hy))
                else:
                    out.append(['nx', n, cs[0], cs[1], hx, hy])
                    P[n] = C.meet(Gc[cs[0]], Gc[cs[1]], (hx, hy))
                if st is C.FITSTEP[0] and n == C.NORB[key][0]:
                    fit = dict(kind='n', node=n, circle=cs[nd])
                break                                   # the mirror is automatic
            for n in C.NORB[key]:
                if n not in P:
                    P[n] = (-P[C.MIR[n]][0], P[C.MIR[n]][1])
        elif kind == 't':
            for p0, ss in F['TIP'][tipkey[key]]:
                j = tipidx[(tipkey[key], p0)]
                src = [s for s in ss if s in Gc]
                if src:
                    out.append(['t', j, src[0], p0[0], p0[1]])
                    T[j] = C.meet_r(Gc[src[0]], F['R'], p0)
                else:
                    out.append(['ta', j, p0[1]])
                    T[j] = C.meet_axis((1.0, 0.0, 0.0, -F['R'] ** 2), p0)
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
                if nd == 2:
                    out.append(['s2', s, refs[0], refs[1]])
                    Gc[s] = C.through_axis2(pts[0], pts[1])
                else:
                    out.append(['s3', s, refs[0], refs[1], refs[2]])
                    Gc[s] = C.through3(pts[0], pts[1], pts[2])
                if st is C.FITSTEP[0] and s == C.SORB[key][0]:
                    fit = dict(kind='s', circle=s, point=refs[nd])
                break                                   # the mirror is automatic
            for s in C.SORB[key]:
                if s not in Gc:
                    g = Gc[C.SMIR[s]]
                    Gc[s] = [g[0], -g[1], g[2], g[3]]
    return dict(R=F['R'], seeds=seeds, tips=TIPS, steps=out, fit=fit,
                solveIdx=solve_idx, mir=C.MIR, smir=C.SMIR)


if __name__ == '__main__':
    # the best-conditioned of two hundred, by measured slider reach
    B = {'ring1': emit('ring1', '../order_ring1_best.json', '../ring1_cubes.py', 0),
         'ring2': emit('ring2', '../order_ring2_best.json', '../ring2_cubes.py', 9)}
    for k, v in B.items():
        npar = sum(2 if s['kind'] == 'free' else 0 if s['kind'] == 'pole' else 1
                   for s in v['seeds'])
        print('%s: %d steps, %d tips, %d seed parameters, fit %s'
              % (k, len(v['steps']), len(v['tips']), npar, v['fit']))
    json.dump(B, open('ringbuild.json', 'w'), separators=(',', ':'))
    print('written ringbuild.json (%d bytes)' % len(open('ringbuild.json').read()))
