#!/usr/bin/env python3
"""Can the one arranged coincidence be DRAWN instead of solved?

The peeling certificate says one coincidence is unavoidable, so no reordering
will ever remove it - that was a false trail, and this file no longer looks for
one. What reordering cannot do, a construction still can: Chiodo does not
remove his coincidence either, he determines the parameter that makes it hold
by intersecting two things he can draw.

Written out, that is a sharp condition. The fit is one incidence too many at
some object. Everything that does not move with the parameter t is already
drawn. So:

    the coincidence is a condition on ONE moving point Q,
    of the form  Q lies on C*,  where C* is built from things that stand still.

C* is therefore constructible. If Q's own locus is constructible too, then

    Q  =  C*  meeting that locus

and t is read off from Q instead of solved for. Q's locus is constructible
when Q is a seed - on the axis, on the rim, or on a line through its fixed
coordinate - and ALSO when Q is a node cut from the axis or the rim, a vertex
on the ring circle, or a crossing of two circles one of which stands still,
since then it simply slides along that circle. The first version of this file
demanded a seed and threw the others away, which was wrong.

Two shapes qualify, and the earlier --d0 sweep asked about only the first, with
a positional proxy for standing still rather than a measurement:

  fit at a circle C, extra point Q
      C is drawn through need(C) points, all still; the coincidence is Q on C.
      C* = C.  Q = C meeting its seed locus.

  fit at a node n confined to a locus L, circles C1 and C2
      the coincidence is that C1, C2 and L are concurrent. If C1 stands still,
      X = C1 meeting L is constructible, and C2 must pass through X. With all
      of C2's other determining points still,
      C* = the circle through those points and X.  Q = C* meeting its locus.

Standing still is measured - move the parameter, run twice, look - not read off
a dependence graph. Every hit is verified by rebuilding C* and checking that
the solved figure really does put Q on it.

    python3 concurrent.py ring1 orders_ring1_f1.json
    python3 concurrent.py ring2 orders_ring2_f1.json
"""
import collections, json, sys
from mpmath import mp, mpf

import dcheck as K
import dscan as Z
import flipexec as Fx
import score_orders as S

mp.dps = 40
RE = lambda z: mpf(z.real) if hasattr(z, 'real') else mpf(z)
STILL = mpf(10) ** -20


def norm(g):
    m = max(abs(RE(x)) for x in g)
    return [RE(x) / m for x in g] if m else None


def snapshots(rb, base, R, slot, h=mpf(1) / 600):
    v = list(base)
    v[slot] = base[slot] + h
    return Fx.run(rb, base, R), Fx.run(rb, v, R)


def circle_still(A, B, i):
    a, b = norm(A[1][i]), norm(B[1][i])
    return a is not None and b is not None and max(abs(x - y) for x, y in zip(a, b)) < STILL


def point_still(A, B, ref):
    a, b = ((A[2][ref[0]], B[2][ref[0]]) if isinstance(ref, list) else (A[0][ref], B[0][ref]))
    if a is None or b is None:
        return False
    return abs(RE(a[0]) - RE(b[0])) + abs(RE(a[1]) - RE(b[1])) < STILL


def step_of(rb, kind, key):
    for st in rb['steps']:
        if st[0] in ('s2', 's3') and kind == 's' and st[1] == key:
            return st
        if st[0] in ('na', 'nr', 'nx') and kind == 'n' and st[1] == key:
            return st
        if st[0] in ('t', 'ta') and kind == 't' and st[1] == key:
            return st
    return None


LOCUS = {'na': 'the axis', 'nr': 'the rim', 't': 'the ring circle', 'ta': 'the axis'}


def q_locus(rb, slot, ref, A, B):
    """Q is drawable if the curve it runs along as the parameter moves is one
    we can draw. A seed runs along the axis, the rim, or a line. But so does a
    node cut from the axis or the rim, a ring vertex on the ring circle, and -
    the case the first version of this file wrongly threw away - a crossing of
    two circles one of which stands still: it slides along that circle."""
    if isinstance(ref, list):
        return 'the ring circle'
    L = seed_locus(rb, slot, ref)
    if L:
        return L
    st = step_of(rb, 'n', ref)
    if st is None:
        return None
    if st[0] == 'na':
        return 'the axis'
    if st[0] == 'nr':
        return 'the rim'
    if st[0] == 'nx':
        for c in (st[2], st[3]):
            if circle_still(A, B, c):
                return 'circle %d, which stands still' % c
    return None


def seed_locus(rb, slot, ref):
    """Q is drawable only if it is the seed carrying this very parameter."""
    if isinstance(ref, list):
        return None
    k = 0
    for s in rb['seeds']:
        if s['kind'] == 'pole':
            continue
        if s['kind'] == 'free':
            if s['i'] == ref and slot in (k, k + 1):
                return ('a line through its fixed coordinate'
                        if slot == k else 'a line through its fixed coordinate')
            k += 2
        else:
            if s['i'] == ref and slot == k:
                return 'the axis' if s['kind'] == 'axis' else 'the rim'
            k += 1
    return None


def examine(rb, slot, base, R, why):
    f = rb['fit']
    A, B = snapshots(rb, base, R, slot)
    if A is None or B is None:
        return None

    if f['kind'] == 's':
        C = f['circle']
        Q = f['point']
        why['fit at a circle'] += 1
        if not circle_still(A, B, C):
            why['  but the fit circle moves'] += 1
            return None
        loc = q_locus(rb, slot, Q, A, B)
        if loc is None:
            why['  circle stands still, but the extra point runs on no drawable curve'] += 1
            return None
        return dict(shape='circle', C=C, Q=Q, locus=loc, through=None, node=None)

    if f['kind'] == 't':
        j = f['tip']
        stn = step_of(rb, 't', j)
        if stn is None:
            return None
        if stn[0] == 'ta':
            # the vertex where the ring circle crosses the axis is (0, R): given
            # outright, nothing to stand still. The coincidence is just "this
            # circle passes through the top of the ring"
            why['fit at the axis vertex (0, R)'] += 1
            anchors, Cf, X = [], f['circle'], ('t', j)
        else:
            why['fit at a ring vertex'] += 1
            anchors, Cf, X = [stn[2]], f['circle'], ('t', j)
    else:
        n = f['node']
        stn = step_of(rb, 'n', n)
        if stn is None:
            return None
        X = ('n', n)
        if stn[0] in ('na', 'nr'):
            why['fit at a node on a locus'] += 1
            anchors, Cf = [stn[2]], f['circle']
        else:
            why['fit at a plain crossing'] += 1
            anchors, Cf = [stn[2], stn[3]], f['circle']
    # the point X is drawable when everything that makes it stands still: one
    # circle plus a given locus, or two circles meeting each other
    pairs = []
    if all(circle_still(A, B, c) for c in anchors):
        pairs.append((anchors, Cf))
    if len(anchors) == 1 and circle_still(A, B, Cf):
        pairs.append(([Cf], anchors[0]))
    for C1s, C2 in pairs:
        stm = step_of(rb, 's', C2)
        if stm is None:
            continue
        pts = list(stm[2:4]) if stm[0] == 's2' else list(stm[2:5])
        moving = [p for p in pts if not point_still(A, B, p)]
        if len(moving) != 1:
            why['  one circle still, but %d of %d points move' % (len(moving), len(pts))] += 1
            continue
        loc = q_locus(rb, slot, moving[0], A, B)
        if loc is None:
            why['  one circle and all but one point still, but that point runs on no drawable curve'] += 1
            continue
        return dict(shape=X[0], C=C2, Q=moving[0], locus=loc, node=X[1],
                    still=C1s, L=LOCUS.get(stn[0], 'each other'),
                    through=[p for p in pts if p != moving[0]])
    why['  what makes the fit point does not stand still' if anchors
        else '  the fit circle moves through more than one point'] += 1
    return None


def verify(rb, base, R, slot, plan):
    """solve the fit for real, then check the solved figure puts Q on C*"""
    t = S.solve(rb, [float(x) for x in base], float(R), slot, span=0.06)
    if t is None:
        return None
    v = list(base)
    v[slot] = mpf(t)
    r = Fx.run(rb, v, R)
    if r is None:
        return None
    P, G, T = r[0], r[1], r[2]      # run() also returns the residual
    pt = lambda q: (T[q[0]] if isinstance(q, list) else P[q])
    if plan['shape'] == 'circle':
        g = G[plan['C']]
    else:
        anchor = T[plan['node']] if plan['shape'] == 't' else P[plan['node']]
        pts = [pt(q) for q in plan['through']] + [anchor]
        if any(p is None for p in pts):
            return None
        pts = [(RE(p[0]), RE(p[1])) for p in pts]
        g = (Fx.through_axis2(*pts) if len(pts) == 2 else Fx.through3(*pts))
    q = pt(plan['Q'])
    if g is None or q is None:
        return None
    x, y = RE(q[0]), RE(q[1])
    g = [RE(z) for z in g]
    m = max(abs(z) for z in g)
    return abs(g[0]*(x*x + y*y) + g[1]*x + g[2]*y + g[3]) / m


def main():
    which, path = sys.argv[1], sys.argv[2]
    why = collections.Counter()
    hits, scanned = [], 0
    for oi, od in enumerate(json.load(open(path))):
        rb = S.compile_one(which, [tuple(o) for o in od], Z.tips_of(which))
        if rb is None:
            continue
        base, R = K.seeds_of(rb), mpf(rb['R'])
        back = K.backward_cone(rb)
        for slot in range(len(base)):
            if K.depth_measured(rb, slot, back, base, R) is None:
                continue
            scanned += 1
            plan = examine(rb, slot, base, R, why)
            if not plan:
                continue
            err = verify(rb, base, R, slot, plan)
            plan.update(order=oi, slot=slot, err=None if err is None else mp.nstr(err, 4))
            hits.append(plan)
            print('  HIT order %d, parameter %d' % (oi, slot), flush=True)
    print('\n%s: %d (order, parameter) pairs examined, %d drawable' % (which, scanned, len(hits)))
    print('where it stops:')
    for k, v in sorted(why.items()):
        print('   %-62s %d' % (k, v))
    for h in hits:
        print('\n  order %d, parameter %d' % (h['order'], h['slot']))
        if h['shape'] == 'circle':
            print('    circle %d stands still; the coincidence is that node %s lies on it'
                  % (h['C'], h['Q']))
        else:
            print('    %s and circle %d meet %s at %s %d; %s stand(s) still'
                  % (', '.join('circle %d' % c for c in h['still']), h['C'], h['L'],
                     'ring vertex' if h['shape'] == 't' else 'node', h['node'],
                     ', '.join('circle %d' % c for c in h['still'])))
            print('    so node %d is drawable, and circle* runs through %s and it'
                  % (h['node'], ', '.join(str(x) for x in h['through'])))
        print('    the parameter is then node %s = circle* meeting %s - drawn, not solved'
              % (h['Q'], h['locus']))
        print('    check: the solved figure puts it off circle* by %s' % h['err'])
    if hits:
        json.dump(hits, open('drawable_%s.json' % which, 'w'), indent=1)
        print('\nwritten drawable_%s.json' % which)


if __name__ == '__main__':
    main()
