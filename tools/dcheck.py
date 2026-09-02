#!/usr/bin/env python3
"""d(t), measured rather than reasoned about - and the two compared.

degree.py answers the question by walking the dependence graph, which is fast
enough to run over hundreds of orders. It is also the kind of code that can be
quietly wrong: its first version closed the dependence set under the figure's
mirror only after the forward pass, and so missed every square root that
becomes t-dependent through a mirrored object. It reported d = 1 where the
truth is 7.

The measured version cannot make that mistake. Move the parameter by 1e-12,
run the construction twice, and see which nodes moved. Anything that moves
depends on t; anything in the residual's backward cone feeds the equation. The
square roots in both sets are exactly the ones that stand between the two.

    python3 dcheck.py <order.json> ring1|ring2      one order, every parameter
    python3 dcheck.py --all                         every order file to hand
"""
import glob, json, os, re, sys
from mpmath import mp, mpf

import degree as D
import flipexec as F
import ring_construct as C
import score_orders as S

mp.dps = 40
SQRT_OPS = ('na', 'nr', 'nx', 't')


def compile_order(which, path):
    order = json.load(open(path))
    tips = json.loads(re.search(r'DATA = json.loads\(r"""(.*?)"""\)',
                                open('%s_cubes.py' % which).read(), re.S).group(1))['tips']
    return S.compile_one(which, order, tips)


def backward_cone(rb):
    f, mir, smir = rb['fit'], rb['mir'], rb['smir']
    ins_of = {}
    for st in rb['steps']:
        op = st[0]
        if op in ('s2', 's3'):
            a = st[2:4] if op == 's2' else st[2:5]
            ins_of[('s', st[1])] = [(('t', x[0]) if isinstance(x, list) else ('n', x))
                                    for x in a]
        elif op in ('na', 'nr'):
            ins_of[('n', st[1])] = [('s', st[2])]
        elif op == 'nx':
            ins_of[('n', st[1])] = [('s', st[2]), ('s', st[3])]
        else:
            ins_of[('t', st[1])] = [('s', st[2])] if op == 't' else []
    tgt = [('s', f['circle'])]
    tgt.append(('t', f['point'][0]) if f['kind'] == 's' else ('n', f['node']))
    back, stack = set(), list(tgt)
    while stack:
        o = stack.pop()
        if o in back:
            continue
        back.add(o)
        stack += ins_of.get(o, [])
        if o[0] == 'n' and mir[o[1]] != o[1]:
            stack.append(('n', mir[o[1]]))
        if o[0] == 's' and smir[o[1]] != o[1]:
            stack.append(('s', smir[o[1]]))
    return back


def seeds_of(rb):
    v = []
    for s in rb['seeds']:
        if s['kind'] == 'free':
            v += [mpf(s['x']), mpf(s['y'])]
        elif s['kind'] == 'rim':
            v.append(mpf(s['a']))
        elif s['kind'] == 'axis':
            v.append(mpf(s['y']))
    return v


def depth_measured(rb, slot, back, base, R, h=mpf('1e-12'), tol=mpf('1e-25')):
    """d(t) by finite differences: the ground truth"""
    def snap(t):
        v = list(base)
        v[slot] = t
        return F.run(rb, v, R)
    A = snap(base[slot])
    B = snap(base[slot] + h)
    if A is None or B is None:
        return None
    def moved(a, b):
        return a and b and abs(a[0] - b[0]) + abs(a[1] - b[1]) > tol
    mn = {i for i in range(69) if moved(A[0][i], B[0][i])}
    mt = {i for i in range(len(A[2])) if moved(A[2][i], B[2][i])}
    if not (('s', rb['fit']['circle']) in back):
        return None
    n, where = 0, []
    for j, st in enumerate(rb['steps']):
        if st[0] not in SQRT_OPS:
            continue
        o = ('t', st[1]) if st[0] == 't' else ('n', st[1])
        hit = (st[1] in mt) if st[0] == 't' else (st[1] in mn)
        if hit and o in back:
            n += 1
            where.append((j, st[0]))
    # does the residual see this parameter at all?
    if abs(A[3] - B[3]) < tol:
        return None
    return n, where


def report(rb, label):
    back = backward_cone(rb)
    base = seeds_of(rb)
    R = mpf(rb['R'])
    spi = D.seed_param_index(rb)
    rows = []
    for slot in range(len(base)):
        g = D.depth_for(rb, spi[slot])
        m = depth_measured(rb, slot, back, base, R)
        rows.append((slot, g, None if m is None else m[0]))
    bad = [r for r in rows if r[1] != r[2]]
    live = [r[2] for r in rows if r[2] is not None]
    print('%-28s  measured d: %s   min %s%s'
          % (label, ' '.join('%s' % ('-' if r[2] is None else r[2]) for r in rows),
             min(live) if live else '-',
             '' if not bad else '   GRAPH DISAGREES at slots %s'
             % ','.join(str(r[0]) for r in bad)))
    return rows


if __name__ == '__main__':
    if sys.argv[1:2] == ['--all']:
        for which in ('ring1', 'ring2'):
            for p in sorted(glob.glob('../order_%s_*.json' % which)):
                try:
                    rb = compile_order(which, p)
                except Exception as e:
                    print('%-28s  will not compile: %s' % (os.path.basename(p), e))
                    continue
                if rb is None:
                    print('%-28s  will not compile' % os.path.basename(p))
                    continue
                report(rb, os.path.basename(p))
    else:
        rb = compile_order(sys.argv[2], sys.argv[1])
        report(rb, os.path.basename(sys.argv[1]))
