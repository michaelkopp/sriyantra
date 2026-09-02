#!/usr/bin/env python3
"""d(t) over a whole collection of orders, measured rather than reasoned about.

Every order in the file is compiled, and for each of its free parameters the
depth is taken twice: once by walking the dependence graph (degree.py, fast
but the kind of code that can be quietly wrong) and once by moving the
parameter and seeing what actually moved. Only the second is trusted; the
first is printed beside it so a disagreement shows.

    python3 dscan.py ring1 orders_ring1_f1.json
"""
import json, os, re, sys, time
from multiprocessing import Pool

from mpmath import mp, mpf

import degree as D
import dcheck as K
import score_orders as S

TIPS = {}


def tips_of(which):
    if which not in TIPS:
        TIPS[which] = json.loads(re.search(
            r'DATA = json.loads\(r"""(.*?)"""\)',
            open('%s_cubes.py' % which).read(), re.S).group(1))['tips']
    return TIPS[which]


def one(job):
    which, idx, order = job
    mp.dps = 30
    try:
        rb = S.compile_one(which, order, tips_of(which))
    except Exception:
        return idx, None
    if rb is None:
        return idx, None
    back = K.backward_cone(rb)
    base = K.seeds_of(rb)
    R = mpf(rb['R'])
    spi = D.seed_param_index(rb)
    out = []
    for slot in range(len(base)):
        m = K.depth_measured(rb, slot, back, base, R)
        if m is None:
            continue
        g = D.depth_for(rb, spi[slot])
        out.append((slot, m[0], g))
    return idx, (out, rb['fit'])


def main():
    which, path = sys.argv[1], sys.argv[2]
    orders = json.load(open(path))
    jobs = [(which, i, [tuple(o) for o in od]) for i, od in enumerate(orders)]
    t0 = time.time()
    res = {}
    with Pool(int(os.environ.get('J', '6'))) as pool:
        for idx, r in pool.imap_unordered(one, jobs, chunksize=4):
            res[idx] = r
    ok = {i: r for i, r in res.items() if r}
    print('%s: %d orders, %d compiled, %.0f s'
          % (which, len(orders), len(ok), time.time() - t0))

    rows = []
    dis = 0
    for i, (rowset, fit) in sorted(ok.items()):
        if not rowset:
            continue
        best = min(r[1] for r in rowset)
        rows.append((best, i, [r for r in rowset if r[1] == best], fit))
        dis += sum(1 for r in rowset if r[1] != r[2])
    rows.sort()
    tot = sum(len(r[0]) for r in [([x for x in rs], ) for _, _, rs, _ in rows])
    print('graph disagrees with measurement on %d of %d (order, parameter) pairs'
          % (dis, sum(len(r[0]) for r in ok.values() if r)))
    print('\nsmallest measured d, best fifteen orders:')
    for best, i, at, fit in rows[:15]:
        print('   d = %-3d order %-4d slots %-16s fit %s'
              % (best, i, ','.join(str(a[0]) for a in at), fit))
    hist = {}
    for best, _, _, _ in rows:
        hist[best] = hist.get(best, 0) + 1
    print('\nhow many orders reach each minimum:')
    for d in sorted(hist):
        print('   d = %-3d %s (%d)' % (d, '#' * min(hist[d], 60), hist[d]))
    json.dump([[b, i, [a[0] for a in at]] for b, i, at, _ in rows],
              open('dscan_%s.json' % which, 'w'))


if __name__ == '__main__':
    main()
