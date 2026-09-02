#!/usr/bin/env python3
"""concurrent.py over many collections at once, with the redundant work removed.

The slow part was running the construction four times per (order, parameter):
twice to see whether the residual notices the parameter at all, and twice more
to see what stands still. One pair of runs answers both, so this does half the
work and prints a line per file.

    python3 scanall.py ring1 orders_ring1_f1_s0.json orders_ring1_f1_t4.json ...
"""
import collections, json, sys, time
from mpmath import mp, mpf

import concurrent as CC
import dcheck as K
import dscan as Z
import flipexec as Fx
import score_orders as S

mp.dps = 40


def scan(which, path):
    why = collections.Counter()
    hits, scanned, orders_ok = [], 0, 0
    for oi, od in enumerate(json.load(open(path))):
        rb = S.compile_one(which, [tuple(o) for o in od], Z.tips_of(which))
        if rb is None:
            continue
        orders_ok += 1
        base, R = K.seeds_of(rb), mpf(rb['R'])
        for slot in range(len(base)):
            A, B = CC.snapshots(rb, base, R, slot)
            if A is None or B is None:
                continue
            if abs(CC.RE(A[3]) - CC.RE(B[3])) < mpf(10) ** -25:
                continue                      # the residual never sees this one
            scanned += 1
            plan = CC.examine(rb, slot, base, R, why)
            if not plan:
                continue
            err = CC.verify(rb, base, R, slot, plan)
            plan.update(order=oi, slot=slot, err=None if err is None else mp.nstr(err, 4))
            hits.append(plan)
    return scanned, orders_ok, hits, why


if __name__ == '__main__':
    which = sys.argv[1]
    allhits = []
    for path in sys.argv[2:]:
        t0 = time.time()
        scanned, ok, hits, why = scan(which, path)
        top = sorted(why.items(), key=lambda kv: -kv[1])[:3]
        print('%-30s %4d orders compiled, %5d pairs, %3d DRAWABLE   (%.0f s)'
              % (path.replace('orders_', '').replace('.json', ''), ok, scanned,
                 len(hits), time.time() - t0), flush=True)
        for k, v in top:
            print('     %-64s %d' % (k.strip(), v))
        for h in hits[:3]:
            print('     hit: order %d slot %d, %s, off by %s'
                  % (h['order'], h['slot'], h['locus'], h['err']))
        allhits += [dict(h, file=path) for h in hits]
    if allhits:
        out = 'drawable_all_%s.json' % which
        json.dump(allhits, open(out, 'w'), indent=1)
        print('\n%d hits written to %s' % (len(allhits), out))
