#!/usr/bin/env python3
"""
ring2_cubes.py - the same question as ring1_cubes.py, for the SECOND grown ring.
("grown twice, outer ring held", but asked as a construction rather than a fit.)

WHAT IS DIFFERENT ABOUT RING 2
------------------------------
Ring 2's fourteen tips reduce to EIGHT independent orbits rather than seven, and
two of those lie on the axis - they are where a side crosses its own mirror
(U3.left x U3.right, and D1.right x D1.left).  That is a real structural
advantage: an axis tip is just E2 meeting the axis, so it is available before
any circle exists, and the side involved then has to be built through it.  The
other six behave like ring 1's: cut the tip out of E2 with one of its circles,
then build the other circle through it.

The counting: 8 orbit radii equal is 7 conditions on 18 parameters, so the
family is 11-dimensional - 10 seed parameters plus E2's radius.  (Ring 1 came
out at 12.)  Both numbers agree with what the solver is asked to find.

Everything else - the sequencing formulation, the cube-and-conquer split, the
buildability rules - is exactly as in ring1_cubes.py.

    pip install z3-solver
    python3 ring2_cubes.py             # all cores, level 4
    python3 ring2_cubes.py -j 8 -k 5
    python3 ring2_cubes.py --any       # drop the 10-seed pin
"""
import json, sys, time, argparse
from multiprocessing import Pool, cpu_count

DATA = json.loads(r"""
{
 "nodes": {
  "0": {
   "rep": 0,
   "axis": false,
   "rim": false,
   "need": 2,
   "seed": 2,
   "on": [
    12,
    13,
    25
   ]
  },
  "1": {
   "rep": 1,
   "axis": false,
   "rim": false,
   "need": 2,
   "seed": 2,
   "on": [
    13,
    22
   ]
  },
  "2": {
   "rep": 2,
   "axis": false,
   "rim": false,
   "need": 2,
   "seed": 2,
   "on": [
    13,
    24
   ]
  },
  "3": {
   "rep": 3,
   "axis": true,
   "rim": false,
   "need": 1,
   "seed": 1,
   "on": [
    13,
    15
   ]
  },
  "7": {
   "rep": 7,
   "axis": false,
   "rim": false,
   "need": 2,
   "seed": 2,
   "on": [
    1,
    12,
    22
   ]
  },
  "9": {
   "rep": 9,
   "axis": false,
   "rim": true,
   "need": 1,
   "seed": 1,
   "on": [
    0,
    1
   ]
  },
  "10": {
   "rep": 10,
   "axis": false,
   "rim": false,
   "need": 2,
   "seed": 2,
   "on": [
    0,
    16
   ]
  },
  "13": {
   "rep": 13,
   "axis": false,
   "rim": false,
   "need": 2,
   "seed": 2,
   "on": [
    1,
    6,
    16
   ]
  },
  "14": {
   "rep": 14,
   "axis": false,
   "rim": false,
   "need": 2,
   "seed": 2,
   "on": [
    1,
    3,
    19
   ]
  },
  "15": {
   "rep": 15,
   "axis": false,
   "rim": false,
   "need": 2,
   "seed": 2,
   "on": [
    1,
    9,
    25
   ]
  },
  "16": {
   "rep": 16,
   "axis": true,
   "rim": false,
   "need": 1,
   "seed": 1,
   "on": [
    1,
    24
   ]
  },
  "20": {
   "rep": 20,
   "axis": false,
   "rim": false,
   "need": 2,
   "seed": 2,
   "on": [
    7,
    24,
    25
   ]
  },
  "21": {
   "rep": 21,
   "axis": false,
   "rim": false,
   "need": 2,
   "seed": 2,
   "on": [
    10,
    22,
    24
   ]
  },
  "24": {
   "rep": 24,
   "axis": false,
   "rim": false,
   "need": 2,
   "seed": 2,
   "on": [
    9,
    10,
    19
   ]
  },
  "25": {
   "rep": 25,
   "axis": false,
   "rim": false,
   "need": 2,
   "seed": 2,
   "on": [
    10,
    25
   ]
  },
  "26": {
   "rep": 26,
   "axis": false,
   "rim": false,
   "need": 2,
   "seed": 2,
   "on": [
    10,
    15
   ]
  },
  "27": {
   "rep": 27,
   "axis": true,
   "rim": false,
   "need": 1,
   "seed": 1,
   "on": [
    10,
    18
   ]
  },
  "29": {
   "rep": 29,
   "axis": false,
   "rim": false,
   "need": 2,
   "seed": 2,
   "on": [
    9,
    22
   ]
  },
  "31": {
   "rep": 31,
   "axis": false,
   "rim": true,
   "need": 1,
   "seed": 1,
   "on": [
    6,
    7
   ]
  },
  "32": {
   "rep": 32,
   "axis": false,
   "rim": false,
   "need": 2,
   "seed": 2,
   "on": [
    3,
    7,
    16
   ]
  },
  "33": {
   "rep": 33,
   "axis": false,
   "rim": false,
   "need": 2,
   "seed": 2,
   "on": [
    7,
    19
   ]
  },
  "34": {
   "rep": 34,
   "axis": false,
   "rim": false,
   "need": 2,
   "seed": 2,
   "on": [
    7,
    15,
    22
   ]
  },
  "35": {
   "rep": 35,
   "axis": false,
   "rim": false,
   "need": 2,
   "seed": 2,
   "on": [
    7,
    18
   ]
  },
  "36": {
   "rep": 36,
   "axis": true,
   "rim": false,
   "need": 1,
   "seed": 1,
   "on": [
    7,
    21
   ]
  },
  "38": {
   "rep": 38,
   "axis": false,
   "rim": false,
   "need": 2,
   "seed": 2,
   "on": [
    6,
    19
   ]
  },
  "44": {
   "rep": 44,
   "axis": false,
   "rim": true,
   "need": 1,
   "seed": 1,
   "on": [
    21,
    22
   ]
  },
  "45": {
   "rep": 45,
   "axis": false,
   "rim": false,
   "need": 2,
   "seed": 2,
   "on": [
    4,
    18,
    22
   ]
  },
  "46": {
   "rep": 46,
   "axis": true,
   "rim": false,
   "need": 1,
   "seed": 1,
   "on": [
    3,
    22
   ]
  },
  "47": {
   "rep": 47,
   "axis": false,
   "rim": true,
   "need": 1,
   "seed": 1,
   "on": [
    15,
    16
   ]
  },
  "48": {
   "rep": 48,
   "axis": false,
   "rim": false,
   "need": 2,
   "seed": 2,
   "on": [
    4,
    15,
    19
   ]
  },
  "52": {
   "rep": 52,
   "axis": false,
   "rim": true,
   "need": 1,
   "seed": 1,
   "on": [
    3,
    4
   ]
  },
  "53": {
   "rep": 53,
   "axis": false,
   "rim": false,
   "need": 2,
   "seed": 2,
   "on": [
    3,
    25
   ]
  },
  "55": {
   "rep": 55,
   "axis": true,
   "rim": true,
   "need": 1,
   "seed": 1,
   "on": [
    4
   ]
  },
  "56": {
   "rep": 56,
   "axis": false,
   "rim": false,
   "need": 2,
   "seed": 2,
   "on": [
    4,
    21
   ]
  },
  "57": {
   "rep": 57,
   "axis": false,
   "rim": false,
   "need": 2,
   "seed": 2,
   "on": [
    4,
    16
   ]
  },
  "61": {
   "rep": 61,
   "axis": true,
   "rim": false,
   "need": 1,
   "seed": 1,
   "on": [
    6,
    25
   ]
  },
  "65": {
   "rep": 65,
   "axis": false,
   "rim": true,
   "need": 1,
   "seed": 1,
   "on": [
    18,
    19
   ]
  },
  "67": {
   "rep": 67,
   "axis": true,
   "rim": false,
   "need": 1,
   "seed": 1,
   "on": [
    0,
    19
   ]
  },
  "68": {
   "rep": 68,
   "axis": true,
   "rim": true,
   "need": 1,
   "seed": 1,
   "on": [
    16
   ]
  }
 },
 "sides": {
  "0": {
   "rep": 0,
   "self": true,
   "need": 2,
   "pts": [
    9,
    10,
    67
   ],
   "tips": [
    5
   ]
  },
  "1": {
   "rep": 1,
   "self": false,
   "need": 3,
   "pts": [
    9,
    13,
    14,
    15,
    7,
    16
   ],
   "tips": [
    7
   ]
  },
  "3": {
   "rep": 3,
   "self": true,
   "need": 2,
   "pts": [
    52,
    32,
    14,
    53,
    46
   ],
   "tips": [
    3
   ]
  },
  "4": {
   "rep": 4,
   "self": false,
   "need": 3,
   "pts": [
    52,
    57,
    48,
    45,
    56,
    55
   ],
   "tips": [
    1,
    5
   ]
  },
  "6": {
   "rep": 6,
   "self": true,
   "need": 2,
   "pts": [
    31,
    13,
    38,
    61
   ],
   "tips": [
    4
   ]
  },
  "7": {
   "rep": 7,
   "self": false,
   "need": 3,
   "pts": [
    31,
    32,
    33,
    20,
    34,
    35,
    36
   ],
   "tips": [
    6
   ]
  },
  "9": {
   "rep": 9,
   "self": true,
   "need": 2,
   "pts": [
    24,
    15,
    29
   ],
   "tips": []
  },
  "10": {
   "rep": 10,
   "self": false,
   "need": 3,
   "pts": [
    24,
    25,
    21,
    26,
    27
   ],
   "tips": []
  },
  "12": {
   "rep": 12,
   "self": true,
   "need": 2,
   "pts": [
    0,
    7
   ],
   "tips": []
  },
  "13": {
   "rep": 13,
   "self": false,
   "need": 3,
   "pts": [
    0,
    1,
    2,
    3
   ],
   "tips": []
  },
  "15": {
   "rep": 15,
   "self": true,
   "need": 2,
   "pts": [
    47,
    48,
    34,
    26,
    3
   ],
   "tips": [
    4
   ]
  },
  "16": {
   "rep": 16,
   "self": false,
   "need": 3,
   "pts": [
    47,
    57,
    32,
    13,
    10,
    68
   ],
   "tips": [
    0,
    6
   ]
  },
  "18": {
   "rep": 18,
   "self": true,
   "need": 2,
   "pts": [
    65,
    45,
    35,
    27
   ],
   "tips": [
    3
   ]
  },
  "19": {
   "rep": 19,
   "self": false,
   "need": 3,
   "pts": [
    65,
    48,
    33,
    24,
    14,
    38,
    67
   ],
   "tips": [
    1
   ]
  },
  "21": {
   "rep": 21,
   "self": true,
   "need": 2,
   "pts": [
    44,
    56,
    36
   ],
   "tips": [
    0
   ]
  },
  "22": {
   "rep": 22,
   "self": false,
   "need": 3,
   "pts": [
    44,
    45,
    34,
    21,
    1,
    7,
    29,
    46
   ],
   "tips": [
    2
   ]
  },
  "24": {
   "rep": 24,
   "self": true,
   "need": 2,
   "pts": [
    20,
    21,
    2,
    16
   ],
   "tips": []
  },
  "25": {
   "rep": 25,
   "self": false,
   "need": 3,
   "pts": [
    20,
    25,
    0,
    15,
    53,
    61
   ],
   "tips": []
  }
 },
 "tips": [
  {
   "i": 0,
   "a": 16,
   "b": 21,
   "axis": false
  },
  {
   "i": 1,
   "a": 4,
   "b": 19,
   "axis": false
  },
  {
   "i": 2,
   "a": 22,
   "b": 22,
   "axis": true
  },
  {
   "i": 3,
   "a": 3,
   "b": 18,
   "axis": false
  },
  {
   "i": 4,
   "a": 6,
   "b": 15,
   "axis": false
  },
  {
   "i": 5,
   "a": 0,
   "b": 4,
   "axis": false
  },
  {
   "i": 6,
   "a": 7,
   "b": 16,
   "axis": false
  },
  {
   "i": 7,
   "a": 1,
   "b": 1,
   "axis": true
  }
 ]
}
""")

nodes = {int(k): v for k, v in DATA['nodes'].items()}
sides = {int(k): v for k, v in DATA['sides'].items()}
tips  = DATA['tips']
for v in nodes.values():
    if v['axis'] and v['rim']:
        v['need'] = 0
        v['seed'] = 0

SID, NID = sorted(sides), sorted(nodes)
obj = [('s', s) for s in SID] + [('n', n) for n in NID] + [('t', t['i']) for t in tips]
IDX = {o: i for i, o in enumerate(obj)}
N = len(obj)
AXTIP = {t['i']: t['axis'] for t in tips}

INC = {}
for s in SID:
    INC[IDX[('s', s)]] = ([IDX[('n', p)] for p in sides[s]['pts']] +
                          [IDX[('t', i)] for i in sides[s]['tips']])
for n in NID:
    INC[IDX[('n', n)]] = [IDX[('s', s)] for s in SID if n in sides[s]['pts']]
for t in tips:
    INC[IDX[('t', t['i'])]] = [IDX[('s', x)] for x in sorted({t['a'], t['b']})]

KIND = [o[0] for o in obj]
NEED = {}
for i, (kind, key) in enumerate(obj):
    if kind == 's':
        NEED[i] = sides[key]['need']
    elif kind == 't':
        NEED[i] = 0 if AXTIP[key] else 1     # an axis tip needs no circle first
    else:
        NEED[i] = nodes[key]['need']
SEEDC = {IDX[('n', n)]: nodes[n]['seed'] for n in NID}

def buildable(i, mask):
    c = sum(1 for j in INC[i] if mask >> j & 1)
    if KIND[i] in ('s', 't'):
        return c == NEED[i]
    return c == NEED[i] or c == 0

def cubes(level):
    cur = {0}
    for _ in range(level):
        nxt = set()
        for m in cur:
            for i in range(N):
                if not (m >> i & 1) and buildable(i, m):
                    nxt.add(m | (1 << i))
        cur = nxt
        if not cur:
            break
    return sorted(cur)

def solve(args):
    mask, level, pin = args
    from z3 import Int, Bool, If, Sum, Distinct, Or, sat, Solver
    pos = [Int('p%d' % i) for i in range(N)]
    S = Solver()
    S.add(Distinct(pos))
    for i, p in enumerate(pos):
        S.add(p >= 0, p < N)
        S.add(p < level if (mask >> i & 1) else p >= level)
    before = lambda i: Sum([If(pos[j] < pos[i], 1, 0) for j in INC[i]])
    seedbits = []
    for i in range(N):
        if KIND[i] in ('s', 't'):
            S.add(before(i) == NEED[i])
        else:
            b = Bool('seed%d' % i)
            S.add(b == (before(i) == 0))
            S.add(Or(b, before(i) == NEED[i]))
            seedbits.append((b, SEEDC[i]))
    if pin is not None:
        S.add(Sum([If(b, c, 0) for b, c in seedbits]) == pin)
    if S.check() == sat:
        m = S.model()
        return sorted(range(N), key=lambda i: m[pos[i]].as_long())
    return None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('-j', type=int, default=cpu_count())
    ap.add_argument('-k', type=int, default=4)
    ap.add_argument('--any', action='store_true')
    a = ap.parse_args()
    t0 = time.time()
    print("cutting at level %d ..." % a.k, flush=True)
    C = cubes(a.k)
    print("  %d cubes, %d workers" % (len(C), a.j), flush=True)
    if not C:
        print("no reachable states at that level - UNSAT outright")
        return
    work = [(m, a.k, None if a.any else 10) for m in C]
    done = 0
    with Pool(a.j) as pool:
        for res in pool.imap_unordered(solve, work, chunksize=8):
            done += 1
            if res is not None:
                pool.terminate()
                print("\nSATISFIABLE after %d cubes (%.0f s).  Build order:\n"
                      % (done, time.time() - t0))
                for k, i in enumerate(res):
                    tag = {'s': 'circle', 'n': 'node  ', 't': 'tip   '}[obj[i][0]]
                    print("  %2d  %s %s" % (k, tag, obj[i][1]))
                return
            if done % 500 == 0 or done == len(C):
                el = time.time() - t0
                print("  %d/%d refuted  %.0f s, ~%.0f s left"
                      % (done, len(C), el, el / done * (len(C) - done)), flush=True)
    print("\nUNSAT over every cube (%.0f s)." % (time.time() - t0))
    print("The second ring cannot be built in either - only held.")

if __name__ == '__main__':
    main()
