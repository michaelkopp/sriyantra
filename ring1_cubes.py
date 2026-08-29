#!/usr/bin/env python3
"""
ring1_cubes.py - the same question as ring1_sat.py, split across your cores.

WHY THIS PARALLELISES WELL
--------------------------
The expected answer is UNSAT, and proving UNSAT means refuting every branch.
That is embarrassingly parallel: cut the search into disjoint pieces, hand one
to each core, and the pieces never need to talk to each other.

The cut is "cube and conquer".  Any valid construction order passes through
exactly one state after its first k objects - the SET of those k, not their
order - and that set has to be reachable under the buildability rules.  So:

    enumerate every reachable state at level k   (a breadth-first sweep)
    -> each is one cube, they are disjoint, and together they are exhaustive
    -> solve each independently: pos[i] < k for i in the cube, >= k otherwise

Cube counts, measured on this problem:

    level 3      9,187 cubes        level 5    613,930
    level 4     84,211 cubes        level 6  3,739,873

Level 4 is the sweet spot for a desktop: enough cubes to keep every core busy
and balance the load, each still big enough that solver startup is not the
cost.  If it comes back SAT the first core to find an order wins and the rest
are cancelled; if UNSAT, the work divides by the number of cores almost exactly.

    pip install z3-solver
    python3 ring1_cubes.py             # all cores, level 4
    python3 ring1_cubes.py -j 8 -k 5   # 8 workers, finer cut
"""
import json, sys, os, argparse, time
from multiprocessing import Pool, cpu_count

DATA = json.loads(r"""
{
 "nodes": {
  "0": {
   "rep": 0,
   "axis": false,
   "rim": false,
   "dof": 2,
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
   "dof": 2,
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
   "dof": 2,
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
   "dof": 1,
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
   "dof": 2,
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
   "dof": 2,
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
   "dof": 2,
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
   "dof": 2,
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
   "dof": 2,
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
   "dof": 2,
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
   "dof": 1,
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
   "dof": 2,
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
   "dof": 2,
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
   "dof": 2,
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
   "dof": 2,
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
   "dof": 2,
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
   "dof": 1,
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
   "dof": 2,
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
   "dof": 2,
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
   "dof": 2,
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
   "dof": 2,
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
   "dof": 2,
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
   "dof": 2,
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
   "dof": 1,
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
   "dof": 2,
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
   "dof": 2,
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
   "dof": 2,
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
   "dof": 1,
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
   "dof": 2,
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
   "dof": 2,
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
   "dof": 2,
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
   "dof": 2,
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
   "dof": 1,
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
   "dof": 2,
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
   "dof": 2,
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
   "dof": 1,
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
   "dof": 2,
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
   "dof": 1,
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
   "dof": 1,
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
    6
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
    2,
    4
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
    5
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
    3
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
    0
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
    1
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
   "b": 18
  },
  {
   "i": 1,
   "a": 21,
   "b": 19
  },
  {
   "i": 2,
   "a": 22,
   "b": 4
  },
  {
   "i": 3,
   "a": 3,
   "b": 15
  },
  {
   "i": 4,
   "a": 4,
   "b": 6
  },
  {
   "i": 5,
   "a": 0,
   "b": 7
  },
  {
   "i": 6,
   "a": 16,
   "b": 1
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

INC = {}
for s in SID:
    INC[IDX[('s', s)]] = ([IDX[('n', p)] for p in sides[s]['pts']] +
                          [IDX[('t', i)] for i in sides[s]['tips']])
for n in NID:
    INC[IDX[('n', n)]] = [IDX[('s', s)] for s in SID if n in sides[s]['pts']]
for t in tips:
    INC[IDX[('t', t['i'])]] = [IDX[('s', x)] for x in sorted({t['a'], t['b']})]

NEED = {}
KIND = [o[0] for o in obj]
for i, (kind, key) in enumerate(obj):
    NEED[i] = sides[key]['need'] if kind == 's' else (1 if kind == 't' else nodes[key]['need'])
SEEDC = {IDX[('n', n)]: nodes[n]['seed'] for n in NID}

def buildable(i, mask):
    c = sum(1 for j in INC[i] if mask >> j & 1)
    if KIND[i] == 's':
        return c == NEED[i]
    if KIND[i] == 't':
        return c == 1
    return c == NEED[i] or c == 0

def cubes(level):
    """every reachable set of `level` objects - disjoint and exhaustive"""
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
    mask, level, pin_cost = args
    from z3 import Int, Bool, If, Sum, Distinct, Or, sat, unsat, Solver
    pos = [Int('p%d' % i) for i in range(N)]
    S = Solver()
    S.add(Distinct(pos))
    for i, p in enumerate(pos):
        S.add(p >= 0, p < N)
        S.add(p < level if (mask >> i & 1) else p >= level)   # the cube
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
    if pin_cost is not None:
        S.add(Sum([If(b, c, 0) for b, c in seedbits]) == pin_cost)
    r = S.check()
    if r == sat:
        m = S.model()
        return sorted(range(N), key=lambda i: m[pos[i]].as_long())
    return None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('-j', type=int, default=cpu_count(), help='workers')
    ap.add_argument('-k', type=int, default=4, help='cube level')
    ap.add_argument('--any', action='store_true',
                    help='drop the 11-seed requirement (an easier question)')
    a = ap.parse_args()
    t0 = time.time()
    print("cutting the search at level %d ..." % a.k, flush=True)
    C = cubes(a.k)
    print("  %d cubes, %d workers" % (len(C), a.j), flush=True)
    if not C:
        print("no reachable states at that level - UNSAT outright")
        return
    work = [(m, a.k, None if a.any else 11) for m in C]
    done = 0
    with Pool(a.j) as pool:
        for res in pool.imap_unordered(solve, work, chunksize=8):
            done += 1
            if res is not None:
                pool.terminate()
                print("\nSATISFIABLE after %d cubes (%.0f s) - an exact construction exists.\n"
                      % (done, time.time() - t0))
                for k, i in enumerate(res):
                    tag = {'s': 'circle', 'n': 'node  ', 't': 'tip   '}[obj[i][0]]
                    print("  %2d  %s %s" % (k, tag, obj[i][1]))
                return
            if done % 500 == 0 or done == len(C):
                el = time.time() - t0
                print("  %d/%d cubes refuted  %.0f s elapsed, ~%.0f s left"
                      % (done, len(C), el, el / done * (len(C) - done)), flush=True)
    print("\nUNSAT over every cube (%.0f s)." % (time.time() - t0))
    print("No exact construction puts the first ring on a circle.")
    print("It can be held there, but not built in - that is a theorem about the figure.")

if __name__ == '__main__':
    main()
