#!/usr/bin/env python3
"""
chiodo_cubes.py - the same question as ring1_cubes.py, asked of the STRAIGHT
                  figure: can Chiodo's Sri Yantra be built by incidence alone?

THE ANSWER IS NO, AND THAT IS NOT A CRITICISM OF CHIODO
-------------------------------------------------------
This script was written expecting SATISFIABLE and returns UNSAT. What that
means is set out at the bottom; read it before quoting the result.

THE QUESTION
------------
Fix the circle and its vertical diameter. Allow exactly three moves:

    draw the line through two points already marked;
    mark the point where two lines already drawn cross;
    mark the point where a line already drawn crosses the circle, or the axis.

Every object of the figure - eighteen side orbits and thirty-nine node orbits,
counted up to the mirror symmetry - must be produced by one of these, once, in
some order, and every one of the eighty-eight incidence orbits must follow from
the order rather than be a coincidence. So a side is drawn at a moment when
EXACTLY need(s) of the points on it exist - one for a base, which is horizontal
and therefore fixed by a single point, two for a general side - and a point is
marked when EXACTLY need(n) of the sides through it exist, or when none do, in
which case it is a seed and costs parameters.

Chiodo's condition (i) needs no separate encoding: the six corners of t3 and t7
are the figure's rim points, and a rim point is marked where one line crosses
the circle. Putting them there is a move, not a wish.

THE BUDGET IS THE FIGURE'S OWN, NOT ONE IMPOSED ON IT
-----------------------------------------------------
This is the identical encoding used for the arc figure. Two numbers move, and
both because a straight side is a LINE rather than a circle:

                                    arcs      straight
    definer slots, sum of need(s)     45            27
    node budget S, sum of need(n)     61            65
    incidence orbits                  88            88

Each incidence is consumed exactly once - by the side, if the point came first,
or by the point, if the side did - so the number of free parameters is forced:

    seeds = S - incidences + slots
          = 61 - 88 + 45 = 18    for the arc figure   (its recipe has 18)
          = 65 - 88 + 27 =  4    for the straight one (Chiodo has 4: a,b,c,d)

The straight figure's own arithmetic reproduces Chiodo's four heights, so the
question is being asked with the right budget. Note also that the pin is
redundant: any complete order has seed cost exactly 4 whether or not it is
asserted. --any drops it and changes nothing.

WHAT COMES BACK
---------------
UNSAT, in seconds rather than the hour the ring questions take, because four
parameters is a hard bound on any prefix and the sweep dies quickly. It dies at
level 32: no order reaches more than 31 of the 57 objects, whatever the four
parameters are spent on. The single branch that gets to 31 spends them on two
free points and leaves eight sides and eighteen nodes unreachable.

THE MACHINERY IS KNOWN TO FIND CONSTRUCTIONS THAT EXIST
--------------------------------------------------------
Two controls, both already passed:

  * the viewer's real arc recipe - twelve seeds and forty-five steps - replays
    through this encoding object by object with 0 violations, 57 objects,
    complete, seed cost 18, exactly as the identity above predicts;
  * classical_cubes.py puts the arc figure through this same cube-and-conquer
    machinery and returns SATISFIABLE.

and the two ways of asking the question here - the reachability sweep and the
constraints the z3 formula asserts - were checked against each other by brute
force over every permutation of 4,200 small random instances: no disagreement.

WHAT THE RESULT DOES AND DOES NOT SAY
-------------------------------------
It does NOT say the straight Sri Yantra is not constructible. Chiodo built it.
His construction is correct and is ruler and compass; its essential step is
Apollonius' circle-line-point problem, solved with an auxiliary circle that is
not part of the figure.

It says the straight figure cannot be propagated. There is no order in which
its own sides and corners produce one another, so any construction of it must
step outside the figure - Chiodo's auxiliary circle, or the two parabolas the
viewer meets to fix one line. The arc figure needs no such step: its forty-five
moves are all "circle through three points" and "where two circles cross", and
they close.

That is the sharper form of what the arc family buys. Not merely that the
classical concurrences survive the change of curve, but that the figure becomes
self-constructing: rigid enough to be exact, loose enough to be built.

    pip install z3-solver
    python3 chiodo_cubes.py            # cube and conquer, all cores
    python3 chiodo_cubes.py --sweep    # the whole proof, no solver, ~5 s
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
   "dof": 2,
   "need": 2,
   "seed": 2,
   "on": [
    12,
    13,
    25
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
    7,
    12,
    22
   ]
  },
  "4": {
   "rep": 4,
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
  "5": {
   "rep": 5,
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
  "6": {
   "rep": 6,
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
  "9": {
   "rep": 9,
   "axis": false,
   "rim": false,
   "dof": 2,
   "need": 2,
   "seed": 2,
   "on": [
    4,
    24,
    25
   ]
  },
  "11": {
   "rep": 11,
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
  "13": {
   "rep": 13,
   "axis": true,
   "rim": false,
   "dof": 1,
   "need": 1,
   "seed": 1,
   "on": [
    7,
    24
   ]
  },
  "14": {
   "rep": 14,
   "axis": true,
   "rim": false,
   "dof": 1,
   "need": 1,
   "seed": 1,
   "on": [
    3,
    25
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
    7,
    9,
    25
   ]
  },
  "16": {
   "rep": 16,
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
  "17": {
   "rep": 17,
   "axis": false,
   "rim": false,
   "dof": 2,
   "need": 2,
   "seed": 2,
   "on": [
    0,
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
    9,
    10,
    19
   ]
  },
  "23": {
   "rep": 23,
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
  "25": {
   "rep": 25,
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
  "28": {
   "rep": 28,
   "axis": false,
   "rim": false,
   "dof": 2,
   "need": 2,
   "seed": 2,
   "on": [
    21,
    22
   ]
  },
  "30": {
   "rep": 30,
   "axis": true,
   "rim": false,
   "dof": 1,
   "need": 1,
   "seed": 1,
   "on": [
    4,
    21
   ]
  },
  "31": {
   "rep": 31,
   "axis": false,
   "rim": false,
   "dof": 2,
   "need": 2,
   "seed": 2,
   "on": [
    1,
    21
   ]
  },
  "33": {
   "rep": 33,
   "axis": true,
   "rim": false,
   "dof": 1,
   "need": 1,
   "seed": 1,
   "on": [
    0,
    22
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
    1,
    18,
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
    4,
    15,
    22
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
    7
   ]
  },
  "40": {
   "rep": 40,
   "axis": true,
   "rim": false,
   "dof": 1,
   "need": 1,
   "seed": 1,
   "on": [
    6,
    19
   ]
  },
  "41": {
   "rep": 41,
   "axis": false,
   "rim": false,
   "dof": 2,
   "need": 2,
   "seed": 2,
   "on": [
    6,
    16
   ]
  },
  "43": {
   "rep": 43,
   "axis": false,
   "rim": false,
   "dof": 2,
   "need": 2,
   "seed": 2,
   "on": [
    0,
    7,
    19
   ]
  },
  "44": {
   "rep": 44,
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
  "47": {
   "rep": 47,
   "axis": false,
   "rim": false,
   "dof": 2,
   "need": 2,
   "seed": 2,
   "on": [
    18,
    19
   ]
  },
  "49": {
   "rep": 49,
   "axis": false,
   "rim": false,
   "dof": 2,
   "need": 2,
   "seed": 2,
   "on": [
    4,
    18
   ]
  },
  "51": {
   "rep": 51,
   "axis": false,
   "rim": false,
   "dof": 2,
   "need": 2,
   "seed": 2,
   "on": [
    3,
    19
   ]
  },
  "52": {
   "rep": 52,
   "axis": false,
   "rim": false,
   "dof": 2,
   "need": 2,
   "seed": 2,
   "on": [
    4,
    19
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
    1,
    15,
    19
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
    3,
    4
   ]
  },
  "59": {
   "rep": 59,
   "axis": false,
   "rim": false,
   "dof": 2,
   "need": 2,
   "seed": 2,
   "on": [
    0,
    4,
    16
   ]
  },
  "61": {
   "rep": 61,
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
  "63": {
   "rep": 63,
   "axis": true,
   "rim": true,
   "dof": 1,
   "need": 0,
   "seed": 0,
   "on": [
    16
   ]
  },
  "64": {
   "rep": 64,
   "axis": false,
   "rim": false,
   "dof": 2,
   "need": 2,
   "seed": 2,
   "on": [
    1,
    16
   ]
  },
  "66": {
   "rep": 66,
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
  "68": {
   "rep": 68,
   "axis": true,
   "rim": true,
   "dof": 1,
   "need": 0,
   "seed": 0,
   "on": [
    1
   ]
  }
 },
 "sides": {
  "0": {
   "rep": 0,
   "self": true,
   "need": 1,
   "pts": [
    66,
    59,
    43,
    17,
    33
   ],
   "tips": []
  },
  "1": {
   "rep": 1,
   "self": false,
   "need": 2,
   "pts": [
    66,
    64,
    53,
    34,
    31,
    68
   ],
   "tips": []
  },
  "3": {
   "rep": 3,
   "self": true,
   "need": 1,
   "pts": [
    57,
    44,
    51,
    14
   ],
   "tips": []
  },
  "4": {
   "rep": 4,
   "self": false,
   "need": 2,
   "pts": [
    57,
    59,
    52,
    9,
    35,
    49,
    30
   ],
   "tips": []
  },
  "6": {
   "rep": 6,
   "self": true,
   "need": 1,
   "pts": [
    38,
    41,
    40
   ],
   "tips": []
  },
  "7": {
   "rep": 7,
   "self": false,
   "need": 2,
   "pts": [
    38,
    44,
    43,
    15,
    2,
    13
   ],
   "tips": []
  },
  "9": {
   "rep": 9,
   "self": true,
   "need": 1,
   "pts": [
    21,
    15,
    23
   ],
   "tips": []
  },
  "10": {
   "rep": 10,
   "self": false,
   "need": 2,
   "pts": [
    21,
    16,
    11,
    26,
    25
   ],
   "tips": []
  },
  "12": {
   "rep": 12,
   "self": true,
   "need": 1,
   "pts": [
    0,
    2
   ],
   "tips": []
  },
  "13": {
   "rep": 13,
   "self": false,
   "need": 2,
   "pts": [
    0,
    6,
    5,
    4
   ],
   "tips": []
  },
  "15": {
   "rep": 15,
   "self": true,
   "need": 1,
   "pts": [
    61,
    53,
    35,
    26,
    4
   ],
   "tips": []
  },
  "16": {
   "rep": 16,
   "self": false,
   "need": 2,
   "pts": [
    61,
    64,
    59,
    44,
    41,
    63
   ],
   "tips": []
  },
  "18": {
   "rep": 18,
   "self": true,
   "need": 1,
   "pts": [
    47,
    34,
    49,
    25
   ],
   "tips": []
  },
  "19": {
   "rep": 19,
   "self": false,
   "need": 2,
   "pts": [
    47,
    53,
    52,
    21,
    43,
    51,
    40
   ],
   "tips": []
  },
  "21": {
   "rep": 21,
   "self": true,
   "need": 1,
   "pts": [
    28,
    31,
    30
   ],
   "tips": []
  },
  "22": {
   "rep": 22,
   "self": false,
   "need": 2,
   "pts": [
    28,
    34,
    35,
    11,
    6,
    2,
    23,
    33
   ],
   "tips": []
  },
  "24": {
   "rep": 24,
   "self": true,
   "need": 1,
   "pts": [
    9,
    11,
    5,
    13
   ],
   "tips": []
  },
  "25": {
   "rep": 25,
   "self": false,
   "need": 2,
   "pts": [
    9,
    16,
    0,
    15,
    17,
    14
   ],
   "tips": []
  }
 },
 "tips": []
}
""")

nodes = {int(k): v for k, v in DATA['nodes'].items()}
sides = {int(k): v for k, v in DATA['sides'].items()}
for v in nodes.values():
    if v['axis'] and v['rim']:      # the two poles are pinned by the figure
        v['need'] = 0
        v['seed'] = 0

SID, NID = sorted(sides), sorted(nodes)
obj = [('s', s) for s in SID] + [('n', n) for n in NID]
IDX = {o: i for i, o in enumerate(obj)}
N = len(obj)

INC = {}
for s in SID:
    INC[IDX[('s', s)]] = [IDX[('n', p)] for p in sides[s]['pts']]
for n in NID:
    INC[IDX[('n', n)]] = [IDX[('s', s)] for s in SID if n in sides[s]['pts']]

KIND = [o[0] for o in obj]
NEED = {i: (sides[k]['need'] if t == 's' else nodes[k]['need'])
        for i, (t, k) in enumerate(obj)}
SEEDC = {i: (nodes[k]['seed'] if t == 'n' else 0) for i, (t, k) in enumerate(obj)}

S_TOT = sum(nodes[n]['need'] for n in NID)
SLOTS = sum(sides[s]['need'] for s in SID)
INCID = sum(len(sides[s]['pts']) for s in SID)
BUDGET = S_TOT - INCID + SLOTS          # forced by the figure, not chosen


def step(i, mask):
    """can object i be built in state mask, and what does it cost?"""
    c = sum(1 for j in INC[i] if mask >> j & 1)
    if KIND[i] == 's':
        return (c == NEED[i], 0)
    if c == NEED[i]:
        return (True, 0)
    return (c == 0, SEEDC[i])           # a seed, if nothing through it yet


def frontier(level=None, verbose=False):
    """the reachable states, level by level.

    A state is the SET built so far together with what has been spent on
    seeds; a branch is dropped once it has spent more than the budget, which
    is sound because the total is forced to exactly BUDGET, so an overspent
    prefix can never complete. Stops at `level`, or when nothing is left."""
    cur = {(0, 0)}
    lev = 0
    while cur and (level is None or lev < level):
        nxt = set()
        for mask, cost in cur:
            for i in range(N):
                if mask >> i & 1:
                    continue
                ok, c = step(i, mask)
                if ok and cost + c <= BUDGET:
                    nxt.add((mask | (1 << i), cost + c))
        lev += 1
        if verbose:
            print('  level %2d  %8d states' % (lev, len(nxt)), flush=True)
        if not nxt:
            return lev - 1, cur
        cur = nxt
    return lev, cur


def cubes(level):
    """the reachable sets of `level` objects - disjoint, and exhaustive"""
    lev, cur = frontier(level)
    return sorted({m for m, _ in cur}) if lev == level else []


def solve(args):
    """z3: is there a complete order whose first `level` objects are the cube?"""
    mask, level, pin = args
    from z3 import Int, Bool, If, Sum, Distinct, Or, sat, Solver
    pos = [Int('p%d' % i) for i in range(N)]
    S = Solver()
    S.add(Distinct(pos))
    for i, p in enumerate(pos):
        S.add(p >= 0, p < N)
        S.add(p < level if (mask >> i & 1) else p >= level)     # the cube
    before = lambda i: Sum([If(pos[j] < pos[i], 1, 0) for j in INC[i]])
    seedbits = []
    for i in range(N):
        if KIND[i] == 's':
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


def header():
    print("Chiodo's straight figure, by incidence alone")
    print("  %d objects: %d side orbits, %d node orbits" % (N, len(SID), len(NID)))
    print("  S = %d   slots = %d   incidences = %d" % (S_TOT, SLOTS, INCID))
    print("  => seeds = %d - %d + %d = %d parameters, which is Chiodo's a, b, c, d"
          % (S_TOT, INCID, SLOTS, BUDGET), flush=True)


def verdict(lev, last):
    print("\nUNSAT. The frontier dies at level %d: no order reaches more than" % (lev + 1))
    print("%d of the %d objects, whatever the four parameters are spent on." % (lev, N))
    m = sorted(last)[0][0]
    miss = [obj[i] for i in range(N) if not (m >> i & 1)]
    print("\nIn the deepest branch these are never reached:")
    print("  sides ", [k for t, k in miss if t == 's'])
    print("  nodes ", [k for t, k in miss if t == 'n'])
    print("\nChiodo's construction exists and is ruler and compass. It simply")
    print("cannot be a construction of THIS kind: it leaves the figure, for an")
    print("auxiliary Apollonius circle, and the result above says it has to.")
    print("The arc figure never leaves - see classical_cubes.py, which returns")
    print("SATISFIABLE on this same machinery.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('-j', type=int, default=cpu_count(), help='workers')
    ap.add_argument('-k', type=int, default=4, help='cube level')
    ap.add_argument('--sweep', action='store_true',
                    help='settle it by exhaustive sweep instead - no z3 needed')
    ap.add_argument('--any', action='store_true',
                    help='drop the seed pin (changes nothing: the count is forced)')
    a = ap.parse_args()
    header()
    if not a.sweep:
        try:
            import z3                                   # noqa: F401
        except ImportError:
            sys.exit("\nneeds z3 for this route:  pip install z3-solver\n"
                     "or run  python3 chiodo_cubes.py --sweep,  which settles the\n"
                     "same question by exhaustion in a few seconds and needs nothing.")

    if a.sweep:
        print("\nsweeping every reachable state:")
        t0 = time.time()
        lev, last = frontier(verbose=True)
        if lev == N:
            print("\nSATISFIABLE - the sweep reached all %d objects (%.1f s)." % (N, time.time() - t0))
            return
        verdict(lev, last)
        print("\n(%.1f s, no solver used)" % (time.time() - t0))
        return

    t0 = time.time()
    print("\ncutting the search at level %d ..." % a.k, flush=True)
    C = cubes(a.k)
    print("  %d cubes, %d workers" % (len(C), a.j), flush=True)
    if not C:
        print("\nno reachable state at that level at all - UNSAT outright.")
        print("Run --sweep to see where the frontier dies.")
        return
    work = [(m, a.k, None if a.any else BUDGET) for m in C]
    done = 0
    with Pool(a.j) as pool:
        for res in pool.imap_unordered(solve, work, chunksize=8):
            done += 1
            if res is not None:
                pool.terminate()
                print("\nSATISFIABLE after %d cubes (%.0f s) - a build order:\n"
                      % (done, time.time() - t0))
                for k, i in enumerate(res):
                    print('  %2d  %s %s' % (k, 'side' if obj[i][0] == 's' else 'node', obj[i][1]))
                return
            if done % 1000 == 0 or done == len(C):
                el = time.time() - t0
                print("  %d/%d cubes refuted  %.0f s elapsed, ~%.0f s left"
                      % (done, len(C), el, el / done * (len(C) - done)), flush=True)
    print("\nUNSAT over every cube (%.0f s)." % (time.time() - t0))
    lev, last = frontier()
    verdict(lev, last)


if __name__ == '__main__':
    main()
