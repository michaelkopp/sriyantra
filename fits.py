#!/usr/bin/env python3
"""
fits.py - how far outside itself must a figure reach to be constructed?

THE QUESTION THE EARLIER SCRIPTS ASKED
--------------------------------------
ring1_cubes.py and ring2_cubes.py refuted constructions in which every
incidence follows from the ORDER: a circle drawn through exactly as many points
as determine it, a node marked where exactly as many circles meet. Call those
INTERNAL: they draw nothing that is not a part of the figure.

Chiodo's construction of the straight figure is not internal. At one point he
chooses a parameter so that a coincidence occurs, and makes the choice by an
auxiliary construction - the circle-line-point problem of Apollonius - rather
than by solving numerically. Nothing in the earlier encoding could express
that move, so the earlier unsatisfiable results say nothing about it.

THE RELAXATION
--------------
A FIT is one arranged coincidence. A circle may be drawn through need(s)+1 of
its points, or a node marked where need(n)+1 of its circles already pass, the
extra incidence holding because a parameter was chosen to make it hold. Each
fit raises the seed weight by one and spends one parameter on the condition, so
the dimension is untouched:

    seed weight + placed tips  ==  BASE + fits
    free parameters            ==  seed weight + placed - fits  ==  BASE

with BASE = 4 for the straight figure, 11 for ring 1, 10 for ring 2. A
construction with F fits therefore still sweeps the whole family; only the
route in changes. F = 0 is the internal question, and is unsatisfiable for all
three figures.

WHAT THIS FOUND
---------------
Randomised rollouts - the same procedure that finds a complete internal order
for the ring-free arc figure in about twenty tries - give these upper bounds,
each order replayed against the rules afterwards by an independent checker,
0 violations and the identity exact:

    figure      BASE   fits found   seed weight   free parameters
    straight       4       4              8            8 - 4 =  4
    ring 1        11       6             17           17 - 6 = 11   (+ r = 12)
    ring 2        10       6             16           16 - 6 = 10   (+ r = 11)

Nothing was found below those, but nothing is refuted below them either: the
rollouts are a heuristic. Use the z3 route to settle a particular F.

WHAT IT MEANS, AND WHAT IS STILL MISSING
----------------------------------------
The straight figure is the calibration. It is constructible - Chiodo built it -
and in this bookkeeping it needs four arranged coincidences. Both rings need
six. They are two steps further out than a figure someone has actually drawn,
not in a different class from it, and the constructions found sweep the full
12- and 11-dimensional families rather than a slice of them.

That is an ordering result and only half a construction. Each arranged
coincidence is one equation in one unknown, and the whole thing is straightedge
and compass only if every one of them has a solution of degree a power of two.
Chiodo's are: his reduce to Apollonius, which is quadratic. Whether the rings'
six are is a question about the field the ring condition generates, and neither
this script nor the earlier ones touch it. It is the obvious next thing to ask.

    python3 fits.py --figure ring1 --fits 1 --late 6 --collect 300
    python3 fits.py --figure ring1 --fits 1 --collect 200  # many orders, to score
    python3 fits.py --figure ring1 --fits 6 --rollout      # find an order, seconds
    python3 fits.py --figure ring1 --fits 5                # decide it, z3, an hour
    python3 fits.py --figure straight --fits 4 --rollout   # the calibration
    python3 fits.py --verify ring1 order.json              # replay a saved order

The rollout route needs nothing installed. The deciding route needs z3.
"""
import json, sys, time, random, argparse
from multiprocessing import Pool, cpu_count

DATA = json.loads(r"""
{
 "straight": {
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
 },
 "ring1": {
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
 },
 "ring2": {
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
}
""")
BASE = {'straight': 4, 'ring1': 11, 'ring2': 10}
VERSION = 'fits.py v3: ring vertices may carry the coincidence; UNSAT writes a record'


def build(which, placed=False, fitmax=0):
    D = DATA[which]
    nodes = {int(k): v for k, v in D['nodes'].items()}
    sides = {int(k): v for k, v in D['sides'].items()}
    tips = D.get('tips', [])
    for v in nodes.values():
        if v['axis'] and v['rim']:              # the poles are pinned by the figure
            v['need'] = 0
            v['seed'] = 0
    SID, NID = sorted(sides), sorted(nodes)
    obj = ([('s', s) for s in SID] + [('n', n) for n in NID] +
           [('t', t['i']) for t in tips])
    IDX = {o: i for i, o in enumerate(obj)}
    N = len(obj)
    INC, AX = {}, {}
    for s in SID:
        INC[IDX[('s', s)]] = ([IDX[('n', p)] for p in sides[s]['pts']] +
                              [IDX[('t', i)] for i in sides[s]['tips']])
    for n in NID:
        INC[IDX[('n', n)]] = [IDX[('s', s)] for s in SID if n in sides[s]['pts']]
    for t in tips:
        INC[IDX[('t', t['i'])]] = [IDX[('s', x)] for x in sorted({t['a'], t['b']})]
        AX[IDX[('t', t['i'])]] = bool(t.get('axis', t['a'] == t['b']))
    KIND = [o[0] for o in obj]
    NEED = {i: (sides[k]['need'] if t == 's' else
                (0 if AX[i] else 1) if t == 't' else nodes[k]['need'])
            for i, (t, k) in enumerate(obj)}
    SEEDC = {i: (nodes[k]['seed'] if t == 'n' else 0) for i, (t, k) in enumerate(obj)}
    return dict(obj=obj, N=N, INC=INC, KIND=KIND, NEED=NEED, SEEDC=SEEDC, AX=AX,
                base=BASE[which], placed=placed, fitmax=fitmax, which=which)


def moves(G, i, mask, fits):
    """[(parameters spent, fits used, what it was)] for producing i now"""
    c = sum(1 for j in G['INC'][i] if mask >> j & 1)
    K, need = G['KIND'][i], G['NEED'][i]
    out = []
    if K == 's':
        if c == need:
            out.append((0, 0, ''))
        elif c == need + 1 and fits < G['fitmax']:
            out.append((1, 1, 'fit: an extra point arranged onto it'))
    elif K == 't':
        if c == need:
            out.append((0, 0, ''))
        elif c == 0 and need == 1 and G['placed']:
            out.append((1, 0, 'placed as a free point on E'))
        elif c == need + 1 and fits < G['fitmax']:
            out.append((1, 1, 'fit: an extra circle arranged through it'))
    else:
        if c == need:
            out.append((0, 0, ''))
        elif c == 0:
            out.append((G['SEEDC'][i], 0, 'seed'))
        elif c == need + 1 and fits < G['fitmax']:
            out.append((1, 1, 'fit: an extra circle arranged through it'))
    return out


def rollout(G, rnd, wseed=0.4):
    """build at random, restart on a dead end; returns a complete order or None"""
    N, FULL = G['N'], (1 << G['N']) - 1
    cap = G['base'] + 2 * G['fitmax']
    mask = cost = fits = 0
    order = []
    while mask != FULL:
        opts = []
        for i in range(N):
            if mask >> i & 1:
                continue
            for c, fd, tag in moves(G, i, mask, fits):
                if cost + c <= cap:
                    opts.append((i, c, fd, tag))
        if not opts:
            return None, len(order)
        free = [o for o in opts if o[1] == 0]
        paid = [o for o in opts if o[1] > 0]
        pool = free if (free and (not paid or rnd.random() > wseed)) else (paid or free)
        i, c, fd, tag = rnd.choice(pool)
        mask |= 1 << i
        cost += c
        fits += fd
        order.append((i, c, fd, tag))
    return (order if cost == G['base'] + 2 * fits else None), len(order)


def verify(G, order):
    """replay an order against the rules, written out again from scratch"""
    N, INC, KIND, NEED, SEEDC = G['N'], G['INC'], G['KIND'], G['NEED'], G['SEEDC']
    pos = {}
    for k, rec in enumerate(order):
        m = [j for j, o in enumerate(G['obj']) if o == (rec[0], rec[1])]
        if len(m) != 1:
            return False, 'unknown object %r' % (rec,)
        pos[m[0]] = k
    if len(pos) != N:
        return False, 'not a permutation: %d of %d' % (len(pos), N)
    bad, seedw, fits, placedn = [], 0, 0, 0
    for i in range(N):
        b = sum(1 for j in INC[i] if pos[j] < pos[i])
        if KIND[i] == 't':
            if b == NEED[i]:                       pass
            elif b == 0 and NEED[i] == 1:          placedn += 1
            else: bad.append((G['obj'][i], b, NEED[i]))
        elif KIND[i] == 'n' and b == 0:            seedw += SEEDC[i]
        elif b == NEED[i]:                         pass
        elif b == NEED[i] + 1:                     fits += 1
        else: bad.append((G['obj'][i], b, NEED[i]))
    if bad:
        return False, 'violations: %r' % (bad[:5],)
    if seedw + placedn != G['base'] + fits:
        return False, ('identity fails: %d + %d != %d + %d'
                       % (seedw, placedn, G['base'], fits))
    return True, ('%d objects, 0 violations, seed weight %d, %d placed, %d fits, '
                  'free parameters %d = the family dimension'
                  % (N, seedw, placedn, fits, seedw + placedn - fits))


def cubes(G, level):
    cap = G['base'] + 2 * G['fitmax']
    cur = {(0, 0, 0)}
    for _ in range(level):
        nxt = set()
        for mask, cost, fits in cur:
            for i in range(G['N']):
                if mask >> i & 1:
                    continue
                for c, fd, _ in moves(G, i, mask, fits):
                    if cost + c <= cap:
                        nxt.add((mask | (1 << i), cost + c, fits + fd))
        cur = nxt
        if not cur:
            break
    return sorted({m for m, _, _ in cur})


def solve(args):
    which, placed, fitmax, mask, level, late, d0, fitat = args
    G = build(which, placed, fitmax)
    N, INC, KIND, NEED, SEEDC = G['N'], G['INC'], G['KIND'], G['NEED'], G['SEEDC']
    from z3 import Int, Bool, If, Sum, Distinct, Or, And, Not, Implies, sat, Solver
    pos = [Int('p%d' % i) for i in range(N)]
    S = Solver()
    S.add(Distinct(pos))
    for i, p in enumerate(pos):
        S.add(p >= 0, p < N)
        if mask is not None:            # None means: no cube, the whole question
            S.add(p < level if (mask >> i & 1) else p >= level)
    before = lambda i: Sum([If(pos[j] < pos[i], 1, 0) for j in INC[i]])
    fit, spend, seedbits = [], [], []       # spend: seeds and placed tips
    fitof, seedof = {}, {}
    IDX = {o: i for i, o in enumerate(G['obj'])}
    for i in range(N):
        b = before(i)
        if KIND[i] == 's':
            f = Bool('f%d' % i)
            S.add(b == NEED[i] + If(f, 1, 0))
            fit.append(f)
            fitof[i] = f
        elif KIND[i] == 't':
            if placed and NEED[i] == 1:
                pl = Bool('pl%d' % i)
                S.add(b == If(pl, 0, 1))
                spend.append((pl, 1))
            else:
                f = Bool('f%d' % i)
                S.add(b == NEED[i] + If(f, 1, 0))
                fit.append(f)
                fitof[i] = f
        else:
            sd, f = Bool('s%d' % i), Bool('f%d' % i)
            S.add(sd == (b == 0))
            S.add(Or(sd, b == NEED[i] + If(f, 1, 0)))
            S.add(Implies(sd, Not(f)))      # a seed is not also a fit
            fit.append(f)
            fitof[i] = f
            seedof[i] = sd
            spend.append((sd, SEEDC[i]))
            seedbits.append((sd, i))
    nfit = Sum([If(f, 1, 0) for f in fit])
    S.add(nfit <= fitmax)
    S.add(Sum([If(b, w, 0) for b, w in spend]) == G['base'] + nfit)
    if d0:
        # THE CONDITION THAT MAKES THE LAST STEP DRAWABLE
        # d0 = (C, s): the coincidence is at circle C, and the seed s is the
        # LAST of C's own points to be placed. Then nothing between s and the
        # residual is an intersection, the equation is rational in s, and if
        # its degree is one the step is: draw the circle through the points
        # already there, and take where it meets the locus s lives on.
        # `--late` asked for a seed near the end of the order, which is a
        # different and much stronger thing: the coincidence may sit anywhere.
        C, sd_i = d0
        iC = IDX[('s', C)]
        S.add(fitof[iC])                       # the coincidence is at this circle
        S.add(pos[sd_i] < pos[iC])
        S.add(seedof[sd_i])
        for q in INC[iC]:
            if q == sd_i:
                continue
            S.add(Or(pos[q] < pos[sd_i], pos[q] > pos[iC]))
    if fitat is not None:
        # PIN THE COINCIDENCE TO ONE OBJECT.
        # The peeling certificate names the only objects the coincidence may
        # sit on - eleven per ring - and among them the nodes confined to a
        # locus are the interesting ones, because a fit there reads "these two
        # circles and that locus are concurrent", which is a thing that can
        # sometimes be drawn. Collect orders with the fit pinned there and let
        # concurrent.py measure which of them are drawable. Note this asks
        # nothing about what stands still: standing still is not a property of
        # the order, it is a property of the geometry, and it is measured
        # afterwards rather than guessed at here.
        if IDX[fitat] not in fitof:
            return None            # nothing may be arranged at that object
        S.add(fitof[IDX[fitat]])

    if late:
        # A SEED PLACED LATE IS WHAT MAKES THE LAST STEP DRAWABLE.
        # The coincidence is solved for some parameter, and the equation is
        # rational in that parameter exactly when no intersection lies
        # between the two. Every node that is not a seed comes from an
        # intersection, so the only way is for the parameter to be a seed
        # that is itself a definer of the coincidence, with everything else
        # in it already drawn. That needs the seed placed near the end -
        # which is one constraint, rather than a hope.
        S.add(Or([And(sd, pos[i] >= N - late) for sd, i in seedbits]))
    if S.check() == sat:
        m = S.model()
        return sorted(range(N), key=lambda i: m[pos[i]].as_long())
    return None


def show(G, order):
    for k, rec in enumerate(order):
        t, key = rec[0], rec[1]
        tag = ''
        if len(rec) > 3 and rec[3]:
            tag = 'FIT'
        elif len(rec) > 2 and rec[2]:
            tag = 'seed(%d)' % rec[2]
        print('  %2d  %-6s %-3s %s' % (k, {'s': 'circle', 'n': 'node', 't': 'tip'}[t],
                                       key, tag))


def run_d0(a, G):
    """Every (circle, seed) pair, asked whole.

    Cube-and-conquer exists to split ONE hard refutation across cores. Here
    there are eighty-eight separate questions, each far more constrained than
    the open one, so the cubes are pure overhead - a full sweep per pair would
    take about a day, and nothing would print until the first finished. One
    solver call per pair, parallel over the pairs instead.
    """
    pairs = []
    for i, (k, key) in enumerate(G['obj']):
        if k != 's':
            continue
        for j in G['INC'][i]:
            if G['KIND'][j] == 'n':
                pairs.append((key, j))
    print('  %d (circle, seed) pairs, one solver call each, %d workers'
          % (len(pairs), a.j), flush=True)
    work = [(a.figure, a.placed, a.fits, None, 0, 0, p, None) for p in pairs]
    t0 = time.time()
    done = 0
    with Pool(a.j) as pool:
        for p, res in zip(pairs, pool.imap(solve, work)):
            done += 1
            lab = 'circle %-3d seed node %-3s' % (p[0], G['obj'][p[1]][1])
            if res is not None:
                pool.terminate()
                print('\n  SATISFIABLE for %s after %.0f s\n' % (lab, time.time() - t0))
                order = [[G['obj'][i][0], G['obj'][i][1]] for i in res]
                fn = 'order_%s_d0.json' % a.figure
                json.dump(order, open(fn, 'w'), indent=1)
                show(G, [G['obj'][i] for i in res])
                print('\nwritten to %s - the coincidence is rational in that seed.' % fn)
                return
            print('  %3d/%d  %s  refuted   %.0f s'
                  % (done, len(pairs), lab, time.time() - t0), flush=True)
    print('\nEvery pair refuted: with one arranged coincidence, no order makes '
          'the equation rational in any parameter.')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--figure', default='ring1', choices=('straight', 'ring1', 'ring2'))
    ap.add_argument('--fits', type=int, default=6)
    ap.add_argument('--placed', action='store_true',
                    help='a ring tip may also be marked as a free point on E')
    ap.add_argument('--rollout', action='store_true', help='search by rollout, no z3')
    ap.add_argument('--seconds', type=float, default=120)
    ap.add_argument('--d0', action='store_true',
                    help='ask instead for an order whose coincidence is rational in '
                         'one parameter: a circle whose last-placed point is a seed. '
                         'Every (circle, seed) pair is tried in turn')
    ap.add_argument('--late', type=int, default=0, metavar='K',
                    help='require a seed among the last K objects built. This is '
                         'what makes the coincidence rational in that seed, and so '
                         'possibly drawable rather than solved')
    ap.add_argument('--collect', type=int, default=0, metavar='N',
                    help='keep going and collect N distinct orders instead of '
                         'stopping at the first; written to orders_<figure>_f<F>.json')
    ap.add_argument('--fitat', metavar='OBJ', default=None,
                    help='pin the one coincidence to an object, e.g. n67 or s0. '
                         'Run with --collect to gather orders for concurrent.py; '
                         'without it, UNSAT means the coincidence can never sit '
                         'there at all.')
    ap.add_argument('--verify', nargs=2, metavar=('FIGURE', 'ORDER.json'))
    ap.add_argument('-j', type=int, default=cpu_count())
    ap.add_argument('-k', type=int, default=4)
    a = ap.parse_args()

    if a.verify:
        G = build(a.verify[0], True, 99)
        ok, msg = verify(G, json.load(open(a.verify[1])))
        print(('VALID   ' if ok else 'INVALID ') + msg)
        return

    G = build(a.figure, a.placed, a.fits)
    print('%s: %d objects, %d free parameters, up to %d fit(s)%s'
          % (a.figure, G['N'], G['base'], a.fits,
             ', tips may be placed' if a.placed else ''), flush=True)

    if a.rollout:
        rnd = random.Random(11)
        t0 = time.time()
        k = best = 0
        while time.time() - t0 < a.seconds:
            k += 1
            o, d = rollout(G, rnd)
            best = max(best, d)
            if o:
                f = sum(x[2] for x in o)
                print('\nSATISFIABLE with %d fit(s), %d rollouts, %.1f s\n' % (f, k, time.time() - t0))
                show(G, [(G['obj'][i][0], G['obj'][i][1], c, fd) for i, c, fd, _ in o])
                out = 'order_%s_f%d.json' % (a.figure, f)
                json.dump([[G['obj'][i][0], G['obj'][i][1], c, fd] for i, c, fd, _ in o],
                          open(out, 'w'), indent=1)
                ok, msg = verify(G, json.load(open(out)))
                print('\nreplayed: ' + ('VALID  ' if ok else 'INVALID  ') + msg)
                print('written to ' + out)
                return
        print('nothing in %d rollouts / %.0f s, deepest %d of %d.'
              % (k, time.time() - t0, best, G['N']))
        print('A rollout failing proves nothing; use the z3 route to decide.')
        return

    try:
        import z3                                       # noqa: F401
    except ImportError:
        sys.exit('needs z3 for this route: pip install z3-solver,'
                 ' or pass --rollout to search without it')
    print(VERSION, flush=True)
    FITAT = None
    if a.fitat:
        k, key = a.fitat[0], int(a.fitat[1:])
        FITAT = (k, key)
        if FITAT not in G['obj']:
            print('no such object: %s' % a.fitat); return
        if a.fits < 1:
            print('--fitat needs --fits at least 1'); return
        print('  the coincidence is pinned to %s' % a.fitat, flush=True)
    globals()['FITAT'] = FITAT
    if a.d0:
        return run_d0(a, G)          # no cubes: each pair is one whole question
    t0 = time.time()
    C = cubes(G, a.k)
    print('  %d cubes at level %d, %d workers' % (len(C), a.k, a.j), flush=True)
    if not C:
        print('no reachable state at that level - unsatisfiable outright.')
        return
    work = [(a.figure, a.placed, a.fits, m, a.k, a.late, None, FITAT) for m in C]
    done = 0
    found = []
    with Pool(a.j) as pool:
        for res in pool.imap_unordered(solve, work, chunksize=8):
            done += 1
            if res is not None:
                order = [[G['obj'][i][0], G['obj'][i][1]] for i in res]
                if not a.collect:
                    pool.terminate()
                    print('\nSATISFIABLE after %d cubes (%.0f s):\n' % (done, time.time() - t0))
                    show(G, [G['obj'][i] for i in res])
                    return
                found.append(order)
                if len(found) % 10 == 0 or len(found) == a.collect:
                    print('  %d orders, %d cubes, %.0f s'
                          % (len(found), done, time.time() - t0), flush=True)
                if len(found) >= a.collect:
                    pool.terminate()
                    break
            if done % 1000 == 0 or done == len(C):
                el = time.time() - t0
                print('  %d/%d cubes refuted  %.0f s, ~%.0f s left'
                      % (done, len(C), el, el / done * (len(C) - done)), flush=True)
    if a.collect:
        out = ('orders_%s_f%d%s.json'
               % (a.figure, a.fits, '_' + a.fitat if a.fitat else ''))
        if not found:
            # An empty list is ambiguous: the old solver produced one in
            # minutes by refuting every cube for a bookkeeping reason, and a
            # real exhaustive UNSAT produces one too. Say which this is.
            json.dump(dict(unsat=True, figure=a.figure, fits=a.fits, fitat=a.fitat,
                           cubes=done, level=a.k, seconds=round(time.time() - t0)),
                      open(out, 'w'), indent=1)
            print('\nUNSAT: no order with %d fit(s)%s exists. %d cubes, %.0f s. '
                  'Written to %s as a record.'
                  % (a.fits, ' pinned at ' + a.fitat if a.fitat else '', done,
                     time.time() - t0, out))
            return
        json.dump(found, open(out, 'w'))
        print('\n%d distinct orders written to %s (%d cubes, %.0f s).'
              % (len(found), out, done, time.time() - t0))
        print('They are all the same mathematics; what differs is how the')
        print('construction behaves in arithmetic. Send the file and I will')
        print('score them for how far the sliders move before the steps')
        print('stop being defined, and keep the best-conditioned one.')
        return
    print('\nUNSAT over every cube (%.0f s): %d arranged coincidence(s) are not enough.'
          % (time.time() - t0, a.fits))


if __name__ == '__main__':
    main()
