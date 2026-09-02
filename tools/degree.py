#!/usr/bin/env python3
"""How hard is the arranged coincidence to SOLVE, as opposed to satisfy?

The construction reaches the residual from the seeds by two kinds of step:

    circle through points   rational in what it is given
    intersection            one square root

So if the coincidence is solved for a parameter t, the field the answer lives
in is a tower over Q(t) with one quadratic floor per intersection that DEPENDS
on t. Intersections that do not involve t contribute constants and cost
nothing. Hence:

    d(t) = the number of t-dependent intersections between t and the residual

    d = 0   the residual is rational in t; the equation is polynomial, and if
            its degree is 1 or 2 the step is a ruler-and-compass one outright
    d <= 2  the answer lies in a tower of at most two quadratics over the
            rational function field, so it is constructible either way
    large   nothing follows; the degree may still be a power of two, but this
            argument does not give it

This counts d for every free parameter of every order, which is a graph
question and costs nothing. A single (order, parameter) with small d turns the
last numerical step of the construction into a drawable one.
"""
import json, math, re, sys
from collections import defaultdict

SQRT_OPS = {'na', 'nr', 'nx', 't'}          # 'ta' is E meeting the axis: exact


def dataflow(rb):
    """producer of each object, and what each step consumes"""
    prod, cons = {}, {}
    seeds = []
    k = 0
    for s in rb['seeds']:
        key = ('n', s['i'])
        prod[key] = ('seed', len(seeds))
        if s['kind'] == 'pole':
            seeds.append((key, 0))
        elif s['kind'] in ('axis', 'rim'):
            seeds.append((key, 1))
        else:
            seeds.append((key, 2))
    for j, st in enumerate(rb['steps']):
        op = st[0]
        if op in ('s2', 's3'):
            out = ('s', st[1])
            args = st[2:4] if op == 's2' else st[2:5]
            ins = [(('t', a[0]) if isinstance(a, list) else ('n', a)) for a in args]
        elif op == 'na':
            out, ins = ('n', st[1]), [('s', st[2])]
        elif op == 'nr':
            out, ins = ('n', st[1]), [('s', st[2])]
        elif op == 'nx':
            out, ins = ('n', st[1]), [('s', st[2]), ('s', st[3])]
        elif op == 't':
            out, ins = ('t', st[1]), [('s', st[2])]
        else:                                   # 'ta': (0, +-R), no square root
            out, ins = ('t', st[1]), []
        prod[out] = (op, j)
        cons[out] = ins
    return prod, cons, seeds


def cone(rb, prod, cons, start):
    """everything downstream of `start`, and the square-root steps among them"""
    seen, stack = set(), [start]
    while stack:
        o = stack.pop()
        if o in seen:
            continue
        seen.add(o)
        for out, ins in cons.items():
            if o in ins and out not in seen:
                stack.append(out)
    # the mirror of a circle is the mirror of its source, so it carries the
    # same dependence; the same for a mirrored node
    return seen


def depth_for(rb, param):
    """d(t) for one free parameter: t-dependent intersections reaching the fit"""
    prod, cons, seeds = dataflow(rb)
    f = rb['fit']
    target = [('s', f['circle'])]
    target.append(('t', f['point'][0]) if f['kind'] == 's' else ('n', f['node']))

    if param == 'R':
        dep = set()
        for st in rb['steps']:
            if st[0] in ('t', 'ta'):
                dep.add(('t', st[1]))
    else:
        dep = {('n', rb['seeds'][param]['i'])}

    # Push dependence forward, closing under the mirror AS WE GO and iterating
    # to a fixpoint. Closing at the end instead is wrong: a step is executed
    # once and its mirror written at the same moment, so a mirrored object
    # that only becomes dependent through the closure still has to be pushed
    # forward again. Doing it in one pass undercounts, and it undercounted
    # badly - it reported d = 1 for an order whose true depth, measured by
    # finite differences, is 7.
    mir, smir = rb['mir'], rb['smir']

    def mark(o, dep):
        if o in dep:
            return False
        dep.add(o)
        k, i = o
        if k == 'n' and mir[i] != i:
            dep.add(('n', mir[i]))
        if k == 's' and smir[i] != i:
            dep.add(('s', smir[i]))
        return True

    for o in list(dep):
        mark(o, dep)
    changed = True
    while changed:
        changed = False
        for st in rb['steps']:
            op = st[0]
            if op in ('s2', 's3'):
                out = ('s', st[1])
                args = st[2:4] if op == 's2' else st[2:5]
                ins = [(('t', a[0]) if isinstance(a, list) else ('n', a)) for a in args]
            elif op == 'na':
                out, ins = ('n', st[1]), [('s', st[2])]
            elif op == 'nr':
                out, ins = ('n', st[1]), [('s', st[2])]
            elif op == 'nx':
                out, ins = ('n', st[1]), [('s', st[2]), ('s', st[3])]
            elif op == 't':
                out, ins = ('t', st[1]), [('s', st[2])]
            else:
                out, ins = ('t', st[1]), []
            if any(i in dep for i in ins) and mark(out, dep):
                changed = True
    if not any(t in dep for t in target):
        return None                    # the residual does not see this parameter
    # only the intersections that BOTH depend on t and feed the residual count:
    # a square root downstream of the residual is not in the equation at all
    ins_of = {}
    for st in rb['steps']:
        op = st[0]
        if op in ('s2', 's3'):
            out = ('s', st[1])
            args = st[2:4] if op == 's2' else st[2:5]
            ins_of[out] = [(('t', a[0]) if isinstance(a, list) else ('n', a)) for a in args]
        elif op in ('na', 'nr'):
            ins_of[('n', st[1])] = [('s', st[2])]
        elif op == 'nx':
            ins_of[('n', st[1])] = [('s', st[2]), ('s', st[3])]
        elif op == 't':
            ins_of[('t', st[1])] = [('s', st[2])]
        else:
            ins_of[('t', st[1])] = []
    back, stack = set(), list(target)
    while stack:
        o = stack.pop()
        if o in back:
            continue
        back.add(o)
        for i in ins_of.get(o, []):
            stack.append(i)
        if o[0] == 'n' and mir[o[1]] != o[1]:
            stack.append(('n', mir[o[1]]))
        if o[0] == 's' and smir[o[1]] != o[1]:
            stack.append(('s', smir[o[1]]))
    n = 0
    for st in rb['steps']:
        if st[0] in SQRT_OPS:
            out = ('t', st[1]) if st[0] == 't' else ('n', st[1])
            if out in dep and out in back:
                n += 1
    return n


def seed_param_index(rb):
    """map each free parameter slot to the seed it belongs to"""
    out, k = [], 0
    for j, s in enumerate(rb['seeds']):
        if s['kind'] == 'pole':
            continue
        if s['kind'] == 'free':
            out.append(j); out.append(j)
        else:
            out.append(j)
    return out
