#!/usr/bin/env python3
"""The certificate that no internal construction places a ring on a circle.

Whatever is drawn LAST has every one of its incidences already present, so it
can have at most `need` of them in all. Remove any such object and repeat. If
the process stops before the figure is exhausted, what remains is a set in
which every object meets more than its need within the set - and then
whichever of them is drawn last, in any order whatever, has one incidence too
many. That set is the certificate: it does not depend on the removal order,
it is printed below, and it can be checked against the labelled figure by
counting.

The same pass, allowing one object a single incidence beyond its need, lists
the objects at which the one unavoidable coincidence may sit.

    python3 core.py              # all three figures
    python3 core.py ring1        # one of straight, ring1, ring2

Needs only fits.py beside it, for the incidence tables.
"""
import importlib.util, os, sys

sp = importlib.util.spec_from_file_location(
    'fits', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fits.py'))
F = importlib.util.module_from_spec(sp)
sp.loader.exec_module(F)


def core(G, excused=()):
    N, INC, NEED = G['N'], G['INC'], G['NEED']
    alive = set(range(N))
    need = {i: NEED[i] + (1 if i in excused else 0) for i in range(N)}
    changed = True
    while changed:
        changed = False
        for v in list(alive):
            if sum(1 for u in INC[v] if u in alive) <= need[v]:
                alive.discard(v)
                changed = True
    return alive


def report(which):
    G = F.build(which, False, 0)
    obj, INC, NEED = G['obj'], G['INC'], G['NEED']
    K = core(G)
    name = lambda i: '%s%d' % obj[i]
    print('%s: %d objects' % (which, G['N']))
    if not K:
        print('   peels completely: an internal construction exists')
        return
    kinds = {k: sum(1 for i in K if obj[i][0] == k) for k in 'nst'}
    print('   the core has %d objects: %d nodes, %d circles, %d ring tips'
          % (len(K), kinds['n'], kinds['s'], kinds['t']))
    inside = sum(1 for v in K for u in INC[v] if u in K) // 2
    print('   %d incidences inside it, needs summing to %d' % (inside, sum(NEED[i] for i in K)))
    print('   every member meets more than its need within the core:')
    for i in sorted(K, key=lambda i: (obj[i][0], obj[i][1])):
        print('      %-5s meets %d, needs %d' % (name(i), sum(1 for u in INC[i] if u in K), NEED[i]))
    sites = [i for i in sorted(K) if not core(G, {i})]
    print('   excusing exactly one object clears it at %d sites: %s'
          % (len(sites), ' '.join(name(i) for i in sorted(sites, key=lambda i: (obj[i][0], obj[i][1])))))
    print()


if __name__ == '__main__':
    for which in (sys.argv[1:] or ('straight', 'ring1', 'ring2')):
        report(which)
