"""Check that the records of the chart A runs close every cell.  python3 summarize.py

The 32,805 cells are the products chart x 9^4: in chart j the j-th of the five known gaps is the largest, and the
four others, in log-ratio to it, fall into nine bins each.  A cell is closed when one of its records in
results2-8.jsonl is closed with nothing delegated, or when it was split into sub-cells (ids CELL/i.j.k.l) that tile
its box and are all closed.  The six cells listed in SIX are closed by the cone certifier restricted to them
(../cones/six_*.jsonl, one record for each of the cones in ../cones/roots_surv.json); that is checked here too."""
import json, os, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
SIX = ['0.0.0.0.0', '0.0.0.0.1', '0.0.0.0.2', '0.0.0.1.2', '0.0.1.2.3', '0.1.2.3.4']
INF = float('inf')
def num(x): return -INF if x == '-inf' else INF if x == 'inf' else float(x)
recs = {}
for i in range(2, 9):
    for line in open(os.path.join(HERE, 'results%d.jsonl' % i)):
        if line.strip():
            r = json.loads(line); recs.setdefault(r['id'], []).append(r)
box = {i: [tuple(num(v) for v in iv) for iv in rs[0]['B']] for i, rs in recs.items()}
chart = {i: rs[0]['chart'] for i, rs in recs.items()}
top = sorted(i for i in recs if '/' not in i)
ids = set(tuple(int(x) for x in i.split('.')) for i in top)
assert ids == set(itertools.product(range(5), *[range(9)] * 4)), 'the cells are not the full grid'
# in every chart, the nine bins of each coordinate tile [-inf, 0]
for j in range(5):
    for m in range(4):
        bins = sorted({box[i][m] for i in top if chart[i] == j})
        assert len(bins) == 9 and bins[0][0] == -INF and bins[-1][1] == 0, (j, m)
        assert all(bins[k][1] == bins[k + 1][0] for k in range(8)), (j, m)
def tiles(parent, kids):
    """the kids' boxes are the 2^4 halves (or a finer product grid) of the parent box"""
    P = box[parent]
    for m in range(4):
        cuts = sorted({box[k][m] for k in kids})
        if cuts[0][0] != P[m][0] or cuts[-1][1] != P[m][1] or any(cuts[a][1] != cuts[a + 1][0] for a in range(len(cuts) - 1)):
            return False
    grid = set(tuple(box[k]) for k in kids)
    axes = [sorted({box[k][m] for k in kids}) for m in range(4)]
    return grid == set(itertools.product(*axes)) and len(kids) == len(grid)
def closed(i):
    if any(r['status'] == 'closed' and r.get('delegated', 0) == 0 for r in recs[i]): return True
    kids = [k for k in recs if k.startswith(i + '/') and '/' not in k[len(i) + 1:]]
    return bool(kids) and tiles(i, kids) and all(closed(k) for k in kids)
bysub = [i for i in top if not any(r['status'] == 'closed' and r.get('delegated', 0) == 0 for r in recs[i]) and closed(i)]
notc = [i for i in top if not closed(i)]
print('cells %d; closed by their own records %d, through sub-cells %d %s; not closed here %d %s'
      % (len(top), len(top) - len(bysub) - len(notc), len(bysub), bysub, len(notc), notc))
assert sorted(notc) == sorted(SIX), 'unexpected open cells'
# the six cells: every cone that survives the root screening has a closed record in ../cones/six_*.jsonl
C = os.path.join(HERE, '..', 'cones')
roots = json.load(open(os.path.join(C, 'roots_surv.json'))); seen = {}
for f in ('six_0.jsonl', 'six_1.jsonl'):
    for line in open(os.path.join(C, f)):
        d = json.loads(line)
        if d['i'] not in seen or seen[d['i']]['open']: seen[d['i']] = d
missing = [i for i in range(len(roots)) if i not in seen or seen[i]['open']]
print('six cells: %d cones, %d boxes, %d not closed' % (len(roots), sum(d['boxes'] for d in seen.values()), len(missing)))
print('CLOSED' if not missing else 'NOT CLOSED')
