"""Summarise cone runs: python3 summarize.py ROOTS.json OUT_*.jsonl"""
import sys, os, json, collections
roots = json.load(open(sys.argv[1])); seen = {}; why = collections.Counter(); boxes = 0; secs = 0
for f in sys.argv[2:]:
    if not os.path.exists(f): continue          # (an unmatched shell pattern)
    for l in open(f):
        d = json.loads(l)
        # a cone is closed when any of its records is closed (each closed record is a complete certificate)
        if d['i'] not in seen or seen[d['i']]['open']: seen[d['i']] = d
for d in seen.values():
    boxes += d['boxes']; secs += d['seconds']
    for k, v in d['why'].items(): why[k] += v
op = sorted(i for i, d in seen.items() if d['open'])
missing = [i for i in range(len(roots)) if i not in seen]
print('cones %d of %d done, %d missing, %d with open boxes; %d boxes, %.0f s CPU' % (len(seen), len(roots), len(missing), len(op), boxes, secs))
print('closing reasons:', {k: v for k, v in why.items() if k not in ('lpcalls', 'mvfcalls', 'cut', 'contracted')})
print('LP calls %d, MVF calls %d, slabs cut off by certified contraction %d' % (why['lpcalls'], why['mvfcalls'], why['cut']))
reg = {k: sorted(i for i, d in seen.items() if not d['open'] and d['why'].get(k)) for k in ('D', 'RB')}
for k, v in reg.items():
    if v: print('cones closed with the help of region %s: %d %s' % (k, len(v), v[:20] + (['...'] if len(v) > 20 else [])))
if op: print('cones with open boxes:', op[:50])
print('CLOSED' if not missing and not op else 'NOT YET CLOSED')
