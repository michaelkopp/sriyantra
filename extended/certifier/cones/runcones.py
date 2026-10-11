"""Run the cone branch and bound over the surviving cones, resumable.
usage: python3 runcones.py ROOTS.json WORKER NWORKERS OUT.jsonl [BUDGET]
Worker w handles the cones with index = w mod NWORKERS; cones already recorded in OUT are skipped."""
import sys, os, json, time, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cone as C, bbcone as BB
roots = json.load(open(sys.argv[1])); w, nw = int(sys.argv[2]), int(sys.argv[3]); out = sys.argv[4]
budget = int(float(sys.argv[5])) if len(sys.argv) > 5 else 5000
done = set()
REDO = os.environ.get('REDO', '0') == '1'      # REDO=1: run again the cones whose last record has open boxes
if os.path.exists(out):
    last = {}
    for l in open(out):
        try: d = json.loads(l); last[d['i']] = d
        except Exception: pass
    done = {i for i, d in last.items() if not (REDO and d['open'])}
# SKIPFILES: other result files (e.g. another machine's); cones closed there are skipped.  REVERSE=1: run from the end.
for sf in [x for x in os.environ.get('SKIPFILES', '').split(',') if x]:
    if os.path.exists(sf) and sf.endswith('.json'):
        done.update(json.load(open(sf)))               # a list of closed cone indices
    elif os.path.exists(sf):
        for l in open(sf):
            try:
                d = json.loads(l)
                if not d['open']: done.add(d['i'])
            except Exception: pass
todo = [i for i in range(len(roots)) if i % nw == w and i not in done]
# CLAIMS=dir (used by run_all.sh with REDO=1): every worker takes from all cones that no result file shows closed;
# a cone is claimed by creating dir/i atomically, so the cones left open are shared evenly between the workers.
CLAIMS = os.environ.get('CLAIMS', '')
if CLAIMS:
    import glob
    closed = set()
    for sf in sorted(set(glob.glob('all_*.jsonl') + glob.glob('cloud_all_*.jsonl') + [x for x in os.environ.get('SKIPFILES', '').split(',') if x.endswith('.jsonl')])):
        for l in open(sf):
            try:
                d = json.loads(l)
                if not d['open']: closed.add(d['i'])
            except Exception: pass
    todo = [i for i in range(len(roots)) if i not in closed]
    os.makedirs(CLAIMS, exist_ok=True)
if os.environ.get('REVERSE', '0') == '1': todo.reverse()
f = open(out, 'a'); t0 = time.time(); last = t0; nopen = 0; nb = 0
for n, i in enumerate(todo):
    if CLAIMS:
        try: os.close(os.open(os.path.join(CLAIMS, str(i)), os.O_CREAT | os.O_EXCL | os.O_WRONLY))
        except FileExistsError: continue
        print('worker %d: cone %d' % (w, i), flush=True)
    T = C.Tree(roots[i]); st = {}; t1 = time.time()
    b, ob = BB.run_cone(T, budget, st); nb += b
    rec = {'i': i, 'boxes': b, 'open': len(ob), 'why': st, 'seconds': round(time.time() - t1, 2)}
    if ob:
        nopen += 1
        rec['parent'] = roots[i]
        rec['openboxes'] = [[[x[0] if math.isfinite(x[0]) else None, x[1]] for x in bx] for bx, _ in ob[:50]]
    f.write(json.dumps(rec) + '\n'); f.flush()
    if CLAIMS:
        print('worker %d: cone %d %s, %d boxes, %.0fs' % (w, i, 'OPEN (%d boxes left)' % len(ob) if ob else 'closed', b, rec['seconds']), flush=True)
    elif time.time() - last > 60:
        last = time.time(); el = last - t0
        print('worker %d: %d/%d cones in %.0fs (%.2fs/cone), %d boxes, %d cones with open boxes; eta %.0f min'
              % (w, n + 1, len(todo), el, el / (n + 1), nb, nopen, el / (n + 1) * (len(todo) - n - 1) / 60), flush=True)
print('worker %d finished: %d cones, %d with open boxes, %.0fs' % (w, len(todo), nopen, time.time() - t0), flush=True)
