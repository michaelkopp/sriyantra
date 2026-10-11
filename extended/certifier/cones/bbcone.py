"""Branch and bound inside the Cartesian-tree cones (cone.py).
usage: python3 bbcone.py roots.json i0 i1 budget_per_cone out.jsonl   (runs cones i0..i1-1 of the list)"""
import sys, os, time, json, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cone as C
INF = math.inf
U = float(os.environ.get('UMAX', '64')); MINW = float(os.environ.get('MINW', '1e-6'))
LPEVERY = os.environ.get('LPEVERY', '1') == '1'
NOREG = os.environ.get('NOREGIONS', '0') == '1'      # 1: do not use the regions certified in chart B
PICK = os.environ.get('PICK', '1') == '1'
USEMVF = os.environ.get('MVF', '1') == '1'           # mean-value relaxation (whole-space mode); tails allowed
MVPICK = os.environ.get('MVPICK', '1') == '1'
COMB = os.environ.get('COMB', '0') == '1'            # 1: one LP with the (s, x) rows and the MVF rows together
CONTRACT = int(os.environ.get('CONTRACT', '4'))      # rounds of LP contraction (branch and contract); 0: plain MVF test        # after a failed MVF, split the coordinate with the largest remainder share
import mvf as MV          # branch on the coordinate that loosens the correction rows most
def width(x):
    """splitting priority: the range of exp over the interval (the corrections depend on s only through exp(s);
    the linear parts are exact at any width)"""
    lo, hi = x
    return math.exp(hi) - (math.exp(lo) if lo > -INF else 0.0)
def split(x):
    lo, hi = x
    if lo == -INF: s = hi - max(2.0, abs(hi)); return [(s, hi), (-INF, s)]
    m = math.log(0.5 * (math.exp(lo) + math.exp(hi)))       # halves the range of exp
    m = min(max(m, lo + 1e-12 * (1 + abs(lo))), hi)
    return [(m, hi), (lo, m)]
def process(T, box, stats):
    """one box: None when it is closed, 'small' when it is too small to split, else its two halves"""
    E = C.Enc(T, box); mp = []
    cs = C.cells_hit(E)
    if not cs: stats['cell'] = stats.get('cell', 0) + 1; return None
    w = C.interval_refute(E)
    if w: stats['iv'] = stats.get('iv', 0) + 1; return None
    if not NOREG and C.in_D(E): stats['D'] = stats.get('D', 0) + 1; return None
    if not NOREG and C.in_RB(E): stats['RB'] = stats.get('RB', 0) + 1; return None
    ws = [width(x) for x in box]; k = max(range(len(box)), key=lambda i: ws[i])
    if LPEVERY or ws[k] < 0.1:
        stats['lpcalls'] = stats.get('lpcalls', 0) + 1
        if len(cs) > 1 and cs[0] is not None:
            # one LP over the box of all cells met; per-cell LPs only when it meets at most two
            if C.lp_refute_sx2(E, C.union_rows(cs)): stats['LP'] = stats.get('LP', 0) + 1; return None
            if len(cs) <= 2 and all(C.lp_refute_sx2(E, C.cell_rows(c)) for c in cs): stats['LP'] = stats.get('LP', 0) + 1; return None
        elif cs[0] is None and USEMVF and COMB:
            if MV.lp_refute_comb(E, T, box): stats['LP'] = stats.get('LP', 0) + 1; return None
        elif C.lp_refute_sx2(E, C.cell_rows(cs[0])): stats['LP'] = stats.get('LP', 0) + 1; return None
        if USEMVF and cs[0] is None:
            stats['mvfcalls'] = stats.get('mvfcalls', 0) + 1; mp = []
            if CONTRACT:
                # branch and contract: slabs of the box shown empty by exact certificates are cut off
                ref, box2, kp = MV.mvf_contract(T, box, CONTRACT, {}, stats)
                if ref: stats['MVF'] = stats.get('MVF', 0) + 1; return None
                if box2 != box:
                    stats['contracted'] = stats.get('contracted', 0) + 1
                    box = box2; E = C.Enc(T, box)
                    if C.interval_refute(E): stats['iv'] = stats.get('iv', 0) + 1; return None
                    ws = [width(x) for x in box]; k = max(range(len(box)), key=lambda i: ws[i])
                if MVPICK: mp = [(kp,)]
            elif MV.mvf_refute(T, box, {}, mp if MVPICK else None): stats['MVF'] = stats.get('MVF', 0) + 1; return None
    if ws[k] < MINW: return 'small', box
    if mp:                    # the MVF was tried: split where its remainder is largest
        k = mp[0][0]
    elif PICK:
        k2, p2 = pick_coord(E, box)
        if p2 > 1e-9: k = k2
    out = []
    for sb in split(box[k]):
        nbx = list(box); nbx[k] = sb; out.append(nbx)
    return out
def run_cone(T, budget, stats):
    stack = [([(-INF, 0.0)] * len(T.nodes), 0)]; nb = 0; openb = []
    t0 = time.time(); last = t0; HB = float(os.environ.get('HB', '0'))
    while stack:
        box, d = stack.pop(); nb += 1
        if HB and time.time() - last > HB:
            last = time.time()
            done = 1.0 - sum(2.0 ** -dd for _, dd in stack) - 2.0 ** -d     # share of the binary tree finished
            print('    cone heartbeat: %.0fs, %d boxes, stack %d, depth %d, open %d, done %.4f' % (last - t0, nb, len(stack), d, len(openb), done), file=sys.stderr, flush=True)
        if nb > budget: openb.append((box, 'budget')); continue
        r = process(T, box, stats)
        if r is None: continue
        if r[0] == 'small': openb.append((r[1], 'small')); continue
        for nbx in r: stack.append((nbx, d + 1))
    return nb, openb
if __name__ == '__main__':
    roots = json.load(open(sys.argv[1])); i0, i1 = int(sys.argv[2]), int(sys.argv[3]); budget = int(float(sys.argv[4]))
    out = open(sys.argv[5], 'a')
    for i in range(i0, min(i1, len(roots))):
        T = C.Tree(roots[i]); st = {}; t0 = time.time()
        nb, openb = run_cone(T, budget, st)
        rec = {'i': i, 'parent': roots[i], 'boxes': nb, 'open': len(openb), 'why': st, 'seconds': round(time.time() - t0, 2)}
        if openb: rec['openboxes'] = [[[x[0] if math.isfinite(x[0]) else None, x[1]] for x in b] for b, _ in openb[:20]]
        out.write(json.dumps(rec) + '\n'); out.flush()
        print(i, nb, len(openb), st, '%.1fs' % (time.time() - t0), flush=True)

# ---- branching on the coordinate that loosens the correction rows most
def _chordgap(a, b):
    """largest gap between the chord of exp over [a, b] and exp (a = -inf: exp(b))"""
    if a == -INF: return math.exp(b)
    if b - a < 1e-12: return 0.0
    s = (math.exp(b) - math.exp(a)) / (b - a); t = math.log(s)
    return math.exp(a) + s * (t - a) - s
def pick_coord(E, box):
    T = E.T; n = len(box)
    lw = [(x[1] - x[0]) if x[0] > -INF else None for x in box]
    pri = [0.0] * n
    for (a_, b_), (m, others) in T.pair.items():
        for _, w in others:
            wl, wh = C.lin_ival(w, box); wh = min(wh, 0.0)
            g = _chordgap(wl, wh)
            if g <= 0: continue
            idx = [i for i in range(n) if w[i] != 0]
            if any(lw[i] is None and box[i][1] > -60 for i in idx):
                for i in idx:
                    if lw[i] is None and box[i][1] > -60: pri[i] = max(pri[i], g)
                continue
            tot = sum(lw[i] for i in idx)
            if tot <= 0: continue
            for i in idx: pri[i] = max(pri[i], g * lw[i] / tot)
    k = max(range(n), key=lambda i: pri[i])
    return k, pri[k]
