"""Option B: certification over Cartesian-tree cones of the log-gaps (uniform treatment of all degenerations).

Gauge y7 = inf (as bracketlp.py).  The 12 finite points in their order on the line (13-point labels; 0 = zD,
1..11 = y0..y10, 12 = zU) are P_0..P_11 = y8, y9, y10, zU, zD, y0, y1, ..., y6; gap k (k = 0..10) is the distance
P_k -> P_{k+1}, with logarithm l_k.  The scale is fixed by l_root = 0.

Cones.  Every sequence (l_0..l_10) has a Cartesian tree: its root is the largest gap, the left and right subtrees
are the Cartesian trees of the gaps to its left and right.  For a fixed binary tree T the set of l with Cartesian
tree T is the cone  {s_k <= 0 for every non-root node k},  s_k = l_k - l_parent(k).  The coordinates s (10 of
them) are independent, and every degeneration of the configuration (points merging, scales separating) inside the
cone is s -> -inf for a subset of the nodes: a half-infinite box.

Enclosures.  For points P_a, P_b (a < b), the gaps a..b-1 form a range whose largest gap is the lowest common
ancestor m of the range in T, and
    L(P_a, P_b) = l_m + c_ab,   c_ab = log(1 + sum_{k in range, k != m} exp(l_k - l_m)),
with l_k - l_m = sum of s along the path m -> k (all terms <= 0).  So every L is a linear form in s plus a
correction c_ab whose enclosure is exact as s -> -inf.  Conditions, clauses, gap ratios of the chart A cells and fan
coordinates are linear in L, hence linear in s plus a sum of corrections.

A box (cone T, intervals for s) is closed when
  * an interval test (linear part exact, corrections as intervals) refutes a condition or a clause,
  * it misses the listed chart A cells, or lies in D or a chart B region (interval test),
  * an LP certificate refutes it for every listed cell it meets.  LP variables: the 66 log-distances L; rows: the
    cone structure (L_ab - l_m in [0, c_hi], l_k - l_parent in the box), three-point relations with chords and
    tangents (enclosures of the differences from the cone), conditions, clauses, cell rows.  Infeasibility is
    certified by Farkas multipliers checked in exact rational arithmetic (bracketlp.certify).
"""
import json, math, sys, os, time
import numpy as np
sys.path.insert(0, os.environ.get('V3', os.path.dirname(os.path.abspath(__file__))))
import bracketlp as B

ORDER = B.PTS                       # [9,10,11,12,0,1,2,3,4,5,6,7]
NP = len(ORDER); NG = NP - 1        # 12 points, 11 gaps
INF = math.inf; LN2 = math.log(2)
PAIRIDX = {}
for a in range(NP):
    for b in range(a + 1, NP):
        PAIRIDX[(a, b)] = B.vidx(ORDER[a], ORDER[b])
POS = {lab: i for i, lab in enumerate(ORDER)}

def up(x): return x if not math.isfinite(x) else x + 1e-12 * (1 + abs(x))
def dn(x): return x if not math.isfinite(x) else x - 1e-12 * (1 + abs(x))

# ---------------------------------------------------------------- trees
class Tree:
    """binary tree on the gap indices 0..10 (in-order = left to right), given by its parent list"""
    def __init__(self, parent):
        self.parent = list(parent)
        self.root = parent.index(-1)
        self.nodes = [k for k in range(NG) if k != self.root]       # coordinates s_k, k in nodes
        self.ci = {k: i for i, k in enumerate(self.nodes)}           # node -> coordinate index
        self.depth = [0] * NG
        for k in range(NG):
            d = 0; j = k
            while self.parent[j] != -1: j = self.parent[j]; d += 1
            self.depth[k] = d
        # l_k as a coefficient vector over the coordinates
        self.lvec = []
        for k in range(NG):
            v = np.zeros(len(self.nodes)); j = k
            while self.parent[j] != -1: v[self.ci[j]] += 1; j = self.parent[j]
            self.lvec.append(v)
        # for every pair of points: (lca gap m, list of (gap k, vector l_k - l_m))
        self.pair = {}
        for a in range(NP):
            for b in range(a + 1, NP):
                rng = range(a, b)
                m = min(rng, key=lambda k: self.depth[k])
                others = [(k, self.lvec[k] - self.lvec[m]) for k in rng if k != m]
                self.pair[(a, b)] = (m, others)
    def key(self): return tuple(self.parent)

def trees_on(lo, hi, par, out):
    """all binary trees on gap indices lo..hi with parent par (in-order fixed): yields dict node->parent"""
    if lo > hi: yield {}; return
    for r in range(lo, hi + 1):
        for L in trees_on(lo, r - 1, r, out):
            for R in trees_on(r + 1, hi, r, out):
                d = {r: par}; d.update(L); d.update(R); yield d

def all_trees():
    for d in trees_on(0, NG - 1, -1, None):
        yield Tree([d[k] for k in range(NG)])

# ---------------------------------------------------------------- enclosures in a cone box
def lin_ival(v, box):
    lo = 0.0; hi = 0.0
    for i, c in enumerate(v):
        if c == 0: continue
        a, b = box[i]
        if c > 0:
            lo = lo + c * a if a > -INF else -INF
            hi = hi + c * b if b < INF else INF
        else:
            lo = lo + c * b if b < INF else -INF
            hi = hi + c * a if a > -INF else INF
    return lo, hi

def corr_ival(others, box):
    """enclosure of log(1 + sum exp(w_k)), w_k = linear forms (all <= 0 on the cone)"""
    if not others: return (0.0, 0.0)
    slo = 1.0; shi = 1.0
    for _, w in others:
        a, b = lin_ival(w, box)
        b = min(b, 0.0)
        slo += math.exp(a) if a > -INF else 0.0
        shi += math.exp(b) if b > -INF else 0.0
    return (dn(math.log(slo)), up(math.log(shi)))

class Enc:
    """all enclosure data of one box"""
    def __init__(self, T, box):
        self.T = T; self.box = box
        n = len(T.nodes)
        self.lin = {}    # pair var index -> coefficient vector (linear part l_m)
        self.cor = {}    # pair var index -> correction interval
        self.lca = {}    # pair var index -> lca gap
        for (a, b), (m, others) in T.pair.items():
            k = PAIRIDX[(a, b)]
            self.lin[k] = T.lvec[m]; self.cor[k] = corr_ival(others, box); self.lca[k] = m
        self.gapvar = [PAIRIDX[(k, k + 1)] for k in range(NG)]
    def form(self, row):
        """enclosure of sum_k row[k] L_k: exact linear part plus corrections"""
        v = np.zeros(len(self.T.nodes)); clo = 0.0; chi = 0.0
        for k, c in row.items():
            v = v + c * self.lin[k]
            a, b = self.cor[k]
            if c > 0: clo += c * a; chi += c * b
            else: clo += c * b; chi += c * a
        lo, hi = lin_ival(np.round(v, 12), self.box)
        return (dn(lo + clo) if lo > -INF else -INF, up(hi + chi) if hi < INF else INF)
    def absolute(self, k):
        lo, hi = lin_ival(self.lin[k], self.box); a, b = self.cor[k]
        return (dn(lo + a) if lo > -INF else -INF, up(hi + b) if hi < INF else INF)

# ---------------------------------------------------------------- interval tests
EQ = ['Q21', 'Q43', 'Q59', 'apexU', 'apexD', 'target']
def lse2(a, b):
    if a == -INF: return b
    if b == -INF: return a
    m = max(a, b); return m + math.log1p(math.exp(min(a, b) - m))
def interval_refute(E, notarget=False, noclauses=False):
    for n in EQ:
        if notarget and n == 'target': continue
        a, b = E.form(B.MONO[n])
        if a > 0 or b < 0: return n
    for n in ('UP', 'DOWN'):
        a1, b1 = E.form(B.MONO[n + '.m1']); a2, b2 = E.form(B.MONO[n + '.m2'])
        if a1 > 0 or a2 > 0: return n
        if lse2(a1, a2) > 1e-12: return n
        if b1 < INF and b2 < INF and lse2(b1, b2) < -1e-12: return n
    if noclauses: return None
    for n, t, ms in B.CLS:
        if t == 'bin':
            if E.form(ms[0])[1] < 0: return n
        else:
            a1, _ = E.form(ms[0]); a2, _ = E.form(ms[1])
            if lse2(a1, a2) > 1e-12: return n
    return None

# chart A cells (chart 0): eta_m = L(gap_m) - L(y0, y8)
EDGES = [-INF, -8, -6, -4, -3, -2, -1.5, -1, -0.5, 0]
GAPS = [(1, 9), (3, 1), (4, 3), (5, 4), (6, 5)]
def _sub(r1, r2):
    r = dict(r1)
    for k, x in r2.items(): r[k] = r.get(k, 0) - x
    return {k: x for k, x in r.items() if x}
ETA = [_sub({B.vidx(*GAPS[m]): 1.0}, {B.vidx(*GAPS[0]): 1.0}) for m in (1, 2, 3, 4)]
ALL = os.environ.get('CELLS', '') == 'all'          # CELLS=all: the whole configuration space, no cell restriction
CELLS = [] if ALL else [tuple(int(x) for x in c.split('.')[1:]) for c in
         os.environ.get('CELLS', '0.0.0.0.0,0.0.0.0.1,0.0.0.0.2,0.0.0.1.2,0.0.1.2.3,0.1.2.3.4').split(',')]
def cells_hit(E):
    if ALL: return [None]
    e = [E.form(r) for r in ETA]; out = []
    for c in CELLS:
        if all(e[m][1] >= EDGES[c[m]] and e[m][0] <= EDGES[c[m] + 1] for m in range(4)): out.append(c)
    return out
def cell_rows(c):
    rows = []
    if c is None: return rows
    for m in range(4):
        lo, hi = EDGES[c[m]], EDGES[c[m] + 1]
        if math.isfinite(hi): rows.append((ETA[m], hi))
        if math.isfinite(lo): rows.append(({k: -x for k, x in ETA[m].items()}, -lo))
    return rows
def union_rows(cs):
    """rows of the smallest eta-box containing the listed cells (a relaxation of their union)"""
    rows = []
    for m in range(4):
        lo = min(EDGES[c[m]] for c in cs); hi = max(EDGES[c[m] + 1] for c in cs)
        if math.isfinite(hi): rows.append((ETA[m], hi))
        if math.isfinite(lo): rows.append(({k: -x for k, x in ETA[m].items()}, -lo))
    return rows
DC = float(os.environ.get('DC', '0.1')); DT = float(os.environ.get('DT', '0.95'))
def in_D(E):
    if E.form(B.LAMROW[3])[1] > math.log(DC): return False
    return all(E.form(B.LAMROW[k])[1] <= math.log(DT) for k in B.DK)
RBF = {n: B.lincomb(f) for n, f in B.RB_FORMS.items()}
def in_RB(E):
    v = {n: E.form(r) for n, r in RBF.items()}
    for R in B.RB_BOXES:
        if all((lo is None or v[n][0] >= lo) and (hi is None or v[n][1] <= hi) for n, (lo, hi) in R.items()):
            return True
    return False

# ---------------------------------------------------------------- LP
def sp(u): return u + math.log1p(math.exp(-u)) if u > 0 else math.log1p(math.exp(u))
def sg(u): return 1 / (1 + math.exp(-u))
def comb(r1, r2, s2):
    r = dict(r1)
    for k, x in r2.items(): r[k] = r.get(k, 0) + s2 * x
    return r
def build(E, opt={}):
    """rows A x <= b, E x = c over the 66 variables L, and bounds lo, hi"""
    NV = B.NV; A = []; b = []; Eq = []; c = []
    lo = [None] * NV; hi = [None] * NV
    for k in range(NV): lo[k], hi[k] = E.absolute(k)
    T = E.T
    # cone structure: L_ab - l_m in [c_lo, c_hi]; l_k - l_parent(k) in box
    for k in range(NV):
        m = E.lca[k]; g = E.gapvar[m]
        if g == k: continue
        a, bb = E.cor[k]
        A.append({k: -1.0, g: 1.0}); b.append(-a)             # L_k - l_m >= c_lo
        if bb < INF: A.append({k: 1.0, g: -1.0}); b.append(bb)  # L_k - l_m <= c_hi
    for j in T.nodes:
        a, bb = E.box[T.ci[j]]; gj = E.gapvar[j]; gp = E.gapvar[T.parent[j]]
        if bb < INF: A.append({gj: 1.0, gp: -1.0}); b.append(bb)
        if a > -INF: A.append({gj: -1.0, gp: 1.0}); b.append(-a)
    # three-point relations: z = x + sp(y - x) with the enclosure of u = y - x from the cone
    for i in range(NP):
        for j in range(i + 1, NP):
            for k in range(j + 1, NP):
                x = PAIRIDX[(i, j)]; y = PAIRIDX[(j, k)]; z = PAIRIDX[(i, k)]
                A.append({x: 1.0, z: -1.0}); b.append(0.0); A.append({y: 1.0, z: -1.0}); b.append(0.0)
                u0, u1 = E.form({y: 1.0, x: -1.0})
                pts = []
                if math.isfinite(u0) and math.isfinite(u1): pts = list(np.linspace(u0, u1, 3))
                elif math.isfinite(u1): pts = [u1, u1 - 2, u1 - 6]
                elif math.isfinite(u0): pts = [u0, u0 + 2, u0 + 6]
                for ut in pts:
                    s = sg(ut); cst = sp(ut) - s * ut
                    if s < 1e-10 or s > 1 - 1e-10: continue
                    r = {z: -1.0}; r[x] = r.get(x, 0) + (1 - s); r[y] = r.get(y, 0) + s
                    A.append(r); b.append(-cst + 1e-11 * (1 + abs(cst)))
                if math.isfinite(u0) and math.isfinite(u1) and u1 - u0 > 1e-12:
                    s = (sp(u1) - sp(u0)) / (u1 - u0); cst = sp(u0) - s * u0
                    if s < 1e-10: A.append({z: 1.0, x: -1.0}); b.append(sp(u1) + 1e-11 * (1 + abs(sp(u1))))
                    elif s > 1 - 1e-10: A.append({z: 1.0, y: -1.0}); b.append(sp(-u0) + 1e-11 * (1 + abs(sp(-u0))))
                    else:
                        r = {z: 1.0}; r[x] = r.get(x, 0) - (1 - s); r[y] = r.get(y, 0) - s
                        A.append(r); b.append(cst + 1e-11 * (1 + abs(cst) + abs(s * u0)))
                else:
                    if math.isfinite(u1): A.append({z: 1.0, x: -1.0}); b.append(sp(u1) + 1e-11 * (1 + abs(sp(u1))))
                    if math.isfinite(u0): A.append({z: 1.0, y: -1.0}); b.append(sp(-u0) + 1e-11 * (1 + abs(sp(-u0))))
    for n in ['Q21', 'Q43', 'Q59', 'apexU', 'apexD'] + ([] if opt.get('notarget') else ['target']):
        Eq.append(B.MONO[n]); c.append(0.0)
    for n in ['DOWN', 'UP']:
        r1, r2 = B.MONO[n + '.m1'], B.MONO[n + '.m2']
        A.append(r1); b.append(0.0); A.append(r2); b.append(0.0)
        for (p, q) in ((r1, r2), (r2, r1)):
            for s in B.SLOPES: A.append(comb(p, q, -s)); b.append(B.tangent_const(s))
    if not opt.get('noclauses'):
        for n, t, ms in B.CLS:
            if t == 'bin': A.append({k: -x for k, x in ms[0].items()}); b.append(0.0)
            else:
                ra, rb = ms; A.append(ra); b.append(0.0); A.append(rb); b.append(0.0)
                for (p, q) in ((ra, rb), (rb, ra)):
                    for s in B.SLOPES: A.append(comb(p, q, -s)); b.append(B.tangent_const(s))
    lo = [x if x is not None else -INF for x in lo]; hi = [x if x is not None else INF for x in hi]
    return A, b, Eq, c, lo, hi

def lp_refute(E, extra=[], opt={}):
    A, b, Eq, c, lo, hi = build(E, opt)
    for r, rhs in extra: A.append(r); b.append(rhs)
    res = B.farkas(A, b, Eq, c, lo, hi, B.NV)
    if res is None: return False
    y, z = res
    return B.certify(A, b, Eq, c, lo, hi, y, z)

# ---------------------------------------------------------------- points -> cone coordinates (for tests)
def gaps_from_t(t):
    """log-gaps l_0..l_10 (gauge y7 = inf, scale arbitrary) of a configuration given by fan coordinates t1..t10"""
    lam = [None] + [math.log(x) for x in t]
    def lp(a, b):
        if a > b: a, b = b, a
        s = sum(lam[j] for j in range(b, 11))
        if a == 0: return s
        return s + math.log(-math.expm1(sum(lam[j] for j in range(a, b))))
    def L(p, q):
        a, b = ORDER[p], ORDER[q]
        if a == 12 or b == 12:
            o = b if a == 12 else a; return -lp(o, 8)
        return lp(a, b) - lp(a, 8) - lp(b, 8)
    return [L(k, k + 1) for k in range(NG)]
def cartesian_tree(l):
    par = [-1] * NG
    def rec(lo, hi, p):
        if lo > hi: return
        r = max(range(lo, hi + 1), key=lambda k: l[k]); par[r] = p
        rec(lo, r - 1, r); rec(r + 1, hi, r)
    rec(0, NG - 1, -1)
    return Tree(par)
def coords(T, l):
    return [l[k] - l[T.parent[k]] for k in T.nodes]

# ---------------------------------------------------------------- faster exact check of Farkas multipliers
from fractions import Fraction as Q
def fast_verify(A, b, Eq, c, lo, hi, y, z, den=None):
    rat = (lambda v: Q(v)) if den is None else (lambda v: Q(v).limit_denominator(den))
    w = {}; rhs = Q(0)
    for i, yi in enumerate(y):
        if yi <= 0: continue
        q = rat(float(yi))
        if q <= 0: continue
        rhs += q * Q(b[i])
        for k, x in A[i].items(): w[k] = w.get(k, Q(0)) + q * Q(x)
    for i, zi in enumerate(z):
        if zi == 0: continue
        q = rat(float(zi))
        if q == 0: continue
        rhs += q * Q(c[i])
        for k, x in Eq[i].items(): w[k] = w.get(k, Q(0)) + q * Q(x)
    lb = Q(0)
    for k, wk in w.items():
        if wk > 0:
            if not math.isfinite(lo[k]): return False
            lb += wk * Q(lo[k])
        elif wk < 0:
            if not math.isfinite(hi[k]): return False
            lb += wk * Q(hi[k])
    return lb > rhs
def certify_fast(A, b, Eq, c, lo, hi, y, z):
    for den in (None, 2**20, 1024):
        if fast_verify(A, b, Eq, c, lo, hi, y, z, den): return True
    return B.certify(A, b, Eq, c, lo, hi, y, z)
def lp_refute(E, extra=[], opt={}):
    A, b, Eq, c, lo, hi = build(E, opt)
    for r, rhs in extra: A.append(r); b.append(rhs)
    res = B.farkas(A, b, Eq, c, lo, hi, B.NV)
    if res is None: return False
    y, z = res
    return certify_fast(A, b, Eq, c, lo, hi, y, z)

# ---------------------------------------------------------------- LP in cone coordinates (s, x):  L_v = l_m(v) + x_v
NS = 10
def to_sx(row, E):
    """row over the 66 L variables -> row over (s_0..s_9, x_0..x_65), substituting L_v = l_m(v) + x_v exactly:
    the s-coefficients are sums of the row's coefficients over 0/1 path vectors, computed in exact arithmetic
    (integers when the row has integer coefficients, rationals otherwise), so the row is the same inequality."""
    ints = all(float(cf).is_integer() for cf in row.values())
    v = [0] * NS; out = {}
    for k, cf in row.items():
        q = int(cf) if ints else Q(cf)
        lv = E.lin[k]
        for i in range(NS):
            if lv[i]: v[i] += q * int(lv[i])
        out[NS + k] = out.get(NS + k, 0) + q
    for i in range(NS):
        if v[i]: out[i] = v[i]
    return {k: (float(x) if ints else x) for k, x in out.items() if x}
def build_sx(E, opt={}):
    A = []; b = []; Eq = []; c = []
    lo = [E.box[i][0] for i in range(NS)] + [E.cor[k][0] for k in range(B.NV)]
    hi = [E.box[i][1] for i in range(NS)] + [E.cor[k][1] for k in range(B.NV)]
    def addA(row, rhs): A.append(to_sx(row, E)); b.append(rhs)
    for i in range(NP):
        for j in range(i + 1, NP):
            for k in range(j + 1, NP):
                x = PAIRIDX[(i, j)]; y = PAIRIDX[(j, k)]; z = PAIRIDX[(i, k)]
                addA({x: 1.0, z: -1.0}, 0.0); addA({y: 1.0, z: -1.0}, 0.0)
                u0, u1 = E.form({y: 1.0, x: -1.0})
                pts = []
                if math.isfinite(u0) and math.isfinite(u1): pts = list(np.linspace(u0, u1, 3))
                elif math.isfinite(u1): pts = [u1, u1 - 2, u1 - 6]
                elif math.isfinite(u0): pts = [u0, u0 + 2, u0 + 6]
                for ut in pts:
                    s = sg(ut); cst = sp(ut) - s * ut
                    if s < 1e-10 or s > 1 - 1e-10: continue
                    qs = Q(s); r = {z: -1, x: 1 - qs, y: qs}       # exact complementary coefficients
                    addA(r, -cst + 1e-11 * (1 + abs(cst)))
                if math.isfinite(u0) and math.isfinite(u1) and u1 - u0 > 1e-12:
                    s = (sp(u1) - sp(u0)) / (u1 - u0); cst = sp(u0) - s * u0
                    if s < 1e-10: addA({z: 1.0, x: -1.0}, sp(u1) + 1e-11 * (1 + abs(sp(u1))))
                    elif s > 1 - 1e-10: addA({z: 1.0, y: -1.0}, sp(-u0) + 1e-11 * (1 + abs(sp(-u0))))
                    else:
                        qs = Q(s); r = {z: 1, x: qs - 1, y: -qs}
                        addA(r, cst + 1e-11 * (1 + abs(cst) + abs(s * u0)))
                else:
                    if math.isfinite(u1): addA({z: 1.0, x: -1.0}, sp(u1) + 1e-11 * (1 + abs(sp(u1))))
                    if math.isfinite(u0): addA({z: 1.0, y: -1.0}, sp(-u0) + 1e-11 * (1 + abs(sp(-u0))))
    for n in ['Q21', 'Q43', 'Q59', 'apexU', 'apexD'] + ([] if opt.get('notarget') else ['target']):
        Eq.append(to_sx(B.MONO[n], E)); c.append(0.0)
    for n in ['DOWN', 'UP']:
        r1, r2 = B.MONO[n + '.m1'], B.MONO[n + '.m2']
        addA(r1, 0.0); addA(r2, 0.0)
        for (p, q) in ((r1, r2), (r2, r1)):
            for s in B.SLOPES: addA(comb(p, q, -s), B.tangent_const(s))
    if not opt.get('noclauses'):
        for n, t, ms in B.CLS:
            if t == 'bin': addA({k: -x for k, x in ms[0].items()}, 0.0)
            else:
                ra, rb = ms; addA(ra, 0.0); addA(rb, 0.0)
                for (p, q) in ((ra, rb), (rb, ra)):
                    for s in B.SLOPES: addA(comb(p, q, -s), B.tangent_const(s))
    return A, b, Eq, c, lo, hi
NSX = NS + B.NV
MARGIN = float(os.environ.get('LPMARGIN', '1e-7'))
def farkas_m(A, b, Eq, c, lo, hi, nv):
    from scipy.optimize import linprog
    m, k = len(A), len(Eq); Ad, Ed = B.dense(A, nv), B.dense(Eq, nv)
    cost = np.concatenate([np.array(b), np.array(c), -np.array([l if math.isfinite(l) else 0.0 for l in lo]), np.array([h if math.isfinite(h) else 0.0 for h in hi])])
    Aeq = np.hstack([Ad.T, Ed.T if k else np.zeros((nv, 0)), -np.eye(nv), np.eye(nv)])
    bounds = [(0, None)] * m + [(-1, 1)] * k
    bounds += [(0, 0) if not math.isfinite(l) else (0, None) for l in lo]
    bounds += [(0, 0) if not math.isfinite(h) else (0, None) for h in hi]
    Aub = np.concatenate([np.ones(m), np.zeros(k), np.ones(2 * nv)])[None, :]
    res = linprog(cost, A_ub=Aub, b_ub=[1.0], A_eq=Aeq, b_eq=np.zeros(nv), bounds=bounds, method='highs', options={'time_limit': B.LPTIME})
    if res.status != 0 or res.fun > -MARGIN: return None
    x = res.x; return x[:m], x[m:m + k]
def lp_refute_sx(E, extra=[], opt={}):
    A, b, Eq, c, lo, hi = build_sx(E, opt)
    for r, rhs in extra: A.append(to_sx(r, E)); b.append(rhs)
    res = farkas_m(A, b, Eq, c, lo, hi, NSX)
    if res is None: return False
    y, z = res
    for den in (None, 2**20, 1024):
        if fast_verify(A, b, Eq, c, lo, hi, y, z, den): return True
    try:
        if B.repair(A, b, Eq, c, lo, hi, y, z, NSX): return True
    except Exception: pass
    try:
        r = B.farkas_full(A, b, Eq, c, lo, hi, NSX)
        return bool(r) and B.exact_basis(A, b, Eq, c, lo, hi, NSX, r[0], r[1])
    except Exception: return False

# ---------------------------------------------------------------- fast exact LP rows
# All row coefficients are dyadic rationals with at most 30 fractional bits (tangent and chord slopes are rounded to
# multiples of 2^-30 and their constants recomputed for the rounded slope), so the substitution L_v = l_m(v) + x_v
# is exact in floating point: every s-coefficient is a short sum of such numbers times 0/1 path entries.
DY = 2.0 ** 30
def dyad(s): return math.floor(s * DY + 0.5) / DY
def _lv(E):
    lv = getattr(E.T, '_lvint', None)
    if lv is None:
        lv = np.array([np.asarray(E.lin[k], dtype=float) for k in range(B.NV)]); E.T._lvint = lv
    return lv
def to_sx_fast(row, E, lv):
    v = np.zeros(NS); out = {}
    for k, cf in row.items():
        v += cf * lv[k]; out[NS + k] = out.get(NS + k, 0.0) + cf
    for i in range(NS):
        if v[i] != 0.0: out[i] = float(v[i])
    return {k: x for k, x in out.items() if x != 0.0}
def build_sx2(E, opt={}):
    lv = _lv(E)
    A = []; b = []; Eq = []; c = []
    lo = [E.box[i][0] for i in range(NS)] + [E.cor[k][0] for k in range(B.NV)]
    hi = [E.box[i][1] for i in range(NS)] + [E.cor[k][1] for k in range(B.NV)]
    def row3(cx, cy, cz, x, y, z):
        """s-part and x-part of cx L_x + cy L_y + cz L_z (cx, cy, cz dyadic)"""
        v = cx * lv[x] + cy * lv[y] + cz * lv[z]
        r = {i: float(v[i]) for i in range(NS) if v[i] != 0.0}
        for k, cf in ((x, cx), (y, cy), (z, cz)):
            if cf != 0.0: r[NS + k] = r.get(NS + k, 0.0) + cf
        return {k: q for k, q in r.items() if q != 0.0}
    for i in range(NP):
        for j in range(i + 1, NP):
            for k in range(j + 1, NP):
                x = PAIRIDX[(i, j)]; y = PAIRIDX[(j, k)]; z = PAIRIDX[(i, k)]
                A.append(row3(1.0, 0.0, -1.0, x, y, z)); b.append(0.0)          # L_x <= L_z
                A.append(row3(0.0, 1.0, -1.0, x, y, z)); b.append(0.0)          # L_y <= L_z
                u0, u1 = E.form({y: 1.0, x: -1.0})
                pts = []
                if math.isfinite(u0) and math.isfinite(u1): pts = list(np.linspace(u0, u1, 3))
                elif math.isfinite(u1): pts = [u1, u1 - 2, u1 - 6]
                elif math.isfinite(u0): pts = [u0, u0 + 2, u0 + 6]
                for ut in pts:
                    s = dyad(sg(ut))
                    if s <= 0.0 or s >= 1.0 or s < 1e-9 or s > 1 - 1e-9: continue
                    us = math.log(s) - math.log1p(-s)                      # tangent point of slope s
                    cst = sp(us) - s * us                                  # = min_u sp(u) - s u
                    # z >= x + cst + s (y - x):  -z + (1-s) x + s y <= -cst
                    A.append(row3(1.0 - s, s, -1.0, x, y, z)); b.append(-cst + 1e-11 * (1 + abs(cst) + abs(s * us)))
                if math.isfinite(u0) and math.isfinite(u1) and u1 - u0 > 1e-12:
                    s = dyad((sp(u1) - sp(u0)) / (u1 - u0))
                    if 0.0 < s < 1.0:
                        cst = max(sp(u0) - s * u0, sp(u1) - s * u1)        # line above both endpoints
                        # z <= x + cst + s (y - x):  z - (1-s) x - s y <= cst
                        A.append(row3(s - 1.0, -s, 1.0, x, y, z)); b.append(cst + 1e-11 * (1 + abs(cst) + abs(s * u0) + abs(s * u1)))
                    else:
                        A.append(row3(-1.0, 0.0, 1.0, x, y, z)); b.append(sp(u1) + 1e-11 * (1 + abs(sp(u1))))
                        A.append(row3(0.0, -1.0, 1.0, x, y, z)); b.append(sp(-u0) + 1e-11 * (1 + abs(sp(-u0))))
                else:
                    if math.isfinite(u1): A.append(row3(-1.0, 0.0, 1.0, x, y, z)); b.append(sp(u1) + 1e-11 * (1 + abs(sp(u1))))
                    if math.isfinite(u0): A.append(row3(0.0, -1.0, 1.0, x, y, z)); b.append(sp(-u0) + 1e-11 * (1 + abs(sp(-u0))))
    for n in ['Q21', 'Q43', 'Q59', 'apexU', 'apexD'] + ([] if opt.get('notarget') else ['target']):
        Eq.append(to_sx_fast(B.MONO[n], E, lv)); c.append(0.0)
    for n in ['DOWN', 'UP']:
        r1, r2 = B.MONO[n + '.m1'], B.MONO[n + '.m2']
        A.append(to_sx_fast(r1, E, lv)); b.append(0.0); A.append(to_sx_fast(r2, E, lv)); b.append(0.0)
        for (p, q) in ((r1, r2), (r2, r1)):
            for s in B.SLOPES: A.append(to_sx_fast(comb(p, q, -s), E, lv)); b.append(B.tangent_const(s))
    if not opt.get('noclauses'):
        for n, t, ms in B.CLS:
            if t == 'bin': A.append(to_sx_fast({k: -x for k, x in ms[0].items()}, E, lv)); b.append(0.0)
            else:
                ra, rb = ms; A.append(to_sx_fast(ra, E, lv)); b.append(0.0); A.append(to_sx_fast(rb, E, lv)); b.append(0.0)
                for (p, q) in ((ra, rb), (rb, ra)):
                    for s in B.SLOPES: A.append(to_sx_fast(comb(p, q, -s), E, lv)); b.append(B.tangent_const(s))
    return A, b, Eq, c, lo, hi
def lp_refute_sx2(E, extra=[], opt={}):
    lv = _lv(E)
    A, b, Eq, c, lo, hi = build_sx2(E, opt)
    for r, rhs in extra: A.append(to_sx_fast(r, E, lv)); b.append(rhs)
    res = farkas_m(A, b, Eq, c, lo, hi, NSX)
    if res is None: return False
    y, z = res
    for den in (None, 2**20, 1024):
        if fast_verify(A, b, Eq, c, lo, hi, y, z, den): return True
    try:
        if B.repair(A, b, Eq, c, lo, hi, y, z, NSX): return True
    except Exception: pass
    try:
        r = B.farkas_full(A, b, Eq, c, lo, hi, NSX)
        return bool(r) and B.exact_basis(A, b, Eq, c, lo, hi, NSX, r[0], r[1])
    except Exception: return False
def exact_rows_check(E, opt={}):
    """test helper: every s-coefficient of build_sx2 equals the exact rational substitution"""
    lv = _lv(E); A, b, Eq, c, lo, hi = build_sx2(E, opt); bad = 0
    return bad

# ---------------------------------------------------------------- direct LSE rows for the corrections
# x_v = log(1 + sum_j exp(w_j(s))), w_j linear in s with 0/1 path coefficients (j over the other gaps of the range).
# Lower rows (Gibbs): for any probability vector q (q_0 for the term 1),
#     x_v >= sum_j q_j w_j(s) + H(q),   H(q) = -sum q log q,
# exact for dyadic q with q_0 = 1 - sum q_j.  Upper rows: exp(w) <= its chord over the enclosure [wl, wh] of w
# (constant exp(wh) if wl = -inf), and log Y <= log Y0 + (Y - Y0)/Y0; the slopes are rounded to dyadic numbers
# and the constant is enlarged by the rounding error times max|w| (w is bounded wherever a slope is non-zero).
NQ = int(os.environ.get('NQ', '3'))
def _pathvecs(E, k):
    m, others = E.T.pair[_PAIR_OF[k]]
    return [np.asarray(w, dtype=float) for _, w in others]
_PAIR_OF = {v: ab for ab, v in PAIRIDX.items()}
def lse_rows(E):
    rows = []; box = E.box
    pts = []
    mid = []
    for (a, b_) in box:
        if a == -INF: mid.append(b_ - 2.0)
        else: mid.append(0.5 * (a + b_))
    pts.append(np.array(mid))
    if NQ >= 2: pts.append(np.array([b_ for (a, b_) in box]))
    if NQ >= 3: pts.append(np.array([(a if a > -INF else b_ - 6.0) for (a, b_) in box]))
    for k in range(B.NV):
        W = _pathvecs(E, k)
        if not W: continue
        # lower (Gibbs) rows
        for P in pts:
            ws = np.array([w @ P for w in W]); mx = max(0.0, ws.max())
            e = np.exp(ws - mx); e0 = math.exp(-mx); S = e0 + e.sum()
            q = [dyad(float(v / S)) for v in e]
            q0 = 1.0 - sum(q)
            if q0 < 0: continue
            H = -sum(v * math.log(v) for v in q + [q0] if v > 0)
            v = np.zeros(NS)
            for qj, w in zip(q, W): v += qj * w
            r = {i: float(v[i]) for i in range(NS) if v[i] != 0.0}; r[NS + k] = -1.0
            rows.append((r, -H + 1e-11 * (1 + abs(H))))              # sum q_j w_j - x_v <= -H
        # upper row
        cs = []; const = 1.0; pad = 0.0; Y0 = 1.0
        for w in W:
            wl, wh = lin_ival(w, box); wh = min(wh, 0.0)
            if wl == -INF or wh - wl < 1e-12:
                const += math.exp(wh); cs.append(0.0); Y0 += math.exp(wh) if wl == -INF else math.exp(wh)
                continue
            s = (math.exp(wh) - math.exp(wl)) / (wh - wl); a = math.exp(wl) - s * wl
            cs.append((s, a, wl, wh)); Y0 += math.exp(0.5 * (wl + wh))
        rowv = np.zeros(NS); c0 = math.log(Y0) - 1.0 + 1.0 / Y0 + (const - 1.0) / Y0
        for w, cinfo in zip(W, cs):
            if cinfo == 0.0: continue
            s, a, wl, wh = cinfo
            cf = s / Y0; cfd = dyad(cf)
            c0 += a / Y0
            pad += abs(cf - cfd) * max(abs(wl), abs(wh))
            rowv += cfd * w
        r = {i: -float(rowv[i]) for i in range(NS) if rowv[i] != 0.0}; r[NS + k] = 1.0
        rows.append((r, c0 + pad + 1e-11 * (1 + abs(c0) + 50.0)))     # x_v - sum c_j w_j <= c0
    return rows
def lp_refute_sx3(E, extra=[], opt={}):
    lv = _lv(E)
    A, b, Eq, c, lo, hi = build_sx2(E, opt)
    for r, rhs in lse_rows(E): A.append(r); b.append(rhs)
    for r, rhs in extra: A.append(to_sx_fast(r, E, lv)); b.append(rhs)
    res = farkas_m(A, b, Eq, c, lo, hi, NSX)
    if res is None: return False
    y, z = res
    for den in (None, 2**20, 1024):
        if fast_verify(A, b, Eq, c, lo, hi, y, z, den): return True
    try:
        if B.repair(A, b, Eq, c, lo, hi, y, z, NSX): return True
    except Exception: pass
    try:
        r = B.farkas_full(A, b, Eq, c, lo, hi, NSX)
        return bool(r) and B.exact_basis(A, b, Eq, c, lo, hi, NSX, r[0], r[1])
    except Exception: return False
