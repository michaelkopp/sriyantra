"""Mean-value (first-order Taylor) relaxation in the cone coordinates s; boxes may have half-infinite coordinates.

Every condition and clause f is a smooth function of the ten cone coordinates: the log-distances are
L_v(s) = l_m(v)(s) + x_v,  x_v = log(1 + sum_j exp(w_j(s))),  l_m and w_j linear with 0/1 path coefficients.
For a half-infinite coordinate s_i in (-inf, a] write u_i = exp(s_i) in [0, e^a]: exp(w_j(s)) is a product of
exp(s_F-part) and of u_i (for every tail coordinate on the path), so x_v is smooth in z = (s_F, u) up to u = 0.
Write f = k.s_T + h(z), k the (integer) coefficients of the tail coordinates in the linear part of f.  On the box
(centre z_c, half-widths r) and for any vector g0 (dyadic, hence exact)
    h(z) - g0.z  in  [h(z_c) - g0.z_c - R, h(z_c) - g0.z_c + R],   R = sum_i max|G_i - g0_i| r_i,
G the interval gradient of h over the box.  Conditions f = 0 give two rows, clauses f >= 0 one row, in the variables
(s, u); clauses of the form -log(e^ma + e^mb) whose two terms have different tail coefficients are omitted (a
weaker system).  The LP is infeasible only if no point of the box satisfies the relaxed system; its Farkas
multipliers are checked in exact rational arithmetic.  h(z_c) and R are enlarged by a relative 1e-12."""
import math, numpy as np
import cone as C, bracketlp as B
INF = math.inf
sp = C.sp

def dy(v): return math.floor(v * 2.0 ** 30 + 0.5) / 2.0 ** 30
UW = 1.0

class MVF:
    def __init__(self, T, box):
        self.T = T; self.box = box; n = len(box); self.contrib = None
        # u-coordinates: the half-infinite ones and the wide ones (width > UW), where exp(s) varies too much
        self.TL = [i for i in range(n) if box[i][0] == -INF or box[i][1] - box[i][0] > UW]
        self.FN = [i for i in range(n) if i not in self.TL]
        self.ulo = np.array([(math.exp(box[i][0]) * (1 - 1e-14) if box[i][0] > -INF else 0.0) for i in self.TL])
        self.uhi = np.array([math.exp(box[i][1]) * (1 + 1e-14) for i in self.TL])
        self.uc = 0.5 * (self.ulo + self.uhi); self.rU = 0.5 * (self.uhi - self.ulo)
        self.sc = np.array([0.5 * (a + b) for a, b in box])
        for t, i in enumerate(self.TL): self.sc[i] = math.log(self.uc[t])
        self.rF = np.array([0.5 * (box[i][1] - box[i][0]) for i in self.FN])
        nF, nT = len(self.FN), len(self.TL)
        self.px = {}
        for (a_, b_), (m, others) in T.pair.items():
            k = C.PAIRIDX[(a_, b_)]
            W = [np.asarray(w, dtype=float) for _, w in others]
            assert all(np.all((w == 0) | (w == 1)) for w in W)        # path vectors
            if not W:
                self.px[k] = (0.0, np.zeros(nF), np.zeros(nT), np.zeros(nF), np.zeros(nF), np.zeros(nT), np.zeros(nT), 0.0); continue
            wc = np.array([w @ self.sc for w in W]); mx = max(0.0, wc.max())
            e = np.exp(wc - mx); S = math.exp(-mx) + e.sum(); x = mx + math.log(S); p = e / S
            gF = np.zeros(nF); gU = np.zeros(nT)
            for pj, w in zip(p, W):
                gF += pj * w[self.FN]
            for t, i in enumerate(self.TL):            # d x / d u_i = sum_{j: i on path} exp(w_j - s_i) / S
                gU[t] = sum(pj for pj, w in zip(p, W) if w[i]) / math.exp(self.sc[i])
            # intervals over the box
            wl = []; wh = []
            for w in W:
                a1, b1 = C.lin_ival(w, box); wl.append(a1); wh.append(min(b1, 0.0))
            el = [math.exp(v) if v > -INF else 0.0 for v in wl]; eh = [math.exp(v) for v in wh]
            Slo = 1 + sum(el); Shi = 1 + sum(eh)
            # d x / d s_i = E_J / (1 + E_J + E_N), J the paths through i (0/1 coefficients), E = sum of exp(w_j):
            # increasing in E_J and decreasing in E_N, so the enclosures of E_J, E_N bound it
            sEl = sum(el); sEh = sum(eh)
            gFl = np.zeros(nF); gFh = np.zeros(nF)
            for t, i in enumerate(self.FN):
                J = [j for j, w in enumerate(W) if w[i] != 0]
                if not J: continue
                EJl = sum(el[j] for j in J); EJh = sum(eh[j] for j in J); ENl = max(0.0, sEl - EJl); ENh = max(0.0, sEh - EJh)
                gFl[t] = EJl / (1 + EJl + ENh) * (1 - 1e-12); gFh[t] = min(1.0, EJh / (1 + EJh + ENl) * (1 + 1e-12))
            # d x / d u_i = F_J / (1 + u_i F_J + E_N), F_J = sum over J of exp(w_j - s_i): increasing in F_J,
            # decreasing in u_i and E_N
            gUl = np.zeros(nT); gUh = np.zeros(nT)
            for t, i in enumerate(self.TL):
                J = [j for j, w in enumerate(W) if w[i] != 0]
                if not J: continue
                FJl = 0.0; FJh = 0.0
                for j in J:
                    w2 = W[j].copy(); w2[i] = 0.0
                    a2, b2 = C.lin_ival(w2, box)
                    FJl += math.exp(a2) if a2 > -INF else 0.0; FJh += math.exp(min(b2, 0.0))
                ENl = max(0.0, sEl - sum(el[j] for j in J)); ENh = max(0.0, sEh - sum(eh[j] for j in J))
                gUl[t] = FJl / (1 + self.uhi[t] * FJl + ENh) * (1 - 1e-12)
                gUh[t] = FJh / (1 + self.ulo[t] * FJh + ENl) * (1 + 1e-12)
            self.px[k] = (x, gF, gU, gFl, gFh, gUl, gUh, abs(x))
    def lin(self, row):
        """for f = sum row[v] L_v: tail coefficients k, and h(z_c), grad h, interval of grad h (in z = (s_F, u))"""
        nF, nT = len(self.FN), len(self.TL)
        linv = np.zeros(len(self.box)); h = 0.0; gF = np.zeros(nF); gU = np.zeros(nT)
        gFl = np.zeros(nF); gFh = np.zeros(nF); gUl = np.zeros(nT); gUh = np.zeros(nT); av = 0.0
        for k, cf in row.items():
            m = C._PAIR_OF and self.T.pair[C._PAIR_OF[k]][0]
            linv += cf * np.asarray(self.T.lvec[m], dtype=float)
            x, a, b, c, d, e, f, ax = self.px[k]
            h += cf * x; av += abs(cf) * (ax + 1)
            gF += cf * a; gU += cf * b
            if cf > 0: gFl += cf * c; gFh += cf * d; gUl += cf * e; gUh += cf * f
            else: gFl += cf * d; gFh += cf * c; gUl += cf * f; gUh += cf * e
        kT = np.array([linv[i] for i in self.TL]); lF = np.array([linv[i] for i in self.FN])
        h += float(lF @ self.sc[self.FN]); av += float(np.abs(lF) @ np.abs(self.sc[self.FN]))
        return kT, h, gF + lF, gU, gFl + lF, gFh + lF, gUl, gUh, av
    def lse(self, ra, rb, sign):
        """f = sign * log(e^ma + e^mb); None when ma, mb have different tail coefficients"""
        A_ = self.lin(ra); B_ = self.lin(rb)
        if np.any(A_[0] != B_[0]): return None
        kT = A_[0]; ha, hb = A_[1], B_[1]
        m = max(ha, hb); h = m + math.log(math.exp(ha - m) + math.exp(hb - m))
        qa = math.exp(ha - h); qb = math.exp(hb - h)
        E = C.Enc(self.T, self.box)
        dl, dh = E.form(C._sub(rb, ra))              # enclosure of mb - ma (tail parts cancel)
        qal = 1 / (1 + math.exp(dh)) if dh < 700 else 0.0
        qah = 1 / (1 + math.exp(dl)) if dl > -700 else 1.0
        qbl, qbh = 1 - qah, 1 - qal
        def imul(lo1, hi1, lo2, hi2):
            c = [lo1 * lo2, lo1 * hi2, hi1 * lo2, hi1 * hi2]; return np.minimum.reduce(c), np.maximum.reduce(c)
        gF = qa * A_[2] + qb * B_[2]; gU = qa * A_[3] + qb * B_[3]
        l1, h1 = imul(qal, qah, A_[4], A_[5]); l2, h2 = imul(qbl, qbh, B_[4], B_[5])
        l3, h3 = imul(qal, qah, A_[6], A_[7]); l4, h4 = imul(qbl, qbh, B_[6], B_[7])
        out = [kT, h, gF, gU, l1 + l2, h1 + h2, l3 + l4, h3 + h4, A_[8] + B_[8] + abs(h)]
        if sign < 0:
            out = [-kT, -h, -gF, -gU, -(h1 + h2), -(l1 + l2), -(h3 + h4), -(l3 + l4), out[8]]
        return out
    def row(self, kT, h, gF, gU, gFl, gFh, gUl, gUh, av):
        """coefficients (s_F, s_T, u) and [lo, hi] with  f - coef.(s, u) in [lo, hi] on the box"""
        g0F = np.array([dy(v) for v in gF]); g0U = np.array([dy(v) for v in gU])
        cF = np.maximum(np.abs(gFl - g0F), np.abs(gFh - g0F)) * self.rF if len(g0F) else np.zeros(0)
        R = float(np.sum(cF))
        cU = np.zeros(len(g0U))
        if len(g0U):
            if not np.all(np.isfinite(gUh)): return None
            cU = np.maximum(np.abs(gUl - g0U), np.abs(gUh - g0U)) * self.rU
            R += float(np.sum(cU))
        if self.contrib is not None:            # share of each coordinate in the remainder R (branching)
            for t, i in enumerate(self.FN): self.contrib[i] += cF[t]
            for t, i in enumerate(self.TL): self.contrib[i] += cU[t]
        D = h - float(g0F @ self.sc[self.FN]) - (float(g0U @ self.uc) if len(g0U) else 0.0)
        pad = 1e-12 * (1 + av + R + float(np.abs(g0F) @ np.abs(self.sc[self.FN])))
        coef = {}
        for t, i in enumerate(self.FN):
            if g0F[t] != 0.0: coef[i] = float(g0F[t])
        for t, i in enumerate(self.TL):
            if kT[t] != 0.0: coef[i] = float(kT[t])
        n = len(self.box)
        for t in range(len(self.TL)):
            if g0U[t] != 0.0: coef[n + t] = float(g0U[t])
        return coef, D - R - pad, D + R + pad

def _system(T, box, opt={}, M=None):
    """rows A z <= b of the relaxation over z = (s_0..s_{n-1}, u_t for the tail coordinates) and the bounds"""
    if M is None: M = MVF(T, box)
    n = len(box); nT = len(M.TL)
    A = []; b = []
    def add(r, sense):
        if r is None: return
        coef, lo, hi = r
        if sense == 'eq':
            A.append(dict(coef)); b.append(-lo); A.append({k: -v for k, v in coef.items()}); b.append(hi)
        else:
            A.append({k: -v for k, v in coef.items()}); b.append(hi)
    for nme in ['Q21', 'Q43', 'Q59', 'apexU', 'apexD'] + ([] if opt.get('notarget') else ['target']):
        add(M.row(*M.lin(B.MONO[nme])), 'eq')
    # extra variables (delta, w) for the trinomials whose two terms differ on the u-coordinates
    xlo = []; xhi = []; E = None
    def lse_split(ra, rb, sense):
        """log(e^ma + e^mb) = ma + sp(delta), delta = mb - ma: an MVF row for ma, an MVF row for delta (as a new
        variable with its enclosure as bounds) and w = sp(delta) relaxed by tangents and a chord of sp"""
        nonlocal E
        if E is None: E = C.Enc(T, box)
        dl, dh = E.form(C._sub(rb, ra))
        if dh == INF and dl > -INF: ra, rb = rb, ra; dl, dh = -dh, -dl
        r_a = M.row(*M.lin(ra)); r_d = M.row(*M.lin(C._sub(rb, ra)))
        if r_a is None or r_d is None: return
        jd = n + nT + len(xlo); jw = jd + 1
        xlo.append(dl); xhi.append(dh)
        wl = 0.0 if dl == -INF else max(0.0, sp(dl) * (1 - 1e-12) - 1e-300)
        wh = INF if dh == INF else sp(dh) * (1 + 1e-12) + 1e-300
        xlo.append(wl); xhi.append(wh)
        coef, lo_, hi_ = r_d                    # delta - coef.z in [lo_, hi_]
        A.append({**coef, jd: -1.0}); b.append(-lo_)
        A.append({**{k: -v for k, v in coef.items()}, jd: 1.0}); b.append(hi_)
        if dl > -INF and dh < INF: pts = [dl + f * (dh - dl) for f in (0.0, 0.25, 0.5, 0.75, 1.0)]
        elif dh < INF: pts = [dh - d for d in (0.0, 1.0, 2.0, 4.0, 8.0)]
        else: pts = [-4.0, -1.0, 0.0, 1.0, 4.0]
        for t0 in pts:                          # w >= sg delta + cst (tangent of sp with dyadic slope sg)
            sg_ = dy(1 / (1 + math.exp(-t0)))
            if not 0 < sg_ < 1: continue
            us = math.log(sg_) - math.log1p(-sg_); cst = sp(us) - sg_ * us
            cst -= 1e-12 * (1 + abs(cst) + abs(sg_ * us))
            A.append({jd: sg_, jw: -1.0}); b.append(-cst)
        if dl > -INF and dh < INF and dh - dl > 1e-9:     # chord: w - sc delta <= c
            sc = dy((sp(dh) - sp(dl)) / (dh - dl))
            c = max(sp(dl) - sc * dl, sp(dh) - sc * dh); c += 1e-12 * (1 + abs(c) + abs(sc) * max(abs(dl), abs(dh)))
            A.append({jw: 1.0, jd: -sc}); b.append(c)
        coef, lo_, hi_ = r_a                    # ma - coef.z in [lo_, hi_];  f = ma + w
        A.append({**coef, jw: 1.0}); b.append(-lo_)            # eq: 0 = ma + w;  ge: ma + w <= 0
        if sense == 'eq': A.append({**{k: -v for k, v in coef.items()}, jw: -1.0}); b.append(hi_)
    for nme in ('UP', 'DOWN'):
        o = M.lse(B.MONO[nme + '.m1'], B.MONO[nme + '.m2'], +1)
        if o is not None: add(M.row(*o), 'eq')
        elif not opt.get('nosplit'): lse_split(B.MONO[nme + '.m1'], B.MONO[nme + '.m2'], 'eq')
    if not opt.get('noclauses'):
        for nme, t, ms in B.CLS:
            if t == 'bin': add(M.row(*M.lin(ms[0])), 'ge')
            else:
                o = M.lse(ms[0], ms[1], -1)
                if o is not None: add(M.row(*o), 'ge')
                elif not opt.get('nosplit'): lse_split(ms[0], ms[1], 'ge')
    # u_t = exp(s_i) >= tangent lines of exp (dyadic slope al, the line through (log al, al) rounded down)
    for t, i in enumerate(M.TL):
        a_, b_ = box[i]
        pts = (b_, b_ - 1.0, b_ - 2.0, b_ - 4.0) if a_ == -INF else (b_, b_ - 1.0, 0.5 * (a_ + b_), a_)
        for c0 in pts:
            al = dy(math.exp(c0))
            if al <= 0: continue
            cst = al * (1 - math.log(al)); cst -= 1e-12 * (1 + abs(cst)) + 1e-300
            A.append({i: al, n + t: -1.0}); b.append(-cst)        # al s - u <= -cst  (tangent at log al)
        if a_ > -INF:                                              # chord: u - sg s <= c on [a, b]
            sg = dy((math.exp(b_) - math.exp(a_)) / (b_ - a_))
            c = max(math.exp(a_) - sg * a_, math.exp(b_) - sg * b_); c += 1e-12 * (1 + abs(c) + abs(sg * a_))
            A.append({i: -sg, n + t: 1.0}); b.append(c)
    lo = [a for a, _ in box] + list(M.ulo) + xlo; hi = [bb for _, bb in box] + list(M.uhi) + xhi
    return M, A, b, lo, hi

def _certified_empty(A, b, lo, hi, nv):
    res = B.farkas(A, b, [], [], lo, hi, nv)
    if res is None: return False
    y, z = res
    for den in (None, 2**20, 1024):
        if C.fast_verify(A, b, [], [], lo, hi, y, z, den): return True
    return False

def mvf_contract(T, box, rounds=4, opt={}, stats=None):
    """Branch-and-contract step.  Returns (refuted, newbox, pick):  refuted = the box holds no point of the system;
    otherwise every point of the box that satisfies the system lies in newbox (each slab cut off the box is shown
    empty by its own exact Farkas certificate, over the relaxation of the box it was cut from).  pick: the coordinate
    with the largest remainder share in the last relaxation."""
    from scipy.optimize import linprog
    box = list(box); n = len(box); pick = None
    for rd in range(rounds):
        mp = []
        M = MVF(T, box); M.contrib = np.zeros(n)
        M, A, b, lo, hi = _system(T, box, opt, M)
        pick = int(np.argmax(M.contrib))
        nT = len(M.TL); nv = len(lo)
        if _certified_empty(A, b, lo, hi, nv): return True, box, pick
        Ad = B.dense(A, nv); bb = np.array(b)
        lo2 = list(lo); hi2 = list(hi); changed = False
        for j in list(M.FN) + [n + t for t in range(nT)] + list(M.TL):
            for sgn in (1.0, -1.0):          # min x_j, max x_j
                c = np.zeros(nv); c[j] = sgn
                bnds = [(l if math.isfinite(l) else None, h if math.isfinite(h) else None) for l, h in zip(lo2, hi2)]
                r = linprog(c, A_ub=Ad, b_ub=bb, bounds=bnds, method='highs')
                if r.status == 2:                       # infeasible: certify the whole box
                    if _certified_empty(A, b, lo2, hi2, nv): return True, box, pick
                    break
                if r.status != 0: continue
                v = sgn * r.fun; L, H = lo2[j], hi2[j]
                wj = (H - L) if math.isfinite(L) else (abs(H) + 1.0)
                mg = 1e-6 * (1 + abs(v)) + 1e-3 * wj
                if sgn > 0:
                    beta = v - mg
                    if not (beta > L + 0.02 * wj or (not math.isfinite(L) and beta > -1e300)): continue
                    if j >= n and beta <= 0: continue
                    h3 = list(hi2); h3[j] = beta
                    if _certified_empty(A, b, lo2, h3, nv):
                        lo2[j] = beta; changed = True
                        if stats is not None: stats['cut'] = stats.get('cut', 0) + 1
                else:
                    beta = v + mg
                    if not (beta < H - 0.02 * wj): continue
                    l3 = list(lo2); l3[j] = beta
                    if _certified_empty(A, b, l3, hi2, nv):
                        hi2[j] = beta; changed = True
                        if stats is not None: stats['cut'] = stats.get('cut', 0) + 1
        if not changed: break
        # the new box: s-bounds, and for tail coordinates the bounds implied by those on u = exp(s)
        nb = []
        for i in range(n):
            a, b_ = lo2[i], hi2[i]
            if i in M.TL:
                t = M.TL.index(i)
                ul, uh = lo2[n + t], hi2[n + t]
                if ul > 0: a = max(a, math.log(ul) - 1e-12 * (1 + abs(math.log(ul))))   # u >= ul  =>  s >= log ul
                if uh < math.exp(b_): b_ = min(b_, math.log(uh) + 1e-12 * (1 + abs(math.log(uh))))
            nb.append((a, b_))
        if any(a > b_ for a, b_ in nb): return True, box, pick     # the certified slabs cover the whole box
        shrink = max((1 - (nb[i][1] - nb[i][0]) / max(box[i][1] - box[i][0], 1e-300)) if math.isfinite(box[i][0]) else (1.0 if math.isfinite(nb[i][0]) else 0.0) for i in range(n))
        box = nb
        if shrink < 0.1: break
    return False, box, pick

def mvf_refute(T, box, opt={}, pick=None):
    """True when the relaxed system is infeasible on the box (exact certificate).  If pick is a list, the
    coordinate with the largest total share in the remainders R of the rows is appended to it."""
    M = MVF(T, box); n = len(box)
    if pick is not None: M.contrib = np.zeros(n)
    M, A, b, lo, hi = _system(T, box, opt, M)
    if pick is not None: pick.append((int(np.argmax(M.contrib)), float(np.max(M.contrib)), M.contrib))
    return _certified_empty(A, b, lo, hi, len(lo))

def lp_refute_comb(E, T, box, opt={}):
    """One LP with the (s, x) rows of cone.build_sx2 and the rows of the mean-value relaxation together (both are
    valid for every point of the system in the box, with x the corrections, u = exp(s_U) and delta, w as in
    lse_split); its infeasibility is certified as in cone.lp_refute_sx2."""
    A, b, Eq, c, lo, hi = C.build_sx2(E, opt)
    n = len(box); NSX = len(lo)
    M, A2, b2, lo2, hi2 = _system(T, box, opt)
    sh = NSX - n
    for r, rhs in zip(A2, b2):
        A.append({(k if k < n else k + sh): v for k, v in r.items()}); b.append(rhs)
    lo = lo + lo2[n:]; hi = hi + hi2[n:]; nv = len(lo)
    res = C.farkas_m(A, b, Eq, c, lo, hi, nv)
    if res is None: return False
    y, z = res
    for den in (None, 2**20, 1024):
        if C.fast_verify(A, b, Eq, c, lo, hi, y, z, den): return True
    try:
        if B.repair(A, b, Eq, c, lo, hi, y, z, nv): return True
    except Exception: pass
    try:
        r = B.farkas_full(A, b, Eq, c, lo, hi, nv)
        return bool(r) and B.exact_basis(A, b, Eq, c, lo, hi, nv, r[0], r[1])
    except Exception: return False
