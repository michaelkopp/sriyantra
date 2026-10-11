"""LP certificates for chart A boxes (gauge y7 = inf), verified in exact rational arithmetic.

Variables: L_ab = log|p_a - p_b| for the 66 pairs of finite points.  For three points i < j < k on the
line, L_ik = LSE(L_ij, L_jk) = L_ij + sp(L_jk - L_ij), sp(u) = log(1 + e^u) convex; this is relaxed
by tangent lines of sp (lower bounds) and its chord over the enclosure of u (upper bound).  Every
condition and clause is linear in L (pair vectors of tri.json).  Bounds on L: the enclosures that
chart A computes for the box (after its solve steps).  Trinomials e^m1 + e^m2 = 1 and clauses
e^ma + e^mb <= 1 are relaxed with tangent lines of fixed dyadic slope."""
import json, math, numpy as np
from fractions import Fraction as Q
from scipy.optimize import linprog
import os
D=json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'tri.json')))
PAIRS=[tuple(p) for p in D['pairs']]
INF=float('inf'); INFPT=8
LPTIME=float(os.environ.get('LPTIME','20'))
ORDER=[9,10,11,12,0,1,2,3,4,5,6,7]
PTS=ORDER
VI={}
for i in range(len(PTS)):
    for j in range(i+1,len(PTS)): VI[(PTS[i],PTS[j])]=len(VI)
NV=len(VI)
def vidx(a,b):
    return VI[(a,b)] if (a,b) in VI else VI[(b,a)]
def row_of(vec):
    r={}
    for (a,b),e in zip(PAIRS,vec):
        if e==0 or a==INFPT or b==INFPT: continue
        k=vidx(a,b); r[k]=r.get(k,0)+float(e)
    return r
MONO={n:row_of(v) for n,v in D['mono'].items()}
CLS=[(n,c['type'],[row_of(m) for m in c['m']]) for n,c in D['clauses'].items()]
def sp(u): return u+math.log1p(math.exp(-u)) if u>0 else math.log1p(math.exp(u))
def sg(u): return 1/(1+math.exp(-u))
def f_(s):
    if s>=0: return -INF
    return math.log(-math.expm1(s)) if s>-0.693 else math.log1p(-math.exp(s))
SLOPES=[-1/64,-1/16,-1/4,-1/2,-1,-2,-4,-16,-64]
def tangent_const(sig):
    es=sig/(sig-1); s=math.log(es); c=f_(s)-sig*s
    return c+1e-12*(1+abs(c))
def build(bounds,opt={}):
    lo=[bounds[k][0] for k in range(NV)]; hi=[bounds[k][1] for k in range(NV)]
    A=[];b=[];E=[];c=[]
    P=PTS
    for i in range(len(P)):
        for j in range(i+1,len(P)):
            for k in range(j+1,len(P)):
                x=vidx(P[i],P[j]); y=vidx(P[j],P[k]); z=vidx(P[i],P[k])
                # z = x + sp(y - x);  z >= x, z >= y
                if opt.get('maxrows',True): A.append({x:1.0,z:-1.0}); b.append(0.0); A.append({y:1.0,z:-1.0}); b.append(0.0)
                u0=lo[y]-hi[x]; u1=hi[y]-lo[x]
                if math.isnan(u0): u0=-INF
                if math.isnan(u1): u1=INF
                pts=[]
                NT=opt.get('ntan',3)
                if math.isfinite(u0) and math.isfinite(u1): pts=list(np.linspace(u0,u1,NT)) if NT>1 else [0.5*(u0+u1)]
                elif math.isfinite(u1): pts=[u1,u1-2,u1-6]
                elif math.isfinite(u0): pts=[u0,u0+2,u0+6]
                else: pts=[-2,0,2]
                for ut in pts:
                    s=sg(ut); cst=sp(ut)-s*ut
                    # slopes within 1e-10 of 0 or 1: the exact rows z >= x, z >= y instead (no tiny coefficients)
                    if s<1e-10: A.append({x:1.0,z:-1.0}); b.append(0.0); continue
                    if s>1-1e-10: A.append({y:1.0,z:-1.0}); b.append(0.0); continue
                    # z >= x + cst + s (y - x)  ->  -z + (1-s) x + s y <= -cst
                    r={z:-1.0}; r[x]=r.get(x,0)+(1-s); r[y]=r.get(y,0)+s
                    A.append(r); b.append(-cst+1e-11*(1+abs(cst)))
                if math.isfinite(u0) and math.isfinite(u1) and u1-u0>1e-12:
                    s=(sp(u1)-sp(u0))/(u1-u0); cst=sp(u0)-s*u0
                    if s<1e-10: A.append({z:1.0,x:-1.0}); b.append(sp(u1)+1e-11*(1+abs(sp(u1))))
                    elif s>1-1e-10: A.append({z:1.0,y:-1.0}); b.append(sp(-u0)+1e-11*(1+abs(sp(-u0))))
                    else:
                        r={z:1.0}; r[x]=r.get(x,0)-(1-s); r[y]=r.get(y,0)-s
                        A.append(r); b.append(cst+1e-11*(1+abs(cst)+abs(s*u0)))
                elif opt.get('halfinf',True):
                    # one-sided enclosure of u: sp(u) <= sp(u1) for u <= u1, sp(u) <= u + sp(-u0) for u >= u0
                    if math.isfinite(u1): A.append({z:1.0,x:-1.0}); b.append(sp(u1)+1e-11*(1+abs(sp(u1))))
                    if math.isfinite(u0): A.append({z:1.0,y:-1.0}); b.append(sp(-u0)+1e-11*(1+abs(sp(-u0))))
    for n in ['Q21','Q43','Q59','apexU','apexD']+([] if opt.get('notarget') else ['target']):
        E.append(MONO[n]); c.append(0.0)
    def comb(r1,r2,s2):
        r=dict(r1)
        for k,x in r2.items(): r[k]=r.get(k,0)+s2*x
        return r
    for n in ['DOWN','UP']:
        r1,r2=MONO[n+'.m1'],MONO[n+'.m2']
        A.append(r1); b.append(0.0); A.append(r2); b.append(0.0)
        for (x,y) in ((r1,r2),(r2,r1)):
            for s in SLOPES: A.append(comb(x,y,-s)); b.append(tangent_const(s))
    if not opt.get('noclauses'):
        for n,t,ms in CLS:
            if t=='bin': A.append({k:-x for k,x in ms[0].items()}); b.append(0.0)
            else:
                ra,rb=ms; A.append(ra); b.append(0.0); A.append(rb); b.append(0.0)
                for (x,y) in ((ra,rb),(rb,ra)):
                    for s in SLOPES: A.append(comb(x,y,-s)); b.append(tangent_const(s))
    return A,b,E,c,lo,hi
def dense(rows,nv):
    M=np.zeros((len(rows),nv))
    for i,r in enumerate(rows):
        for k,x in r.items(): M[i,k]=x
    return M
def farkas_full(A,b,E,c,lo,hi,nv):
    """the dual (Farkas) LP: returns its solution x = (y, z, p, q) and the data, or None"""
    m,k=len(A),len(E); Ad,Ed=dense(A,nv),dense(E,nv)
    cost=np.concatenate([np.array(b),np.array(c),-np.array([l if math.isfinite(l) else 0.0 for l in lo]),np.array([h if math.isfinite(h) else 0.0 for h in hi])])
    Aeq=np.hstack([Ad.T,Ed.T if k else np.zeros((nv,0)),-np.eye(nv),np.eye(nv)])
    bounds=[(0,None)]*m+[(-1,1)]*k
    bounds+=[(0,0) if not math.isfinite(l) else (0,None) for l in lo]
    bounds+=[(0,0) if not math.isfinite(h) else (0,None) for h in hi]
    Aub=np.concatenate([np.ones(m),np.zeros(k),np.ones(2*nv)])[None,:]
    res=linprog(cost,A_ub=Aub,b_ub=[1.0],A_eq=Aeq,b_eq=np.zeros(nv),bounds=bounds,method='highs',options={'time_limit':LPTIME})
    if res.status!=0 or res.fun>-1e-9: return None
    return res.x,bounds
def exact_basis(A,b,E,c,lo,hi,nv,x,bounds,tol=1e-11):
    """Exact rational certificate from the support of the floating dual solution x: variables at a bound are
    fixed there exactly, the others solve the equations A^T y + E^T z - p + q = 0 exactly (a basis of them chosen by
    pivoted QR, the rest fixed at rounded values); then the certificate is checked exactly."""
    m,k=len(A),len(E); n=m+k+2*nv
    cols=[]                                   # column of each dual variable as {row j: coef}
    for r in A: cols.append(dict(r))
    for r in E: cols.append(dict(r))
    for j in range(nv): cols.append({j:-1.0})
    for j in range(nv): cols.append({j:1.0})
    val=[None]*n; free=[]
    for i in range(n):
        lb,ub=bounds[i]; xi=float(x[i])
        if lb is not None and abs(xi-lb)<=tol: val[i]=Q(lb)
        elif ub is not None and abs(xi-ub)<=tol: val[i]=Q(ub)
        else: free.append(i)
    if not free: return False
    Mf=np.zeros((nv,len(free)))
    for t,i in enumerate(free):
        for j,v in cols[i].items(): Mf[j,t]=v
    from scipy.linalg import qr
    _,R_,piv=qr(Mf,pivoting=True,mode='economic')
    d=np.abs(np.diag(R_)) if R_.size else np.array([])
    rank=int(np.sum(d>1e-9*max(1.0,d[0] if d.size else 1.0)))
    if rank==0: return False                    # (LAPACK rejects an empty pivoted QR below)
    basic=[free[piv[t]] for t in range(rank)]; nonb=[i for i in free if i not in basic]
    for i in nonb: val[i]=Q(float(x[i])).limit_denominator(2**40)
    # rows: choose rank independent equations
    Mb=np.zeros((nv,rank))
    for t,i in enumerate(basic):
        for j,v in cols[i].items(): Mb[j,t]=v
    _,R2,piv2=qr(Mb.T,pivoting=True,mode='economic')
    eqs=[int(piv2[t]) for t in range(rank)]
    rhs=[Q(0)]*nv
    for i in range(n):
        if i in basic or val[i]==0: continue
        for j,v in cols[i].items(): rhs[j]-=val[i]*Q(v)
    Aug=[[Q(0)]*rank+[rhs[j]] for j in eqs]
    for t,i in enumerate(basic):
        for j,v in cols[i].items():
            if j in eqs: Aug[eqs.index(j)][t]=Q(v)
    nn=rank
    for col in range(nn):
        pv=next((r_ for r_ in range(col,nn) if Aug[r_][col]!=0),None)
        if pv is None: return False
        Aug[col],Aug[pv]=Aug[pv],Aug[col]
        for r_ in range(nn):
            if r_!=col and Aug[r_][col]!=0:
                f=Aug[r_][col]/Aug[col][col]; Aug[r_]=[a_-f*b_ for a_,b_ in zip(Aug[r_],Aug[col])]
    for t,i in enumerate(basic): val[i]=Aug[t][nn]/Aug[t][t]
    # all equations must hold exactly, and the sign constraints
    for i in range(n):
        lb,ub=bounds[i]
        if lb is not None and val[i]<lb: return False
        if ub is not None and val[i]>ub: return False
    w=[Q(0)]*nv
    for i in range(n):
        if val[i]:
            for j,v in cols[i].items(): w[j]+=val[i]*Q(v)
    if any(v!=0 for v in w): return False
    # exact objective: b.y + c.z - lo.p + hi.q < 0  (p, q are zero where the bound is infinite)
    obj=Q(0)
    for i in range(m): obj+=val[i]*Q(b[i])
    for i in range(k): obj+=val[m+i]*Q(c[i])
    for j in range(nv):
        if val[m+k+j]: obj-=val[m+k+j]*Q(lo[j])
        if val[m+k+nv+j]: obj+=val[m+k+nv+j]*Q(hi[j])
    return obj<0
def farkas(A,b,E,c,lo,hi,nv):
    m,k=len(A),len(E); Ad,Ed=dense(A,nv),dense(E,nv)
    cost=np.concatenate([np.array(b),np.array(c),-np.array([l if math.isfinite(l) else 0.0 for l in lo]),np.array([h if math.isfinite(h) else 0.0 for h in hi])])
    Aeq=np.hstack([Ad.T,Ed.T if k else np.zeros((nv,0)),-np.eye(nv),np.eye(nv)])
    bounds=[(0,None)]*m+[(-1,1)]*k
    bounds+=[(0,0) if not math.isfinite(l) else (0,None) for l in lo]
    bounds+=[(0,0) if not math.isfinite(h) else (0,None) for h in hi]
    Aub=np.concatenate([np.ones(m),np.zeros(k),np.ones(2*nv)])[None,:]
    res=linprog(cost,A_ub=Aub,b_ub=[1.0],A_eq=Aeq,b_eq=np.zeros(nv),bounds=bounds,method='highs',options={'time_limit':LPTIME})
    if res.status!=0 or res.fun>-1e-9: return None
    x=res.x; return x[:m],x[m:m+k]
def verify(A,b,E,c,lo,hi,y,z,den):
    rat=(lambda v: Q(float(v)).limit_denominator(den)) if den else (lambda v: Q(float(v)))
    yr=[max(Q(0),rat(v)) for v in y]; zr=[rat(v) for v in z]
    w={}; rhs=Q(0)
    for yi,r,bi in zip(yr,A,b):
        if yi:
            rhs+=yi*Q(bi)
            for k,x in r.items(): w[k]=w.get(k,Q(0))+yi*Q(x)
    for zi,r,ci in zip(zr,E,c):
        if zi:
            rhs+=zi*Q(ci)
            for k,x in r.items(): w[k]=w.get(k,Q(0))+zi*Q(x)
    lb=Q(0)
    for k,wk in w.items():
        if wk>0:
            if not math.isfinite(lo[k]): return False
            lb+=wk*Q(lo[k])
        elif wk<0:
            if not math.isfinite(hi[k]): return False
            lb+=wk*Q(hi[k])
    return lb>rhs
def exact_check(A,b,E,c,lo,hi,yr,zr):
    """exact Farkas check with rational multipliers yr >= 0, zr"""
    if any(v<0 for v in yr): return False
    w={}; rhs=Q(0)
    for yi,r,bi in zip(yr,A,b):
        if yi:
            rhs+=yi*Q(bi)
            for k,x in r.items(): w[k]=w.get(k,Q(0))+yi*Q(x)
    for zi,r,ci in zip(zr,E,c):
        if zi:
            rhs+=zi*Q(ci)
            for k,x in r.items(): w[k]=w.get(k,Q(0))+zi*Q(x)
    lb=Q(0)
    for k,wk in w.items():
        if wk>0:
            if not math.isfinite(lo[k]): return False
            lb+=wk*Q(lo[k])
        elif wk<0:
            if not math.isfinite(hi[k]): return False
            lb+=wk*Q(hi[k])
    return lb>rhs
def repair(A,b,E,c,lo,hi,y,z,nv,den=2**30,tol=1e-9):
    """Round the multipliers to rationals and correct a few of them exactly so that the columns of
    variables without the needed finite bound cancel exactly; then check exactly."""
    rows=[(r,True) for r in A]+[(r,False) for r in E]
    m=[float(v) for v in y]+[float(v) for v in z]
    w=np.zeros(nv)
    for mi,(r,_) in zip(m,rows):
        if mi:
            for k,x in r.items(): w[k]+=mi*x
    K=[]
    for j in range(nv):
        fl,fh=math.isfinite(lo[j]),math.isfinite(hi[j])
        if fl and fh: continue
        if abs(w[j])<tol: K.append(j)
        elif (w[j]>0 and not fl) or (w[j]<0 and not fh): return False
    mr=[Q(v).limit_denominator(den) if abs(v)>1e-15 else Q(0) for v in m]
    for i in range(len(A)):
        if mr[i]<0: mr[i]=Q(0)
    if K:
        S=[i for i,v in enumerate(m) if (v>1e-7 if rows[i][1] else abs(v)>1e-7)]
        M=np.zeros((len(K),len(S))); kp={j:a for a,j in enumerate(K)}
        for c_,i in enumerate(S):
            for k,x in rows[i][0].items():
                if k in kp: M[kp[k],c_]=x
        if len(S)<len(K): return False
        from scipy.linalg import qr
        _,R_,piv=qr(M,pivoting=True,mode='economic')
        rank=int(np.sum(np.abs(np.diag(R_))>1e-9*max(1.0,abs(R_[0,0])))) if R_.size else 0
        if rank==0: return False                # (LAPACK rejects an empty pivoted QR below)
        cols=[S[piv[t]] for t in range(rank)]
        # rows of K that are independent: solve the least-norm exact system on the chosen columns
        Mi=[[Q(0)]*len(cols) for _ in K]
        for c_,i in enumerate(cols):
            for k,x in rows[i][0].items():
                if k in kp: Mi[kp[k]][c_]=Q(x)
        wk=[Q(0)]*len(K)
        for mi,(r,_) in zip(mr,rows):
            if mi:
                for k,x in r.items():
                    if k in kp: wk[kp[k]]+=mi*Q(x)
        # choose rank independent rows of Mi
        Mf=np.array([[float(v) for v in row] for row in Mi])
        _,R2,piv2=qr(Mf.T,pivoting=True,mode='economic')
        rws=[piv2[t] for t in range(rank)]
        sq=[Mi[r_] for r_ in rws]; rhs=[-wk[r_] for r_ in rws]
        n=rank; Aug=[list(sq[i])+[rhs[i]] for i in range(n)]
        for col in range(n):
            pv=next((r_ for r_ in range(col,n) if Aug[r_][col]!=0),None)
            if pv is None: return False
            Aug[col],Aug[pv]=Aug[pv],Aug[col]
            for r_ in range(n):
                if r_!=col and Aug[r_][col]!=0:
                    f=Aug[r_][col]/Aug[col][col]; Aug[r_]=[a_-f*b_ for a_,b_ in zip(Aug[r_],Aug[col])]
        d=[Aug[i][n]/Aug[i][i] for i in range(n)]
        for t,i in enumerate(cols): mr[i]+=d[t]
        # the remaining rows of K must now cancel as well
        wk2=[Q(0)]*len(K)
        for mi,(r,_) in zip(mr,rows):
            if mi:
                for k,x in r.items():
                    if k in kp: wk2[kp[k]]+=mi*Q(x)
        if any(v!=0 for v in wk2): return False
    return exact_check(A,b,E,c,lo,hi,mr[:len(A)],mr[len(A):])
def certify(A,b,E,c,lo,hi,y,z,nv=None):
    for den in (None,1024,2**20,64,8,1):
        if verify(A,b,E,c,lo,hi,y,z,den): return True
    try:
        if repair(A,b,E,c,lo,hi,y,z,nv or NV): return True
    except Exception: pass
    try:
        r=farkas_full(A,b,E,c,lo,hi,nv or NV)
        return bool(r) and exact_basis(A,b,E,c,lo,hi,nv or NV,r[0],r[1])
    except Exception: return False
def refute_bounds(bounds,opt={}):
    A,b,E,c,lo,hi=build(bounds,opt)
    r=farkas(A,b,E,c,lo,hi,NV)
    if r is None: return False
    y,z=r
    return certify(A,b,E,c,lo,hi,y,z)
def bounds_from_dump(rec):
    bd=[None]*NV
    for key,(l,h) in rec['br'].items():
        a,bb=map(int,key.split(','))
        bd[vidx(a,bb)]=(-INF if l is None else l, INF if h is None else h)
    return bd
if __name__=='__main__':
    import sys, time
    from collections import Counter
    d=json.load(open(sys.argv[1])); n=int(sys.argv[2]) if len(sys.argv)>2 else len(d['boxes'])
    cnt=Counter(); t0=time.time(); why=Counter()
    for rec in d['boxes'][:n]:
        ok=refute_bounds(bounds_from_dump(rec),{'noclauses':'--nocl' in sys.argv})
        cnt[ok]+=1; why[(rec['why'] or 'none',ok)]+=1
    print(dict(cnt),'%.3f s per box'%((time.time()-t0)/n)); print(why.most_common(12))
# ---------- handoff to chart B: the certified region R_B in fan coordinates ----------
# lambda_k = log t_k = L(k,0) + L(k+1,12) - L(k+1,0) - L(k,12)  (13-point labels; brackets with y7 drop out)
def lamrow(k):
    r={}
    for (a,bb),s in (((k,0),1),((k+1,12),1),((k+1,0),-1),((k,12),-1)):
        if a==INFPT or bb==INFPT: continue
        i=vidx(a,bb); r[i]=r.get(i,0)+s
    return {i:x for i,x in r.items() if x}
LAMROW={k:lamrow(k) for k in range(1,11)}
def lincomb(spec):
    r={}
    for k,s in spec.items():
        for i,x in LAMROW[k].items(): r[i]=r.get(i,0)+s*x
    return {i:x for i,x in r.items() if x}
RB_FORMS={'e':{3:1},'d':{6:1,3:-1},'r2':{2:1,3:-1},'r4':{4:1,6:-1},'r5':{5:1,6:1,3:-1},'r7':{7:1,6:1,3:-1},'r8':{8:1,3:-1},'r9':{9:1,3:-1}}
def complement_halfspaces(RB):
    """RB = {'e':(lo,hi),...} -> list of (row, rhs) with row.L <= rhs describing each piece of the complement"""
    out=[]
    for n,(lo,hi) in RB.items():
        row=lincomb(RB_FORMS[n])
        if hi is not None: out.append(({i:-x for i,x in row.items()},-hi))     # form >= hi
        if lo is not None: out.append((row,lo))                                # form <= lo
    return out
# certified in chartB (fanRBB.js): {t3<=0.1, d>=2, |r|<=2}, {t3<=0.1, 0.5<=d<=2.2, |r|<=0.4}, {0.1<=t3<=0.2, d>=1, |r|<=0.25}
RB_DEFAULT={'e':(None,math.log(0.2)),'d':(1.0,None),'r2':(-0.25,0.25),'r4':(-0.25,0.25),'r5':(-0.25,0.25),'r7':(-0.25,0.25),'r8':(-0.25,0.25),'r9':(-0.25,0.25)}
HALF=complement_halfspaces(RB_DEFAULT)
# the three boxes certified in chart B (their union contains RB_DEFAULT)
RB_BOXES=[{'e':(None,-2.3026),'d':(2.0,None),**{r:(-2.0,2.0) for r in ('r2','r4','r5','r7','r8','r9')}},
          {'e':(None,-2.3026),'d':(0.5,2.2),**{r:(-0.4,0.4) for r in ('r2','r4','r5','r7','r8','r9')}},
          {'e':(-2.3026,-1.6095),'d':(1.0,None),**{r:(-0.25,0.25) for r in ('r2','r4','r5','r7','r8','r9')}}]
HALFS=[HALF]+[complement_halfspaces(R) for R in RB_BOXES[:2]]
def refute_rows(bounds,extra,opt={}):
    A,b,E,c,lo,hi=build(bounds,opt)
    for r,rhs in extra: A.append(r); b.append(rhs)
    res=farkas(A,b,E,c,lo,hi,NV)
    if res is None: return False
    y,z=res
    return certify(A,b,E,c,lo,hi,y,z)
DK=(2,5,7,8,9)
def in_region_D(bounds,c,theta,opt={}):
    """True if every point of the box on the relaxed V7 lies in D = {t3 <= c, t_k <= theta (k in DK)},
    the region certified in chart B (fanD.js): each half-space t3 >= c, t_k >= theta is refuted"""
    rows=[({i:-x for i,x in LAMROW[3].items()},-math.log(c))]+[({i:-x for i,x in LAMROW[k].items()},-math.log(theta)) for k in DK]
    return all(refute_rows(bounds,[r],opt) for r in rows)
def handoff(bounds,opt={}):
    """True if every point of the box on the relaxed V7 lies in one of the regions certified in chart B
    (RB_DEFAULT, or the first or second of RB_BOXES): each piece of that region's complement is refuted"""
    for H in (HALFS if opt.get('regions',True) else [HALF]):
        if all(refute_rows(bounds,[(r,rhs)],opt) for r,rhs in H): return True
    return False
