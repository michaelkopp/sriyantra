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
                    # z >= x + cst + s (y - x)  ->  -z + (1-s) x + s y <= -cst
                    r={z:-1.0}; r[x]=r.get(x,0)+(1-s); r[y]=r.get(y,0)+s
                    A.append(r); b.append(-cst+1e-11*(1+abs(cst)))
                if math.isfinite(u0) and math.isfinite(u1) and u1-u0>1e-12:
                    s=(sp(u1)-sp(u0))/(u1-u0); cst=sp(u0)-s*u0
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
def farkas(A,b,E,c,lo,hi,nv):
    m,k=len(A),len(E); Ad,Ed=dense(A,nv),dense(E,nv)
    cost=np.concatenate([np.array(b),np.array(c),-np.array([l if math.isfinite(l) else 0.0 for l in lo]),np.array([h if math.isfinite(h) else 0.0 for h in hi])])
    Aeq=np.hstack([Ad.T,Ed.T if k else np.zeros((nv,0)),-np.eye(nv),np.eye(nv)])
    bounds=[(0,None)]*m+[(-1,1)]*k
    bounds+=[(0,0) if not math.isfinite(l) else (0,None) for l in lo]
    bounds+=[(0,0) if not math.isfinite(h) else (0,None) for h in hi]
    Aub=np.concatenate([np.ones(m),np.zeros(k),np.ones(2*nv)])[None,:]
    res=linprog(cost,A_ub=Aub,b_ub=[1.0],A_eq=Aeq,b_eq=np.zeros(nv),bounds=bounds,method='highs')
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
def refute_bounds(bounds,opt={}):
    A,b,E,c,lo,hi=build(bounds,opt)
    r=farkas(A,b,E,c,lo,hi,NV)
    if r is None: return False
    y,z=r
    for den in (None,1024,2**20,64,8,1):
        if verify(A,b,E,c,lo,hi,y,z,den): return True
    return False
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
    for den in (None,1024,2**20,64,8,1):
        if verify(A,b,E,c,lo,hi,y,z,den): return True
    return False
def handoff(bounds,opt={}):
    """True if every point of the box on the relaxed V7 lies in one of the regions certified in chart B
    (RB_DEFAULT, or the first or second of RB_BOXES): each piece of that region's complement is refuted"""
    for H in (HALFS if opt.get('regions',True) else [HALF]):
        if all(refute_rows(bounds,[(r,rhs)],opt) for r,rhs in H): return True
    return False
