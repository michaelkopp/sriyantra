"""Maximise the smallest clearance over the family, from the laminar points of the certifier. Writes maximin.json."""
import os, json,sys,math,numpy as np, warnings
warnings.filterwarnings('ignore')
from scipy.optimize import minimize
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
from model import *
EQ=['E2','E3','E4','CD','CU']
def ylad(u):
    g=np.exp(np.clip(np.concatenate([[0.0],u]),-40,40)); c=np.concatenate([[0.0],np.cumsum(g)]); return list(-1+2*c/c[-1])
def uof(y):
    g=np.diff(y); return np.log(g[1:]/g[0])
def marg(y):
    W=Wtab(y); M=margins(y,W)
    zD,zU=ext(y,W)
    M['tip U']=math.log1p(max(zU-y[10],1e-300)/2) if zU>y[10] else -1
    M['tip D']=math.log1p(max(y[0]-zD,1e-300)/2) if zD<y[0] else -1
    return M
def normalise(t,lam):
    """fan-chart point t (zD=0, zU=inf, y10=1) -> levels (zD, y0..y10, zU) after x->x/(1+lam x), y0=-1, y10=1"""
    p=[0.0]*12; p[11]=1.0
    for k in range(10,0,-1): p[k]=t[k-1]*p[k+1]
    mu=lambda x: x/(1+lam*x)
    q=[mu(x) for x in p]+[1/lam]
    A=2/(q[11]-q[1]); B=-1-A*q[1]
    return [A*x+B for x in q]
KEYS=None
def run(y0,tipw=1.0):
    global KEYS
    u0=uof(y0); KEYS=sorted(marg(y0))
    def mv(v):
        y=ylad(v[:9]); M=marg(y); return np.array([M[k]*(tipw if k.startswith('tip') else 1) for k in KEYS])
    def eqf(v):
        y=ylad(v[:9]); r=resid(y); return np.array([r[k] for k in EQ])
    t0=min(mv(np.concatenate([u0,[0]])))
    x0=np.concatenate([u0,[t0]])
    cons=[{'type':'eq','fun':eqf},{'type':'ineq','fun':lambda v: mv(v)-v[9]}]
    res=minimize(lambda v:-v[9],x0,method='SLSQP',constraints=cons,options={'maxiter':500,'ftol':1e-12})
    y=ylad(res.x[:9]); return y,res
if __name__=='__main__':
    # starting points: the laminar points of the certifier (../../certifier/cones/pts_t.json)
    T=json.load(open(os.path.join(HERE,'..','..','certifier','cones','pts_t.json')))
    import random; random.seed(1)
    out=[]
    idx=random.sample(range(len(T)),40)
    for i in idx:
        for lam in (0.3,3.0):
            y=normalise(T[i],lam)[1:12]
            try: y,res=run(y)
            except Exception as e: print('fail',i,e); continue
            M=marg(y); r=resid(y)
            print(i,lam,'t=%.4f'%res.x[9],'eqres %.1e'%max(abs(r[k]) for k in EQ),'gap %.3f'%gap(y),'argmin',min(M,key=M.get),flush=True)
            out.append(dict(i=i,lam=lam,t=float(res.x[9]),y=y,gap=gap(y)))
    json.dump(out,open(os.path.join(HERE,'maximin.json'),'w'))
