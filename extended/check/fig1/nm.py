"""Ring-uniformity search on the family (Nelder-Mead in a local chart of the 4-dimensional family, constraints as penalties).
Usage: RINGOBJ=range python3 nm.py 0.06 0.17,0.14,0.12,0.10 nmD.json"""
"""Nelder-Mead on the 4-dimensional family (local chart + projection)."""
import os, json,sys,math,numpy as np, warnings
warnings.filterwarnings('ignore')
from scipy.optimize import minimize
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
from model import *
from maximin import ylad, uof, EQ
from opt import lam_margins, objective
def eqv(u):
    r=resid(ylad(u)); return np.array([r[k] for k in EQ])
def jac(u,h=1e-7):
    r0=eqv(u); J=np.zeros((5,len(u)))
    for i in range(len(u)):
        e=np.zeros(len(u)); e[i]=h; J[:,i]=(eqv(u+e)-r0)/h
    return J
def project(u,tol=1e-13,it=30):
    u=np.array(u,float)
    for _ in range(it):
        r=eqv(u)
        if not np.all(np.isfinite(r)): return None
        if max(abs(r))<tol: return u
        J=jac(u); u=u-np.linalg.lstsq(J,r,rcond=None)[0]
    r=eqv(u)
    return u if np.all(np.isfinite(r)) and max(abs(r))<1e-10 else None
def score(y,m0,g0,gmin=0.03,tau=0.02,taumax=0.10,wm=0.0):
    L,tips,W=lam_margins(y); g=gap(y,W)
    pen=0.0
    mm=min(L.values())
    pen+=max(0,m0-mm)*100
    for t in tips: pen+=max(0,tau-t)*100+max(0,t-taumax)*100
    pen+=max(0,g-g0)*100+max(0,gmin-g)*100
    try: o=objective(y,W)
    except Exception: o=None
    if o is None: return 1e3+pen
    return o+pen-wm*mm
def nm(u0,m0,g0,rounds=6,scale=0.05,**kw):
    u=project(u0)
    for r in range(rounds):
        J=jac(u); _,_,Vt=np.linalg.svd(J); N=Vt[5:].T
        def f(s):
            v=project(u+N@s)
            if v is None: return 1e4
            return score(ylad(v),m0,g0,**kw)
        res=minimize(f,np.zeros(4),method='Nelder-Mead',options={'initial_simplex':np.vstack([np.zeros(4),scale*np.eye(4)]),'maxfev':400,'xatol':1e-5,'fatol':1e-6})
        v=project(u+N@res.x)
        if v is not None: u=v
        scale*=0.6
    return u,score(ylad(u),m0,g0,**kw)
if __name__=='__main__':
    m0=float(sys.argv[1]); g0s=[float(x) for x in sys.argv[2].split(',')]; out=sys.argv[3]
    F=json.load(open(os.path.join(HERE,'front.json'))); fr={r['g0']:r for r in F}
    u=uof(fr[0.17]['y']); res=[]
    for g0 in g0s:
        u,s=nm(u,m0,g0)
        y=ylad(u); L,tips,W=lam_margins(y)
        print('g0 %.2f score %.3f obj %s gap %.3f minm %.3f tips %.3f %.3f'%(g0,s,objective(y),gap(y),min(L.values()),tips[0],tips[1]),flush=True)
        res.append(dict(g0=g0,m0=m0,score=s,y=y,gap=gap(y),minm=min(L.values()),tips=tips))
        json.dump(res,open(os.path.join(HERE,out),'w'))
