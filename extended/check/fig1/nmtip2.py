"""Step 5. Ring search on the family with the tip clearances tied to the end rungs:
   y0 - zD = y2 - y1 and zU - y10 = y10 - y9 (two more equations; the family is then two-dimensional).
Usage: RINGOBJ=range python3 nmtip2.py   (writes nmtip2.json)"""
import os, json, sys, math, numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
import nm
from model import *
from maximin import uof, ylad, EQ
from opt import lam_margins, objective
from scipy.optimize import minimize
EQ2=EQ+['TD','TU']
def eqv(u):
    y=ylad(u); r=resid(y); W=Wtab(y)
    try:
        zD,zU=ext(y,W)
        r['TD']=math.log((y[0]-zD)/(y[2]-y[1])) if y[0]>zD else float('nan')
        r['TU']=math.log((zU-y[10])/(y[10]-y[9])) if zU>y[10] else float('nan')
    except Exception: return np.full(len(EQ2),np.nan)
    return np.array([r[k] for k in EQ2])
def jac(u,h=1e-7):
    r0=eqv(u); J=np.zeros((len(EQ2),len(u)))
    for i in range(len(u)):
        e=np.zeros(len(u)); e[i]=h; J[:,i]=(eqv(u+e)-r0)/h
    return J
nm.eqv=eqv; nm.jac=jac
def nm2(u0,m0,g0,rounds=8,scale=0.05):
    u=nm.project(u0)
    if u is None: raise RuntimeError('projection failed')
    for r in range(rounds):
        J=jac(u); _,_,Vt=np.linalg.svd(J); N=Vt[len(EQ2):].T
        def f(s):
            v=nm.project(u+N@s)
            if v is None: return 1e4
            return nm.score(ylad(v),m0,g0,tau=0.0,taumax=10.0)
        res=minimize(f,np.zeros(N.shape[1]),method='Nelder-Mead',options={'initial_simplex':np.vstack([np.zeros(N.shape[1]),scale*np.eye(N.shape[1])]),'maxfev':300,'xatol':1e-6,'fatol':1e-7})
        v=nm.project(u+N@res.x)
        if v is not None: u=v
        scale*=0.6
    return u,nm.score(ylad(u),m0,g0,tau=0.0,taumax=10.0)
if __name__=='__main__':
    starts={'nmtip':json.load(open(os.path.join(HERE,'nmtip.json')))[0]['y']}
    for r in json.load(open(os.path.join(HERE,'nmD.json'))): starts['nmD%.2f'%r['g0']]=r['y']
    out=[]
    for name,y0 in starts.items():
        for g0 in (0.14,0.17,0.25):
            try: u,s=nm2(uof(y0),0.06,g0)
            except Exception as e: print(name,g0,'fail',e); continue
            y=ylad(u); L,tips,W=lam_margins(y); zD,zU=ext(y,W)
            sh=rings(y); rr={k:round(float(max(sh[k])/min(sh[k])),1) for k in RINGS}
            print('%-8s g0 %.2f score %.3f gap %.3f minm %.3f  y0-zD %.3f y2-y1 %.3f  zU-y10 %.3f y10-y9 %.3f'%(name,g0,s,gap(y),min(L.values()),y[0]-zD,y[2]-y[1],zU-y[10],y[10]-y[9]),rr,flush=True)
            out.append(dict(start=name,g0=g0,score=s,y=y,gap=gap(y),minm=min(L.values()),rr=rr))
            json.dump(out,open(os.path.join(HERE,'nmtip2.json'),'w'))
