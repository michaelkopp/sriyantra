"""Clearances and the ring objective used by the figure search."""
import json,sys,math,os,numpy as np, warnings
warnings.filterwarnings('ignore')
from scipy.optimize import minimize
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
from model import *
from maximin import ylad, uof, EQ
LAMK=None
def lam_margins(y):
    W=Wtab(y); M=margins(y,W); zD,zU=ext(y,W)
    L={k:v for k,v in M.items() if not k.startswith('tip')}
    tips=((zU-y[10])/2, (y[0]-zD)/2)
    return L,tips,W
def objective(y,W=None):
    sh=rings(y,W)
    cnt={k:len(v) for k,v in sh.items()}
    if cnt!=EXPECT: return None
    s=0.0; mode=os.environ.get('RINGOBJ','std')
    for k in RINGS:
        a=np.log(np.array(sh[k]))
        s+= a.std() if mode=='std' else (a.max()-a.min())
    return s
def solve(y0,m0=0.05,tau=0.03,taumax=0.15,g0=0.12,gmin=0.03,maxiter=300):
    u0=uof(y0)
    L0,_,_=lam_margins(y0); keys=sorted(L0)
    def ineq(u):
        y=ylad(u); L,tips,W=lam_margins(y); g=gap(y,W)
        return np.array([L[k]-m0 for k in keys]+[tips[0]-tau,tips[1]-tau,taumax-tips[0],taumax-tips[1],g0-g,g-gmin])
    def eqf(u):
        r=resid(ylad(u)); return np.array([r[k] for k in EQ])
    def f(u):
        y=ylad(u)
        try: o=objective(y)
        except Exception: o=None
        return 50.0 if o is None else o
    res=minimize(f,u0,method='SLSQP',constraints=[{'type':'eq','fun':eqf},{'type':'ineq','fun':ineq}],
                 options={'maxiter':maxiter,'ftol':1e-9,'eps':1e-6})
    y=ylad(res.x)
    return y,res,float(min(ineq(res.x))),float(max(abs(eqf(res.x))))
