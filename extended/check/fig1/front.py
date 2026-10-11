"""Pareto front: max min laminarity margin subject to gap <= g0 and tau <= tip clearance <= taumax."""
import os, json,sys,math,numpy as np, warnings
warnings.filterwarnings('ignore')
from scipy.optimize import minimize
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
from model import *
from maximin import ylad, uof, EQ
from opt import lam_margins, objective
def maxmargin(y0,g0,tau=0.02,taumax=0.12):
    u0=uof(y0); L0,_,_=lam_margins(y0); keys=sorted(L0)
    def parts(v):
        y=ylad(v[:9]); L,tips,W=lam_margins(y); return np.array([L[k] for k in keys]),tips,gap(y,W)
    def ineq(v):
        m,tips,g=parts(v)
        return np.concatenate([m-v[9],[tips[0]-tau,tips[1]-tau,taumax-tips[0],taumax-tips[1],g0-g]])
    def eqf(v):
        r=resid(ylad(v[:9])); return np.array([r[k] for k in EQ])
    m,_,_=parts(np.concatenate([u0,[0]]))
    x0=np.concatenate([u0,[min(m)]])
    res=minimize(lambda v:-v[9],x0,method='SLSQP',constraints=[{'type':'eq','fun':eqf},{'type':'ineq','fun':ineq}],
                 options={'maxiter':400,'ftol':1e-12})
    viol=min(min(ineq(res.x)),0); eq=max(abs(eqf(res.x)))
    return ylad(res.x[:9]),float(res.x[9]),float(viol),float(eq)
if __name__=='__main__':
    MM=json.load(open(os.path.join(HERE,'maximin.json'))); y=max(MM,key=lambda r:r['t'])['y']
    out=[]
    for tm in (0.5,0.3,0.2,0.12):
        y,t,viol,eq=maxmargin(y,1.0,taumax=tm); print('taumax',tm,'t %.4f viol %.1e eq %.1e gap %.3f'%(t,viol,eq,gap(y)),flush=True)
    for g0 in (0.3,0.25,0.2,0.17,0.15,0.13,0.11,0.1,0.09,0.08,0.07,0.06,0.05,0.04):
        y,t,viol,eq=maxmargin(y,g0)
        L,tips,W=lam_margins(y)
        o=objective(y)
        print('g0 %.2f t %.4f viol %.1e eq %.1e gap %.3f tips %.3f %.3f ringobj %s argmin %s'%(g0,t,viol,eq,gap(y),tips[0],tips[1],None if o is None else round(o,3),min(L,key=L.get)),flush=True)
        out.append(dict(g0=g0,t=t,viol=viol,eq=eq,y=y,gap=gap(y),obj=o))
        json.dump(out,open(os.path.join(HERE,'front.json'),'w'))
