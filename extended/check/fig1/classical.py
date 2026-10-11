"""Step 8. A classical straight yantra (all twelve concurrences, laminar): the Figure-1 point projected onto E1..E4 (min-norm Newton); writes
../classical_pt.json, which ../drawkites.py draws."""
import json,math,numpy as np
import os, sys
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
from model import *
from maximin import ylad, uof
import nm
EQC=['E1','E2','E3','E4']
def eqv(u):
    r=resid(ylad(u)); return np.array([r[k] for k in EQC])
def jac(u,h=1e-7):
    r0=eqv(u); J=np.zeros((4,len(u)))
    for i in range(len(u)):
        e=np.zeros(len(u)); e[i]=h; J[:,i]=(eqv(u+e)-r0)/h
    return J
nm.eqv=eqv; nm.jac=jac
d=json.load(open(os.path.join(HERE,'..','straight_pt.json'))); y0=d['q'][1:12]
u=nm.project(uof(y0),tol=1e-15,it=60); y=ylad(u); W=Wtab(y); zD,zU=ext(y,W)
R=resid(y,W); M=margins(y,W)
nc={k:math.expm1(v) for k,v in M.items() if k.startswith('nc')}
print('resid',{k:'%.1e'%v for k,v in R.items()})
print('min nc clearance %.3f (%s)'%(min(nc.values()),min(nc,key=nc.get)))
sh=rings(y,W,with11=False)
print({k:len(v) for k,v in sorted(sh.items())}, {k:round(max(sh[k])/min(sh[k]),2) for k in (1,3,5,7)})
# collinearity misses of the middle corners
lU=lambda t: W['D1']+(W['D3']-W['D1'])*(t-y[7])/(y[9]-y[7])
lD=lambda t: W['U3']+(W['U1']-W['U3'])*(t-y[1])/(y[3]-y[1])
print('D2 corner vs line D1-D3: %.4f   U2 corner vs line U1-U3: %.4f (relative)'%(W['D2']/lU(y[8])-1,W['U2']/lD(y[2])-1))
json.dump({'y':y,'W':W,'zD':zD,'zU':zU},open(os.path.join(HERE,'..','classical_pt.json'),'w'),indent=1)
