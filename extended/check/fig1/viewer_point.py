"""Step 8b. The yantra of Figure 1: Chiodo's straight construction at the viewer parameters of preset Fig1paper2
(lower corners on one line), converted to the ladder model."""
import json, sys, math, numpy as np
import os
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
from model import *
T=json.load(open(os.path.join(HERE,'T9.json')))
lev={}
for k,v in T.items(): lev.setdefault(round(v['apex'][1],12),[]).append(k+' apex'); lev.setdefault(round(v['right'][1],12),[]).append(k+' base')
ys=sorted(lev); print(len(ys)); 
for v in ys: print('%.6f'%v, lev[v])
y=ys
Wv={k:v['right'][0] for k,v in T.items()}
for k in Wv: assert abs(T[k]['apex'][1]-y[apex[k]])<1e-9 and abs(T[k]['right'][1]-y[base[k]])<1e-9, k
Wt=Wtab(y); kap=Wv['D1']/Wt['D1']
print('kappa %.6f'%kap, 'max width mismatch', max(abs(Wv[k]/kap-Wt[k]) for k in Wt))
R=resid(y,Wt); print({k:'%.1e'%v for k,v in R.items()})
M=margins(y,Wt); nc={k:math.expm1(v) for k,v in M.items() if k.startswith('nc')}
print('min nc %.3f %s'%(min(nc.values()),min(nc,key=nc.get)))
sh=rings(y,Wt,with11=False); print({k:len(v) for k,v in sorted(sh.items())},{k:round(max(sh[k])/min(sh[k]),2) for k in (1,3,5,7)})
lU=lambda t: Wt['D1']+(Wt['D3']-Wt['D1'])*(t-y[7])/(y[9]-y[7]); lD=lambda t: Wt['U3']+(Wt['U1']-Wt['U3'])*(t-y[1])/(y[3]-y[1])
print('D2 corner vs line D1-D3 %.4f   U2 corner vs line U1-U3 %.2e'%(Wt['D2']/lU(y[8])-1,Wt['U2']/lD(y[2])-1))
zD,zU=ext(y,Wt); print('zD %.4f zU %.4f'%(zD,zU))
json.dump({'source':'viewer preset Fig1paper2: Chiodo construction at exactp, lower corners collinear','exactp':json.load(open(os.path.join(HERE,'fig1paper2_preset.json')))['presets']['Fig1paper2']['exactp'],
           'y':y,'W':Wt,'kappa':kap,'zD':zD,'zU':zU},open(os.path.join(HERE,'..','classical_pt.json'),'w'),indent=1)
