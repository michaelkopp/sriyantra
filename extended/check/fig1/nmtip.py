"""Final step: re-run the ring search at miss <= 0.14 with tip clearances in [0.05,0.10] (and two variants). Usage: RINGOBJ=range python3 nmtip.py"""
import os, json,sys,os
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
from nm import *
from maximin import uof, ylad
r=[r for r in json.load(open(os.path.join(HERE,'nmD.json'))) if r['g0']==0.14][0]
u=uof(r['y']); res=[]
for tau,tmax in ((0.05,0.10),(0.07,0.10),(0.05,0.07)):
    u2,s=nm(u,0.06,0.14,tau=tau,taumax=tmax)
    y=ylad(u2); L,tips,W=lam_margins(y)
    sh=rings(y); rr={k:round(float(max(sh[k])/min(sh[k])),1) for k in RINGS}
    print('tau %.2f-%.2f score %.3f gap %.3f minm %.3f tips %.3f %.3f'%(tau,tmax,s,gap(y),min(L.values()),tips[0],tips[1]),rr,flush=True)
    res.append(dict(tau=tau,taumax=tmax,score=s,y=y))
json.dump(res,open(os.path.join(HERE,'nmtip.json'),'w'))
