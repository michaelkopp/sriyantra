"""Step 6. Refine the tied-tip search at misses 0.14, 0.15, 0.16 (writes tip2b.json; the figure uses 0.15).
Usage: RINGOBJ=range python3 tip2b.py"""
import json,math
import os, sys
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
from nmtip2 import *
y0=[r for r in json.load(open(os.path.join(HERE,'nmtip2.json'))) if r['g0']==0.14][0]['y']
out=[]
for g0 in (0.14,0.15,0.16):
    u,s=nm2(uof(y0),0.06,g0,rounds=10)
    y=ylad(u); L,tips,W=lam_margins(y); sh=rings(y); rr={k:round(float(max(sh[k])/min(sh[k])),2) for k in RINGS}
    print('g0 %.2f score %.3f gap %.4f minrel %.4f'%(g0,s,gap(y),math.expm1(min(L.values()))),rr,min(L,key=L.get),flush=True)
    out.append(dict(g0=g0,y=y,score=s))
json.dump(out,open(os.path.join(HERE,'tip2b.json'),'w'))
