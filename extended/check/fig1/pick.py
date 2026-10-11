"""Step 7. Polish the chosen point (miss 0.15 in tip2b.json) onto the family with the tied tips and write
../straight_pt.json, which ../drawstraight.py draws and ../straightfig.py checks."""
import os, json, sys
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
from nmtip2 import *
r=[r for r in json.load(open(os.path.join(HERE,'tip2b.json'))) if r['g0']==0.15][0]
u=nm.project(uof(r['y']),tol=1e-15,it=60); y=ylad(u); W=Wtab(y); zD,zU=ext(y,W)
json.dump({'source':'check/fig1 (ring-uniformity search on the family with y0-zD = y2-y1 and zU-y10 = y10-y9; miss <= 0.15, clearances >= 0.06)',
           'q':[zD]+list(y)+[zU],'W':W,'kappa':0.6},open(os.path.join(HERE,'..','straight_pt.json'),'w'),indent=1)
print('written')
