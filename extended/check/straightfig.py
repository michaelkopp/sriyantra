"""Check the point of Figure 2 (straight_pt.json): the equations that hold, the one that fails, every laminarity
inequality with its relative margin, the tip clearances, and the area spread in each ring of cells.
The point itself comes from the search in fig1/ (see fig1/README.md)."""
import os, json, sys, math
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,os.path.join(HERE,'fig1'))
from model import Wtab, resid, margins, ext, gap, rings, RINGS, EXPECT
d=json.load(open(os.path.join(HERE,'straight_pt.json'))); q=d['q']; y=q[1:12]; zD,zU=q[0],q[12]
W=Wtab(y)
assert all(a<b for a,b in zip(q,q[1:])), 'order'
z=ext(y,W); assert abs(z[0]-zD)<1e-12 and abs(z[1]-zU)<1e-12, 'apex relations'
R=resid(y,W)
print('log residuals: '+', '.join('%s %.1e'%(k,v) for k,v in R.items()))
M=margins(y,W)
cl=sorted(((k,math.expm1(v)) for k,v in M.items() if not k.startswith('tip')),key=lambda kv:kv[1])
print('25 laminarity inequalities, smallest relative margins:',', '.join('%s %.3f'%kv for kv in cl[:4]))
print('tips: z_U - y_10 = %.4f (y_10 - y_9 = %.4f), y_0 - z_D = %.4f (y_2 - y_1 = %.4f); yantra height 2'%(zU-y[10],y[10]-y[9],y[0]-zD,y[2]-y[1]))
print('miss of (t2,t6,t6): %.4f of the half-width of U4'%gap(y,W))
sh=rings(y,W); assert {k:len(v) for k,v in sh.items()}==EXPECT, 'cell structure'
names={2:'cells between yantra and fifth kite (14)',3:'ring of 14',5:'outer ring of 10',7:'inner ring of 10',9:'ring of 8'}
print('largest/smallest cell area: '+', '.join('%s %.2f'%(names[k],max(sh[k])/min(sh[k])) for k in RINGS))
