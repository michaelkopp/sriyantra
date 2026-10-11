import sys, time, json, math
import os; os.environ['CELLS']='all'
sys.path.insert(0,'.')
import cone as C
from collections import Counter
t0=time.time(); cnt=Counter(); surv=[]
for T in C.all_trees():
    E=C.Enc(T,[(-math.inf,0.0)]*10)
    cs=C.cells_hit(E)
    if not cs: cnt['cell']+=1; continue
    w=C.interval_refute(E)
    if w: cnt['iv:'+w]+=1; continue
    pass
    pass
    cnt['survive']+=1; surv.append(T.parent)
print(dict(cnt), '%.0fs'%(time.time()-t0))
json.dump(surv,open('roots_all_surv.json','w'))
