"""Direct check of the gradient enclosures in mvf.py: at random points of random boxes (finite, wide and
half-infinite coordinates) in random cones, the derivatives of every correction x_v with respect to s_i (finite
coordinates; central differences) and u_i = exp(s_i) (the others; closed form) must lie in the enclosures."""
import json, math, sys, random
sys.path.insert(0,'.')
import numpy as np
import cone as C, mvf as M
R=json.load(open('roots_all_surv.json')); random.seed(int(sys.argv[1]) if len(sys.argv)>1 else 1)
INF=math.inf; bad=0; chk=0
def xval(W, s):
    v=[w@s for w in W]; mx=max([0.0]+v); return mx+math.log(math.exp(-mx)+sum(math.exp(x-mx) for x in v))
for trial in range(int(sys.argv[2]) if len(sys.argv)>2 else 40):
    T=C.Tree(R[random.randrange(len(R))]); box=[]
    for i in range(10):
        r=random.random(); b=-random.expovariate(0.7)
        box.append((-INF,b) if r<0.25 else ((b-random.uniform(1.0,6.0),b) if r<0.5 else (b-random.uniform(0.001,0.9),b)))
    Mv=M.MVF(T,box)
    for (a_,b_),(m,others) in T.pair.items():
        W=[np.asarray(w,dtype=float) for _,w in others]
        if not W: continue
        k=C.PAIRIDX[(a_,b_)]; x,gF,gU,gFl,gFh,gUl,gUh,_=Mv.px[k]
        for rep in range(3):
            s=np.array([ (b-random.expovariate(1.0)*3 if a==-INF else random.uniform(a,b)) for a,b in box])
            h=1e-6
            for t,i in enumerate(Mv.FN):
                e=np.zeros(10); e[i]=h; g=(xval(W,s+e)-xval(W,s-e))/(2*h); chk+=1
                if not (gFl[t]-1e-6<=g<=gFh[t]+1e-6): bad+=1
            S=1+sum(math.exp(w@s) for w in W)
            for t,i in enumerate(Mv.TL):     # d x / d u_i = sum over the paths through i of exp(w_j - s_i) / S
                g=sum(math.exp(w@s-s[i]) for w in W if w[i])/S; chk+=1
                if not (gUl[t]-1e-6*(1+abs(g))<=g<=gUh[t]+1e-6*(1+abs(g))): bad+=1
print('derivatives checked',chk,'outside the enclosure',bad); print('PASS' if bad==0 else 'FAIL')
