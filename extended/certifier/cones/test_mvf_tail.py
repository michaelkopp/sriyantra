"""Soundness of mvf.py on boxes with half-infinite coordinates: boxes around points of V8 (target, no clauses) and
of the laminar V7 (clauses, no target) must never be refuted."""
import json, math, sys, random
sys.path.insert(0,'.')
import cone as C, mvf as M
V=[t for t in json.load(open('v8_t_polished.json')) if min(t)>0]
L=[t for t in json.load(open('pts_t_polished.json')) if min(t)>0]
random.seed(int(sys.argv[1]) if len(sys.argv)>1 else 1)
N=int(sys.argv[2]) if len(sys.argv)>2 else 60
bad=0; n=0; nt=0
for name,pts,opt in (('V8',V,{'noclauses':True}),('laminar V7',L,{'notarget':True})):
    nr=0
    for t in random.sample(pts,min(N,len(pts))):
        l=C.gaps_from_t(t); T=C.cartesian_tree(l); sh=l[T.root]; l=[x-sh for x in l]; s=C.coords(T,l)
        for w in (1e-3,0.05,0.3):
            for ptail in (0.2,0.5):
                box=[]
                for v in s:
                    b=min(0.0,v+random.uniform(0,w))
                    if random.random()<ptail: box.append((-math.inf,b)); continue
                    if random.random()<ptail: box.append((v-random.uniform(1.0,6.0),b)); continue     # wide
                    a=v-random.uniform(0,w)
                    if a>b: a=b
                    box.append((a,b))
                n+=1; nt+=sum(1 for a,_ in box if a==-math.inf)
                if M.mvf_refute(T,box,opt): nr+=1
    print('%s: refuted %d (must be 0)'%(name,nr)); bad+=nr
print('boxes',n,'tail coords',nt); print('PASS' if bad==0 else 'FAIL')
