"""Soundness of the contraction (mvf.mvf_contract): boxes (some with half-infinite coordinates) around points of V8
(target, no clauses) and of the laminar V7 (clauses, no target).  The box must not be refuted and the contracted box
must still contain the point."""
import json, math, sys, random, time
sys.path.insert(0,'.')
import cone as C, mvf as M
V=[t for t in json.load(open('v8_t_polished.json')) if min(t)>0]
L=[t for t in json.load(open('pts_t_polished.json')) if min(t)>0]
random.seed(int(sys.argv[1]) if len(sys.argv)>1 else 1)
N=int(sys.argv[2]) if len(sys.argv)>2 else 20
bad=0; n=0; st={}; t0=time.time()
for name,pts,opt in (('V8',V,{'noclauses':True}),('laminar V7',L,{'notarget':True})):
    nr=0; nout=0
    for t in random.sample(pts,min(N,len(pts))):
        l=C.gaps_from_t(t); T=C.cartesian_tree(l); sh=l[T.root]; l=[x-sh for x in l]; s=C.coords(T,l)
        for w in (0.02,0.1,0.4):
            for ptail in (0.0,0.3):
                box=[]
                for v in s:
                    b=min(0.0,v+random.uniform(0,w))
                    if random.random()<ptail: box.append((-math.inf,b)); continue
                    if random.random()<ptail: box.append((v-random.uniform(1.0,6.0),b)); continue     # wide
                    a=v-random.uniform(0,w)
                    if a>b: a=b
                    box.append((a,b))
                n+=1
                ref,nb,_=M.mvf_contract(T,box,4,opt,st)
                if ref: nr+=1
                elif not all(a<=v<=b for v,(a,b) in zip(s,nb)): nout+=1
    print('%s: refuted %d, point outside the contracted box %d (both must be 0)'%(name,nr,nout)); bad+=nr+nout
print('boxes',n,'cuts',st.get('cut',0),'%.0fs'%(time.time()-t0)); print('PASS' if bad==0 else 'FAIL')
