"""Soundness of the split trinomials in mvf.py: boxes around points of V8 and of the laminar V7 whose coordinates below -1
are made half-infinite or wide, so that trinomials whose terms differ on the u-coordinates occur; no box may be
refuted and the contracted box must contain the point."""
import json, math, sys, random
sys.path.insert(0,'.')
import cone as C, mvf as M
V=[t for t in json.load(open('v8_t_polished.json')) if min(t)>0]
L=[t for t in json.load(open('pts_t_polished.json')) if min(t)>0]
random.seed(11); nx=0; ns=0; bad=0
for name,pts,opt in (('V8',V,{'noclauses':True}),('laminar V7',L,{'notarget':True})):
    for t in random.sample(pts,40):
        l=C.gaps_from_t(t); T=C.cartesian_tree(l); sh=l[T.root]; l=[x-sh for x in l]; s=C.coords(T,l)
        for w in (0.01,0.1):
            box=[]
            for v in s:
                b=min(0.0,v+random.uniform(0,w))
                if v<-1.0: box.append((-math.inf,b) if random.random()<0.5 else (v-random.uniform(1.0,5.0),b))
                else: box.append((max(v-random.uniform(0,w),-50),b))
            Mo,A,b,lo,hi=M._system(T,box,opt)
            ns+=1; nx+=(len(lo)-len(box)-len(Mo.TL))//2
            if M.mvf_refute(T,box,opt): bad+=1
            ref,nb,_=M.mvf_contract(T,box,4,opt)
            if ref or not all(a<=v<=bb for v,(a,bb) in zip(s,nb)): bad+=1
print('systems',ns,'split trinomials',nx,'failures',bad,'PASS' if bad==0 else 'FAIL')
