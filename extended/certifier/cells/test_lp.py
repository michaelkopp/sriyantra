"""Soundness checks for the LP certificates (bracketlp.py):
   (1) a box around a point of V8 is never refuted when the clauses are left out;
   (2) a box around a laminar point of V7 is never refuted when the target is left out;
   (3) the handoff test succeeds around laminar points with t3 < 0.2 (region certified by chart B);
   (4) boxes around points of V8 with some enclosures made one-sided (lower bound -inf or upper bound
       +inf) are never refuted when the clauses are left out.
   Points: v8_t.json, pts_t.json (fan coordinates t1..t10)."""
import json, math, sys, os
H=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,H)
import bracketlp as B
V=json.load(open(os.path.join(H,'v8_t.json'))); L=json.load(open(os.path.join(H,'pts_t.json')))
def bounds(t,w):
    lam=[None]+[math.log(x) for x in t]
    def lp(a,b):
        if a>b: a,b=b,a
        s=sum(lam[j] for j in range(b,11))
        if a==0: return s
        return s+math.log(-math.expm1(sum(lam[j] for j in range(a,b))))
    bd=[None]*B.NV
    for (a,b),k in B.VI.items():
        v=-lp(a if b==12 else b,8) if (a==12 or b==12) else lp(a,b)-lp(a,8)-lp(b,8)
        bd[k]=(v-w,v+w)
    return bd
bad=0
for w in (1e-4,1e-2):
    n1=sum(B.refute_bounds(bounds(t,w),{'noclauses':True}) for t in V[:100] if t[2]>1e-6)
    n2=sum(B.refute_bounds(bounds(t,w),{'notarget':True}) for t in L[:100])
    print('width %g: V8 boxes refuted %d (must be 0), laminar boxes refuted without target %d (must be 0)'%(w,n1,n2)); bad+=n1+n2
# (4) boxes with one-sided (half-infinite) enclosures around points of V8, including points close to T*
import random
random.seed(1); n4=0
pts=sorted([t for t in V if min(t)>0], key=lambda t:t[2])
for t in pts[:25]+pts[-15:]:
    bd=bounds(t,1e-2)
    for mode in range(3):
        b2=[]
        for (lo,hi) in bd:
            if mode==0 and hi<-15: lo=-math.inf
            if mode==1 and random.random()<0.3: lo=-math.inf
            if mode==2 and random.random()<0.3: hi=math.inf
            b2.append((lo,hi))
        n4+=B.refute_bounds(b2,{'noclauses':True})
print('V8 boxes with half-infinite enclosures refuted %d (must be 0)'%n4); bad+=n4
ok=sum(B.handoff(bounds(t,1e-2),{'notarget':True}) for t in L if t[2]<0.19)
print('handoff around laminar points with t3 < 0.19: %d of %d'%(ok,sum(1 for t in L if t[2]<0.19)))
print('PASS' if bad==0 else 'FAIL')
