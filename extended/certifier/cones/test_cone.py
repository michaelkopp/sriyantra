"""Soundness checks for cone.py: enclosures contain the true values; boxes around V8 points are never refuted
without clauses; boxes around laminar V7 points are never refuted without the target; also half-infinite boxes."""
import json, math, sys, random
sys.path.insert(0,'.')
import cone as C, bracketlp as B
V=json.load(open('v8_t.json')); L=json.load(open('pts_t.json'))
V=[t for t in V if min(t)>0]; L=[t for t in L if min(t)>0]
random.seed(2)
def true_L(t):
    lam=[None]+[math.log(x) for x in t]
    def lp(a,b):
        if a>b: a,b=b,a
        s=sum(lam[j] for j in range(b,11))
        if a==0: return s
        return s+math.log(-math.expm1(sum(lam[j] for j in range(a,b))))
    out={}
    for (p,q),k in C.PAIRIDX.items():
        a,b=C.ORDER[p],C.ORDER[q]
        out[k]= -lp(b if a==12 else a,8) if (a==12 or b==12) else lp(a,b)-lp(a,8)-lp(b,8)
    return out
bad=0; nenc=0
def box_around(T,s,w,mode):
    bx=[]
    for v in s:
        lo,hi=v-w,min(v+w,0.0)
        if mode==1 and v<-12: lo=-math.inf
        if mode==2 and random.random()<0.3: lo=-math.inf
        bx.append((lo,hi))
    return bx
for name,pts,opt in (('V8',V,{'noclauses':True}),('laminar V7',L,{'notarget':True})):
    sel=random.sample(pts,min(60,len(pts)))
    nref=0; niv=0
    for t in sel:
        l=C.gaps_from_t(t); T=C.cartesian_tree(l); sh=l[T.root]; l=[x-sh for x in l]; s=C.coords(T,l)
        tl={k:v-sh for k,v in true_L(t).items()}
        for w,mode in ((1e-3,0),(1e-2,1),(1e-2,2)):
            E=C.Enc(T,box_around(T,s,w,mode))
            for k,v in tl.items():
                a,b=E.absolute(k)
                if not (a<=v+1e-9 and v<=b+1e-9): nenc+=1
            iv=C.interval_refute(E,notarget=opt.get('notarget',False),noclauses=opt.get('noclauses',False))
            if iv: niv+=1
            if C.lp_refute(E,[],opt): nref+=1
            if C.lp_refute_sx(E,[],opt): nref+=1
            if C.lp_refute_sx2(E,[],opt): nref+=1
            if C.lp_refute_sx3(E,[],opt): nref+=1
    print('%s: boxes refuted by interval %d, by LP %d (must be 0)'%(name,niv,nref)); bad+=niv+nref
print('enclosure failures %d (must be 0)'%nenc); bad+=nenc
print('PASS' if bad==0 else 'FAIL')
