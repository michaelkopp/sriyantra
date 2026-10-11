import json, math, random, os
D=json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','certifier','cones','tri.json'))); P=[tuple(p) for p in D['pairs']]
def L(vec, p):
    s=0.0
    for (a,b),e in zip(P,vec):
        if e==0: continue
        if 12 in (a,b): continue      # zU at infinity: weight zero, cancels
        s+=e*math.log(abs(p[a]-p[b]))
    return s
def pts(t):  # t[1..10]
    p=[0.0]*13; p[11]=1.0
    for k in range(10,0,-1): p[k]=t[k]*p[k+1]
    p[0]=0.0; p[12]=float('inf'); return p
def g(t,a,b):
    pr=1.0
    for k in range(a,b): pr*=t[k]
    return 1-pr
M=D['mono']
random.seed(1)
npos=[(1,8),(2,4),(2,5),(3,7),(4,10),(5,8),(7,9),(8,11),(9,10)]
nneg=[(1,4),(2,8),(2,9),(3,5),(4,8),(5,7),(7,10),(8,10),(9,11)]
for trial in range(5):
    t=[None]+[random.uniform(0.05,0.95) for _ in range(10)]
    p=pts(t)
    m1,m2=math.exp(L(M['DOWN.m1'],p)),math.exp(L(M['DOWN.m2'],p))
    Phi=L(M['target'],p)-L(M['Q21'],p)+L(M['Q59'],p)+math.log(m1+m2)
    R=m2/m1
    closed=math.log(1+R)+sum(math.log(g(t,a,b)) for a,b in npos)-sum(math.log(g(t,a,b)) for a,b in nneg)
    print('%.15f %.15f'%(Phi,closed))

print('--- check R on apexU')
def f_apex(t): return L(M['apexU'],pts(t))
for trial in range(4):
    t=[None]+[random.uniform(0.05,0.95) for _ in range(10)]
    # solve apexU = 0 for t9 by bisection on (1e-9, 1-1e-9)
    lo,hi=1e-9,1-1e-9
    t[9]=lo; flo=f_apex(t); t[9]=hi; fhi=f_apex(t)
    if flo*fhi>0: print('no sign change', flo, fhi); continue
    for _ in range(200):
        mid=(lo+hi)/2; t[9]=mid; fm=f_apex(t)
        if fm*flo>0: lo,flo=mid,fm
        else: hi=mid
    p=pts(t); m1,m2=math.exp(L(M['DOWN.m1'],p)),math.exp(L(M['DOWN.m2'],p))
    print('%.12f %.12f'%(m2/m1, t[9]*(1-t[8])/(1-t[9])))
