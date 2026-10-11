"""Exact checks of the straight-figure algebra used in the paper (Fractions only)."""
from fractions import Fraction as F
import random
apex={'D1':0,'D2':1,'D3':4,'D4':2,'D5':3,'U1':10,'U2':9,'U3':7,'U4':8}
base={'D1':7,'D2':8,'D3':9,'D4':6,'D5':5,'U1':3,'U2':2,'U3':1,'U4':4}
def Wtab(y):
    R=lambda b,a: y[a]-y[b]
    return {
 'D1': R(0,7),
 'D2': R(1,8)*R(0,2)*R(3,7)/(R(2,7)*R(1,3)),
 'D3': R(4,9)*R(0,2)*R(5,7)/(R(2,7)*R(4,5)),
 'D4': R(2,6)*R(0,2)*R(4,7)/(R(2,7)*R(2,4)),
 'D5': R(0,2)*R(5,7)*R(4,6)*R(5,8)/(R(2,7)*R(4,5)*R(6,8)),
 'U1': R(3,10)*R(0,2)*R(5,7)*R(4,8)/(R(2,7)*R(4,5)*R(8,10)),
 'U2': R(2,9)*R(0,3)/R(3,9),
 'U3': R(1,7)*R(0,2)/R(2,7),
 'U4': R(4,8)*R(0,2)*R(5,7)*R(4,6)/(R(2,7)*R(4,5)*R(6,8))}
def hw(W,y,T,lev):
    a,b=y[apex[T]],y[base[T]]
    return W[T]*(lev-a)/(b-a)
# Chiodo's twelve triples (t',t'',t''') -> (down leg, base, up leg), straight labelling
C=[('D3','D2','U1'),('D2','D1','U1'),('D3','D1','U2'),('D3','D4','U4'),('D3','D5','U3'),('D4','U4','U3'),
   ('D2','U1','U3'),('D1','U1','U2'),('D1','U2','U3'),('D4','D4','U2'),('D5','D5','U4'),('D2','U4','U4')]
def quart(y):
    R=lambda b,a: y[a]-y[b]
    Q={'Q9':(R(1,3)*R(5,7)*R(4,6)*R(4,8), R(4,5)*R(6,8)*R(3,7)*R(1,4)),
       'Q21':(R(2,7)*R(2,4)*R(0,3)*R(6,9), R(3,9)*R(0,2)*R(4,7)*R(2,6)),
       'Q43':(R(2,7)*R(4,5)*R(0,3)*R(7,9), R(3,9)*R(0,2)*R(5,7)*R(4,7)),
       'Q59':(R(1,3)*R(5,7)*R(4,8)*R(7,10), R(4,5)*R(8,10)*R(3,7)*R(1,7))}
    return Q
random.seed(3)
tree=set(); res={}
for trial in range(20):
    v=sorted(random.sample(range(1,10000),11)); y=[F(x) for x in v]
    W=Wtab(y); Q=quart(y)
    for k,(d,b,u) in enumerate(C):
        lev=y[base[b]]
        r=hw(W,y,d,lev)-hw(W,y,u,lev)
        if r==0: tree.add(k)
        res.setdefault(k,[]).append((r,lev,W,Q,y))
print('identically satisfied by the table:',sorted(tree), [C[k] for k in sorted(tree)])
# for the others: ratio w_d/w_u at the level, compared with quartic ratios P/N
for k in range(12):
    if k in tree: continue
    d,b,u=C[k]
    vals=[]
    for (r,lev,W,Q,y) in res[k]:
        ratio=hw(W,y,d,lev)/hw(W,y,u,lev)
        m=[n for n,(P,N) in Q.items() if ratio in (P/N,N/P)]
        vals.append(tuple(m))
    print(k,C[k],'ratio w_down/w_up equals a quartic ratio:',set(vals))
