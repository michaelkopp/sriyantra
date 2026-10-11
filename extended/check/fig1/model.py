import os
"""Straight extended figure: ladder y0..y10 (y0=-1, y10=1), widths, extension, residuals, margins, rings."""
import math, numpy as np, json, sys
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
from arr import faces, segments, Wtab, apex, base, corners, TRI
DN=['D1','D2','D3','D4','D5']; UPF=['U1','U2','U3','U4']
def yy(y,i,j): return y[j]-y[i]
def ext(y,W):
    zU=(W['D1']*y[9]-W['D3']*y[7])/(W['D1']-W['D3'])
    zD=(W['U1']*y[1]-W['U3']*y[3])/(W['U1']-W['U3'])
    return zD,zU
def _lg(x):
    return math.log(x) if (x==x and x>0) else float('nan')
def resid(y,W=None):
    if W is None: W=Wtab(y)
    Y=lambda i,j: y[j]-y[i]
    E1=_lg(Y(1,4)*Y(3,7)*Y(4,5)*Y(6,8)/(Y(1,3)*Y(4,6)*Y(4,8)*Y(5,7)))
    E2=_lg(Y(0,3)*Y(2,4)*Y(2,7)*Y(6,9)/(Y(0,2)*Y(2,6)*Y(3,9)*Y(4,7)))
    E3=_lg(Y(0,3)*Y(2,7)*Y(4,5)*Y(7,9)/(Y(0,2)*Y(3,9)*Y(4,7)*Y(5,7)))
    E4=_lg(Y(1,3)*Y(4,8)*Y(5,7)*Y(7,10)/(Y(1,7)*Y(3,7)*Y(4,5)*Y(8,10)))
    CD=_lg(W['D2']*Y(7,9)/(W['D1']*Y(8,9)+W['D3']*Y(7,8)))
    CU=_lg(W['U2']*Y(1,3)/(W['U3']*Y(2,3)+W['U1']*Y(1,2)))
    return dict(E1=E1,E2=E2,E3=E3,E4=E4,CD=CD,CU=CU)
def w(y,W,T,h):
    a,b=y[apex[T]],y[base[T]]; return W[T]*(h-a)/(b-a)
DORD=['D1','D2','D4','D5','D3']   # apexes y0<y1<y2<y3<y4
UORD=['U1','U2','U4','U3']         # apexes y10>y9>y8>y7
def _sl(x):
    return math.log(x) if (x==x and x>0) else -10.0
def margins(y,W=None):
    if W is None: W=Wtab(y)
    M={}
    for i,A in enumerate(DORD):
        for B in DORD[i+1:]:
            h=min(y[base[A]],y[base[B]])
            M['nc %s%s'%(A,B)]=_sl(w(y,W,A,h)/w(y,W,B,h))
    for i,A in enumerate(UORD):
        for B in UORD[i+1:]:
            h=max(y[base[A]],y[base[B]])
            M['nc %s%s'%(A,B)]=_sl(w(y,W,A,h)/w(y,W,B,h))
    lU=lambda t: W['D1']+(W['D3']-W['D1'])*(t-y[7])/(y[9]-y[7])
    lD=lambda t: W['U3']+(W['U1']-W['U3'])*(t-y[1])/(y[3]-y[1])
    for T in UPF: M['corner '+T]=_sl(lU(y[base[T]])/W[T])
    for T in DN:  M['corner '+T]=_sl(lD(y[base[T]])/W[T])
    zD,zU=ext(y,W)
    M['tip U']=math.log((zU-y[10])/(y[10]-y[0])) if zU>y[10] else -99
    M['tip D']=math.log((y[0]-zD)/(y[10]-y[0])) if zD<y[0] else -99
    return M
def gap(y,W=None):
    if W is None: W=Wtab(y)
    return w(y,W,'D2',y[4])/W['U4']-1
def segs11(y,W,kappa=1.0):
    zD,zU=ext(y,W)
    S=segments(y,W,kappa)
    s=(y[0]-zU)/(y[7]-zU); xb=kappa*W['D1']*s
    S+=[((-xb,y[0]),(xb,y[0])),((0,zU),(xb,y[0])),((0,zU),(-xb,y[0]))]
    s=(y[10]-zD)/(y[3]-zD); xb=kappa*W['U1']*s
    S+=[((-xb,y[10]),(xb,y[10])),((0,zD),(xb,y[10])),((0,zD),(-xb,y[10]))]
    return S
EXPECT={1:6,2:14,3:14,4:10,5:10,6:8,7:10,8:8,9:8,10:3,11:1}
RINGS=[2,3,5,7,9]   # new ring, 14, 10, 10, 8
def rings(y,W=None,kappa=1.0,with11=True):
    if W is None: W=Wtab(y)
    S=segs11(y,W,kappa) if with11 else segments(y,W,kappa)
    pts,F,depth,outer=faces(S)
    sh={}
    for i,f in enumerate(F):
        if i==outer: continue
        sh.setdefault(depth[i],[]).append(f['area'])
    return sh
def ringcv(sh):
    out={}
    for k in RINGS:
        a=np.array(sh.get(k,[1.0])); out[k]=float(a.std()/a.mean())
    return out
