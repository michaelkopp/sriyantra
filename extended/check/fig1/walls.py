"""Straight extended figures with walls: triangles carried out beyond the fifth kite, as in Section 8.
Positions on each hand, outermost first: downward D3 (y9), D2 (y8), D1 (y7) read down from the top;
upward U3 (y1), U2 (y2), U1 (y3) read up from the bottom."""
import sys, os, json, numpy as np
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from model import Wtab, ext, apex, base, TRI
DPOS=['D3','D2','D1']; UPOS=['U3','U2','U1']
def build(y,W,kap,wD=(),wU=()):
    """returns {triangle: (apex level, base level, half-width at the base)} for the nine triangles and U0, D0"""
    zD,zU=ext(y,W)
    T={t:(y[apex[t]],y[base[t]],kap*W[t]) for t in TRI}
    def carry(fam,walls,tip,zap,through):
        # walls carried out, outermost first; returns the base level left for the new triangle
        lev=tip
        for t in (DPOS if fam=='D' else UPOS):
            if t in walls:
                a,b,w=T[t]; old=b
                T[t]=(a,lev,w*(lev-a)/(b-a)); lev=old
        return lev
    bD=carry('D',wD,y[10],zD,None)     # the base left for D0
    bU=carry('U',wU,y[0],zU,None)      # the base left for U0
    # D0: apex zD, legs through the corners of U1 and U3 (on the line of the upward corners), base at bD
    sD=(kap*W['U1']-kap*W['U3'])/(y[3]-y[1])          # dx/dy of the line of upward corners
    T['D0']=(zD,bD,sD*(bD-zD))
    sU=(kap*W['D3']-kap*W['D1'])/(y[9]-y[7])          # dx/dy of the line of downward corners (negative)
    T['U0']=(zU,bU,sU*(bU-zU))
    return T
def segs(T):
    S=[]; lab=[]
    for t,(a,b,w) in T.items():
        S+=[((-w,b),(w,b)),((0,a),(w,b)),((0,a),(-w,b))]; lab+=[(t,'base'),(t,'R'),(t,'L')]
    return S,lab
