"""Cells of the straight figures, classified into rings of outward-pointing triangles (shaded) and the cells
between them (pale twin of the innermost ring they touch)."""
import sys, os, math, numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
from model import Wtab, ext, segs11, segments, apex, base, TRI
from arr import faces
from matplotlib.path import Path
def cells(S):
    pts,F,depth,outer=faces(S)
    ef={}
    for i,f in enumerate(F):
        for e in f['edges']: ef[e]=i
    nb={i:set() for i in range(len(F))}
    for (a,b),i in ef.items():
        j=ef.get((b,a))
        if j is not None and j!=i: nb[i].add(j)
    return pts,F,depth,outer,nb
def meet(p,u,q,v):
    A=np.array([[u[0],-v[0]],[u[1],-v[1]]]); s=np.linalg.solve(A,np.array(q)-np.array(p)); return np.array(p)+s[0]*np.array(u)
def kite(top,ucorner,bot,dcorner):
    """kite of the legs of an upward triangle (apex top, right leg through ucorner) and a downward one
    (apex bot, right leg through dcorner): returns [top, right vertex, bottom, left vertex]"""
    K=meet((0,top),(ucorner[0],ucorner[1]-top),(0,bot),(dcorner[0],dcorner[1]-bot))
    return [np.array([0,top]),K,np.array([0,bot]),np.array([-K[0],K[1]])]
def classify9(y,W,kap):
    """classical nine-triangle figure: ring index 0 (centre) .. 4 (fourteen) from the depth of the cell"""
    pts,F,depth,outer,nb=cells(segments(y,W,kap))
    ring={}; pale={}
    for i in range(len(F)):
        if i==outer: continue
        d=depth[i]
        if d%2==1: ring[i]=(9-d)//2
    for i in range(len(F)):
        if i==outer or i in ring: continue
        rs=[ring[j] for j in nb[i] if j in ring]
        pale[i]=min(rs)
    return pts,F,ring,pale,outer
def kites9(y,W,kap):
    """the nested kites (U3,D5), (U4,D4), (U2,D2), (U1,D1) and the kite of the extension (U0,D0)"""
    C=lambda T: (kap*W[T],y[base[T]])
    K=[]
    for U,D in (('U3','D5'),('U4','D4'),('U2','D2'),('U1','D1')):
        K.append(kite(y[apex[U]],C(U),y[apex[D]],C(D)))
    zD,zU=ext(y,W)
    K.append(kite(zU,C('D1'),zD,C('U1')))
    return K
def classify11(y,W,kap):
    """the extended figure: the cells at odd depth from the outside are the rings, from the new ring (depth 1:
    the tips of U0 and D0 at their apexes and the four triangles at the corners of their bases) inward to the
    central triangle (depth 11); every other cell takes a pale tint of the innermost ring it touches"""
    S=segs11(y,W,kap); pts,F,depth,outer,nb=cells(S)
    zD,zU=ext(y,W)
    KK=kite(zU,(kap*W['D1'],y[7]),zD,(kap*W['U1'],y[3]))
    ring={}; pale={}
    for i in range(len(F)):
        if i!=outer and depth[i]%2==1: ring[i]=(11-depth[i])//2
    for i in range(len(F)):
        if i==outer or i in ring: continue
        pale[i]=min(ring[j] for j in nb[i] if j in ring)
    return pts,F,ring,pale,outer,KK
