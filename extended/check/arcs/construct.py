"""The constructions ARC_REC2 and ARC_REC3 of the viewer (index.html), evaluated in floating point or in interval
arithmetic. Generalised circles g = (a, b, c, d) stand for a(x^2+y^2) + bx + cy + d = 0."""
import json, os, math
from fractions import Fraction
from iv import I, cos_sin
HERE=os.path.dirname(os.path.abspath(__file__))
REC=json.load(open(os.path.join(HERE,'recipes.json')))
AXIS=(0,1,0,0); RIM=(1,0,0,-1)
def det3(m):
    return (m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])-m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])
            +m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]))
def through3(P,Q,R):
    rows=[[p[0]*p[0]+p[1]*p[1],p[0],p[1],1] for p in (P,Q,R)]
    g=[]
    for j in range(4):
        m=[[r[k] for k in range(4) if k!=j] for r in rows]
        g.append(det3(m)*(1 if j%2==0 else -1))
    return tuple(g)
def axis2(P,Q):
    r1=[P[0]*P[0]+P[1]*P[1],P[1],1]; r2=[Q[0]*Q[0]+Q[1]*Q[1],Q[1],1]
    v=(r1[1]*r2[2]-r1[2]*r2[1], r1[2]*r2[0]-r1[0]*r2[2], r1[0]*r2[1]-r1[1]*r2[0])
    return (v[0],0,v[1],v[2])
def mid(x): return x.mid() if isinstance(x,I) else x
def mag(x): return abs(mid(x))
def meet_both(g1,g2,allow_none=False):
    """both meeting points of two generalised circles (the formulas of the viewer's gcMeet); with allow_none,
    an empty list when the discriminant is certainly negative"""
    a1,b1,c1,e1=g1; a2,b2,c2,e2=g2
    B=b1*a2-b2*a1; C=c1*a2-c2*a1; E=e1*a2-e2*a1
    if mag(C)>=mag(B):
        u=-B/C if not (isinstance(B,int) and B==0) else 0; v=-E/C if not (isinstance(E,int) and E==0) else 0
        A=a1*(1+u*u); Bq=b1+c1*u+2*a1*u*v; Cq=a1*v*v+c1*v+e1
        roots=quad(A,Bq,Cq,allow_none); return [(x,u*x+v) for x in roots]
    else:
        u=-C/B; v=-E/B
        A=a1*(1+u*u); Bq=c1+b1*u+2*a1*u*v; Cq=a1*v*v+b1*v+e1
        roots=quad(A,Bq,Cq,allow_none); return [(u*y+v,y) for y in roots]
def quad(A,B,C,allow_none=False):
    d=B*B-4*A*C
    if isinstance(d,I):
        if allow_none and d.neg(): return []
        if not d.pos(): raise ValueError('discriminant not certainly positive: %r'%d)
        s=d.sqrt()
    else:
        if d<0: raise ValueError('no meeting'); 
        s=math.sqrt(d)
    return [(-B+s)/(2*A),(-B-s)/(2*A)]
def pick(pts,hint):
    """the meeting point nearest the hint; in intervals the choice must be decisive"""
    ds=[(p[0]-hint[0])*(p[0]-hint[0])+(p[1]-hint[1])*(p[1]-hint[1]) for p in pts]
    if isinstance(ds[0],I):
        k=0 if ds[0].mid()<ds[1].mid() else 1
        if not ds[k].hi<ds[1-k].lo: raise ValueError('branch choice not decisive')
        return pts[k]
    return pts[0] if ds[0]<ds[1] else pts[1]
def mirp(p): return (-p[0],p[1])
def mirg(g): return (g[0],-g[1],g[2],g[3])
def values(key,exact=False,digits=None):
    """the parameter values: the viewer's defaults, or those rounded to `digits` decimals as Fractions"""
    v=REC[key]['vals']
    if digits is not None: return [Fraction(round(x*10**digits),10**digits) for x in v]
    return v
def build(key,vals,interval=False):
    R=REC[key]; N=[None]*R['nnodes']; G=[None]*R['nsides']; sv={}; k=0
    num=(lambda x: I(Fraction(x)) if not isinstance(x,I) else x) if interval else float
    for ni,kind in R['par']:
        if kind==0: sv[ni]=(num(vals[k]),num(vals[k+1])); k+=2
        elif kind==1: sv[ni]=(0,num(vals[k])); k+=1
        else:
            if interval: c,s=cos_sin(vals[k])
            else: c,s=math.cos(vals[k]),math.sin(vals[k])
            sv[ni]=(c,s); k+=1
    for op in R['ops']:
        kind,i=op[0],op[1]
        if kind=='seed': N[i]=sv[i]; continue
        if kind=='pole': N[i]=(0,1 if R['hint'][i][1]>0 else -1); continue
        if kind in ('s2','s3'):
            pts=[mirp(N[op[a]]) if op[a+1] else N[op[a]] for a in range(2,len(op),2)]
            G[i]=axis2(*pts) if kind=='s2' else through3(*pts); continue
        gs=[mirg(G[op[a]]) if op[a+1] else G[op[a]] for a in range(2,len(op),2)]
        other={'nx':None,'na':AXIS,'nr':RIM}[kind]
        g2=gs[1] if kind=='nx' else other
        N[i]=pick(meet_both(gs[0],g2),R['hint'][i])
    return N,G
def incidences(key):
    """the node-on-side incidences each step realises: (step, node, side, mirrored)"""
    R=REC[key]; out=[]
    for st,op in enumerate(R['ops']):
        kind,i=op[0],op[1]
        if kind in ('s2','s3'):
            for a in range(2,len(op),2): out.append((st,op[a],i,op[a+1]))
        elif kind=='nx':
            for a in (2,4): out.append((st,i,op[a],op[a+1]))
        elif kind in ('na','nr'): out.append((st,i,op[2],op[3]))
    return out
