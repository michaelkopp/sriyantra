"""Certify the two constructed figures at the corners of Table 4 (the viewer's defaults rounded to four decimals):
  (1) the construction is defined: every meeting it takes exists and its branch is chosen decisively;
  (2) its node-on-side incidences are exactly those its steps realise, each once, with no accidental one;
  (3) the eleven triangles are laminar: no leg meets a leg of another triangle of its family.
Everything is evaluated in interval arithmetic from the exact rational corners.
Usage: python3 certify.py   (writes incidences_22.csv, incidences_36.csv and certify.log)"""
import sys, os, json, math
from fractions import Fraction
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
from iv import I
from construct import REC, build, values, incidences, mirp, mirg, meet_both, AXIS, RIM
PAT={'ARC_REC2':'22','ARC_REC3':'36'}
def gval(g,p): return g[0]*(p[0]*p[0]+p[1]*p[1])+g[1]*p[0]+g[2]*p[1]+g[3]
def iszero_possible(x): return x.lo<=0<=x.hi
def cross(o,a,b): return (a[0]-o[0])*(b[1]-o[1])-(a[1]-o[1])*(b[0]-o[0])
def sgn(x):
    if x.lo>0: return 1
    if x.hi<0: return -1
    return 0
log=[]
def say(*a):
    s=' '.join(str(x) for x in a); print(s); log.append(s)
allok=True
for key in ('ARC_REC2','ARC_REC3'):
    R=REC[key]; pat=PAT[key]
    vals=values(key,digits=4)
    say('== pattern %s: %d steps, %d nodes, %d sides; corners %s'%(pat,len(R['ops']),R['nnodes'],R['nsides'],[str(v) for v in vals]))
    # (1) the construction in interval arithmetic (raises if a meeting or a branch is not certain)
    N,G=build(key,vals,interval=True)
    wmax=max(max(p[0].w(),p[1].w()) for p in N if isinstance(p[0],I))
    say('  (1) construction defined; widest node enclosure %.1e'%wmax)
    # (2) incidences
    inc=incidences(key)
    names=R['sideNames']; legs=[j for j,n in enumerate(names) if n.endswith('.leg')]
    onaxis=set(op[1] for op in R['ops'] if op[0] in ('na','pole'))|set(ni for ni,k in R['par'] if k==1)
    def norm(node,side,m):
        if not names[side].endswith('.leg') or node in onaxis: m=0     # bases are symmetric; an axis node lies on both legs
        return (node,side,m)
    realised=[norm(n,s,m) for _,n,s,m in inc]
    assert len(set(realised))==len(realised), 'an incidence is realised twice'
    rs=set(realised)
    extra=[]
    for node in range(R['nnodes']):
        P=N[node]
        for side in range(R['nsides']):
            for m in ((0,1) if side in legs else (0,)):
                key3=norm(node,side,m)
                g=mirg(G[side]) if m else G[side]
                v=gval(g,P)
                if key3 in rs: continue
                if iszero_possible(v): extra.append(key3)
    say('  (2) %d node-on-side incidences, each realised by exactly one step; pairs not realised and not excluded: %d'%(len(rs),len(extra)))
    if extra: allok=False; say('      not excluded:',extra[:10])
    with open(os.path.join(HERE,'incidences_%s.csv'%pat),'w') as f:
        f.write('step,move,node,side,mirror\n')
        for st,n,s,m in inc: f.write('%d,%s,%s,%s,%d\n'%(st,R['ops'][st][0],R['nodeNames'][n],names[s],m))
    # (3) laminarity: each leg is the arc of its circle from its apex through its nodes to its corner
    fl,_=build(key,values(key,digits=4),interval=False)
    arcs={}
    for side in legs:
        pts=set()
        for n,s,m in rs:
            if s==side: pts.add((n,m))
        # node positions (right leg): mirrored incidences put the mirror image of the node on the right leg
        P=[(N[n] if m==0 else mirp(N[n]), fl[n] if m==0 else mirp(fl[n])) for n,m in pts]
        ax=[p for p in P if abs(p[1][0])<1e-12]
        assert len(ax)>=1, names[side]
        apex=max(ax,key=lambda p: 0) if len(ax)==1 else None
        if apex is None:   # the leg of U0 or D0 starts at a pole; another axis node may lie on its circle
            apex=[p for p in ax if abs(abs(p[1][1])-1)<1e-12][0] if any(abs(abs(p[1][1])-1)<1e-12 for p in ax) else ax[0]
        g=G[side]; gm=[x.mid() for x in g]
        cx,cy=-gm[1]/(2*gm[0]),-gm[2]/(2*gm[0])
        a0=math.atan2(apex[1][1]-cy,apex[1][0]-cx)
        rel=[]
        for p in P:
            if p is apex: continue
            d=(math.atan2(p[1][1]-cy,p[1][0]-cx)-a0)%(2*math.pi)
            rel.append((d,p))
        # the arc runs the shorter way round from the apex to cover all nodes
        span_pos=max(d for d,_ in rel); span_neg=max(2*math.pi-d for d,_ in rel)
        if span_neg<span_pos: rel=[(2*math.pi-d,p) for d,p in rel]
        rel.sort(key=lambda t:t[0])
        corner=rel[-1][1]; inner=rel[len(rel)//2][1] if len(rel)>1 else None
        assert inner is not None and inner is not corner, names[side]
        arcs[side]=dict(g=g,P=apex[0],Q=corner[0],M=inner[0],n=len(P))
    def mirror_arc(a): return dict(g=mirg(a['g']),P=mirp(a['P']),Q=mirp(a['Q']),M=mirp(a['M']),n=a['n'])
    def off_arc(a,X):
        s=sgn(cross(a['P'],a['Q'],X)); sm=sgn(cross(a['P'],a['Q'],a['M']))
        assert sm!=0
        return s==-sm
    fam=lambda side: names[side][0]
    pairs=0; disjoint=0
    for i,A in enumerate(legs):
        for B in legs[i+1:]:
            if fam(A)!=fam(B): continue
            for mb in (0,1):
                a=arcs[A]; b=arcs[B] if mb==0 else mirror_arc(arcs[B])
                pairs+=1
                ga,gb=a['g'],b['g']
                # the circles: disjoint, or both meeting points off one of the two arcs
                X=meet_both(ga,gb,allow_none=True)      # [] when the two circles certainly do not meet
                if not X: disjoint+=1
                for x in X:
                    if not (off_arc(a,x) or off_arc(b,x)):
                        allok=False; say('      NOT CERTIFIED: %s and %s%s may meet near (%.4f, %.4f)'%(names[A],names[B],"'" if mb else '',x[0].mid(),x[1].mid()))
    say('  (3) laminarity: %d pairs of legs of one family (a right leg with a right or a left leg): %d on circles that do not meet,'%(pairs,disjoint)
        +' %d with both meeting points certainly off one of the two arcs'%(pairs-disjoint) if allok else '  (3) laminarity: FAILED')
say('ALL CERTIFIED' if allok else 'NOT ALL CERTIFIED')
open(os.path.join(HERE,'certify.log'),'w').write('\n'.join(log)+'\n')
