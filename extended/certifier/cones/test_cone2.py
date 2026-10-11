"""Stronger soundness check: whole-cone root boxes and random sub-boxes containing known points.
A box containing a V8 point must never be refuted without clauses; a box containing a laminar V7 point never
without the target.  Boxes: the point's whole Cartesian-tree cone, and random boxes around it (half-infinite)."""
import json, math, sys, random
sys.path.insert(0,'.')
import cone as C
V=[t for t in json.load(open('v8_t.json')) if min(t)>0]
L=[t for t in json.load(open('pts_t.json')) if min(t)>0]
random.seed(int(sys.argv[1]) if len(sys.argv)>1 else 3)
bad=0; n=0
def rbox(s):
    bx=[]
    for v in s:
        r=random.random()
        if r<0.3: bx.append((-math.inf,0.0))
        elif r<0.6: bx.append((v-random.uniform(0,5),min(0.0,v+random.uniform(0,5))))
        else:
            lo=v-random.uniform(0,0.5); hi=min(0.0,v+random.uniform(0,0.5)); bx.append((lo if random.random()<0.8 else -math.inf,hi))
    return bx
for name,pts,opt in (('V8',V,{'noclauses':True}),('laminar V7',L,{'notarget':True})):
    sel=random.sample(pts,min(int(sys.argv[2]) if len(sys.argv)>2 else 40,len(pts))); nr=0; nin=0
    for t in sel:
        l=C.gaps_from_t(t); T=C.cartesian_tree(l); sh=l[T.root]; l=[x-sh for x in l]; s=C.coords(T,l)
        for bx in ([(-math.inf,0.0)]*10, rbox(s), rbox(s)):
            assert all(a<=v<=b for v,(a,b) in zip(s,bx))
            E=C.Enc(T,bx); n+=1
            if C.interval_refute(E,notarget=opt.get('notarget',False),noclauses=opt.get('noclauses',False)): nr+=1
            if C.lp_refute_sx2(E,[],opt): nr+=1
            # with the rows of a cell that contains the point
            eta=[sum(x*C.gaps_from_t and 0 for x in []) for _ in range(0)]
            cs=C.cells_hit(E)
            for c in cs:
                # only cells that contain the point itself
                tl=None
            for c in C.CELLS:
                Ept=C.Enc(T,[(v,v) for v in s]); e=[Ept.form(r) for r in C.ETA]
                if all(e[m][1]>=C.EDGES[c[m]] and e[m][0]<=C.EDGES[c[m]+1] for m in range(4)):
                    nin+=1
                    if C.lp_refute_sx2(E,C.cell_rows(c),opt): nr+=1
    print('%s: refuted %d (must be 0); boxes with the point in a listed cell: %d'%(name,nr,nin)); bad+=nr
print('boxes tested',n); print('PASS' if bad==0 else 'FAIL')
