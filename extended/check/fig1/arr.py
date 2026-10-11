"""Planar arrangement of the nine straight triangles: atomic faces, their areas and their depth from the outside."""
import math, numpy as np
apex={'D1':0,'D2':1,'D3':4,'D4':2,'D5':3,'U1':10,'U2':9,'U3':7,'U4':8}
base={'D1':7,'D2':8,'D3':9,'D4':6,'D5':5,'U1':3,'U2':2,'U3':1,'U4':4}
TRI=list(apex)
def Wtab(y):
    R=lambda b,a: y[a]-y[b]
    return {'D1': R(0,7),'D2': R(1,8)*R(0,2)*R(3,7)/(R(2,7)*R(1,3)),'D3': R(4,9)*R(0,2)*R(5,7)/(R(2,7)*R(4,5)),
 'D4': R(2,6)*R(0,2)*R(4,7)/(R(2,7)*R(2,4)),'D5': R(0,2)*R(5,7)*R(4,6)*R(5,8)/(R(2,7)*R(4,5)*R(6,8)),
 'U1': R(3,10)*R(0,2)*R(5,7)*R(4,8)/(R(2,7)*R(4,5)*R(8,10)),'U2': R(2,9)*R(0,3)/R(3,9),'U3': R(1,7)*R(0,2)/R(2,7),
 'U4': R(4,8)*R(0,2)*R(5,7)*R(4,6)/(R(2,7)*R(4,5)*R(6,8))}
def segments(y,W,kappa=1.0):
    S=[]
    for T in TRI:
        a,b,w=y[apex[T]],y[base[T]],kappa*W[T]
        S+= [((-w,b),(w,b)),((0,a),(w,b)),((0,a),(-w,b))]
    return S
def faces(S,tol=1e-10):
    n=len(S); params=[[0.0,1.0] for _ in range(n)]
    P=[np.array(s[0],float) for s in S]; Q=[np.array(s[1],float) for s in S]
    for i in range(n):
        for j in range(i+1,n):
            d1=Q[i]-P[i]; d2=Q[j]-P[j]; den=d1[0]*d2[1]-d1[1]*d2[0]
            if abs(den)<1e-14: continue
            r=P[j]-P[i]; t=(r[0]*d2[1]-r[1]*d2[0])/den; u=(r[0]*d1[1]-r[1]*d1[0])/den
            if -tol<=t<=1+tol and -tol<=u<=1+tol:
                params[i].append(min(1,max(0,t))); params[j].append(min(1,max(0,u)))
    pts=[]; ids={}
    def vid(p):
        key=(round(p[0]/1e-9),round(p[1]/1e-9))
        for dk in ((0,0),):
            pass
        if key in ids: return ids[key]
        # search neighbours for snapping
        for dx in (-1,0,1):
            for dy in (-1,0,1):
                k=(key[0]+dx,key[1]+dy)
                if k in ids: return ids[k]
        ids[key]=len(pts); pts.append(p); return ids[key]
    adj={}
    for i in range(n):
        ts=sorted(set(params[i])); vs=[]
        for t in ts:
            v=vid(P[i]+t*(Q[i]-P[i]))
            if not vs or vs[-1]!=v: vs.append(v)
        for a,b in zip(vs,vs[1:]):
            if a==b: continue
            adj.setdefault(a,set()).add(b); adj.setdefault(b,set()).add(a)
    pts=np.array(pts)
    ang={v:sorted(adj[v],key=lambda w: math.atan2(pts[w][1]-pts[v][1],pts[w][0]-pts[v][0])) for v in adj}
    used=set(); F=[]; edgeface={}
    for u in adj:
        for v in adj[u]:
            if (u,v) in used: continue
            cyc=[]; a,b=u,v
            while (a,b) not in used:
                used.add((a,b)); cyc.append((a,b))
                L=ang[b]; k=L.index(a); c=L[(k-1)%len(L)]   # next edge: most clockwise turn
                a,b=b,c
            area=0.5*sum(pts[x][0]*pts[y_][1]-pts[y_][0]*pts[x][1] for x,y_ in cyc)
            fi=len(F); F.append({'edges':cyc,'area':area,'verts':[x for x,_ in cyc]})
            for e in cyc: edgeface[e]=fi
    outer=min(range(len(F)),key=lambda i:F[i]['area'])   # the one with the most negative area
    # depth by BFS on the dual graph
    from collections import deque
    depth={outer:0}; dq=deque([outer])
    while dq:
        f=dq.popleft()
        for (a,b) in F[f]['edges']:
            g=edgeface.get((b,a))
            if g is not None and g not in depth: depth[g]=depth[f]+1; dq.append(g)
    return pts,F,depth,outer
def shells(y,W,kappa=1.0,amin=1e-9):
    pts,F,depth,outer=faces(segments(y,W,kappa))
    sh={}
    for i,f in enumerate(F):
        if i==outer: continue
        if f['area']<amin: continue
        sh.setdefault(depth[i],[]).append(f['area'])
    return sh
def corners(P,tol=1e-7):
    """number of geometric corners of a polygon (collinear vertices dropped)"""
    n=len(P); c=0
    for i in range(n):
        a,b,d=P[i-1],P[i],P[(i+1)%n]
        u=b-a; v=d-b
        cr=u[0]*v[1]-u[1]*v[0]
        if abs(cr)>tol*np.linalg.norm(u)*np.linalg.norm(v): c+=1
    return c
