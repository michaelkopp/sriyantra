"""The executor again, with two changes that make the algebra visible:

  * no normalisation of the circle vectors. `unit()` divides by a square root
    that has nothing to do with the geometry - the circle is projective - and
    it would destroy the rationality we are trying to measure.
  * the branch at one named intersection can be flipped, which conjugates the
    single t-dependent square root and lets A and B be separated:
        r  = A + B*sqrt(C)      r~ = A - B*sqrt(C)
        A  = (r + r~)/2         r*r~ = A^2 - B^2 C  is rational in t
"""
from mpmath import mp, mpf, mpc, sqrt


def through3(p1, p2, p3):
    (x1, y1), (x2, y2), (x3, y3) = p1, p2, p3
    r1, r2, r3 = x1*x1 + y1*y1, x2*x2 + y2*y2, x3*x3 + y3*y3
    d3 = lambda a1,b1,c1,a2,b2,c2,a3,b3,c3: (a1*(b2*c3-b3*c2) - b1*(a2*c3-a3*c2)
                                             + c1*(a2*b3-a3*b2))
    one = mpf(1)
    return [ d3(x1,y1,one, x2,y2,one, x3,y3,one),
            -d3(r1,y1,one, r2,y2,one, r3,y3,one),
             d3(r1,x1,one, r2,x2,one, r3,x3,one),
            -d3(r1,x1,y1, r2,x2,y2, r3,x3,y3)]


def through_axis2(p1, p2):
    (x1, y1), (x2, y2) = p1, p2
    r1, r2 = x1*x1 + y1*y1, x2*x2 + y2*y2
    return [y1 - y2, mpf(0), r2 - r1, r1*y2 - r2*y1]


def on(g, p):
    x, y = p
    return g[0]*(x*x + y*y) + g[1]*x + g[2]*y + g[3]


def _roots(a, b, c):
    """over C: the conjugate branch of one intersection sends later ones
       complex, and refusing them would hide exactly what we want to see"""
    if a == 0:
        return [] if b == 0 else [-c/b]
    s = sqrt(mpc(b*b - 4*a*c))
    return [(-b + s)/(2*a), (-b - s)/(2*a)]


def meet_line(g, L, hint, flip=False):
    A, B, C = L
    a, b, c, d = g
    if abs(B) > abs(A):
        p, q = -A/B, -C/B
        rs = _roots(a*(1 + p*p), b + c*p + 2*a*p*q, a*q*q + c*q + d)
        pts = [(x, p*x + q) for x in rs]
    else:
        p, q = -B/A, -C/A
        rs = _roots(a*(1 + p*p), c + b*p + 2*a*p*q, a*q*q + b*q + d)
        pts = [(p*y + q, y) for y in rs]
    if not pts:
        return None
    dist = lambda P: abs(P[0]-hint[0])**2 + abs(P[1]-hint[1])**2
    pts.sort(key=dist)
    if flip:
        return pts[1] if len(pts) > 1 else None
    return pts[0]


def meet(g1, g2, hint, flip=False):
    L = (g1[1]*g2[0] - g2[1]*g1[0], g1[2]*g2[0] - g2[2]*g1[0], g1[3]*g2[0] - g2[3]*g1[0])
    if max(abs(x) for x in L) == 0:
        return None
    return meet_line(g1, L, hint, flip)


def run(rb, vals, R, flipstep=None):
    P = [None]*69
    G = [None]*27
    T = [None]*len(rb['tips'])
    mir, smir = rb['mir'], rb['smir']

    def put(i, q):
        P[i] = q
        if mir[i] != i:
            P[mir[i]] = (-q[0], q[1])
    k = 0
    for s in rb['seeds']:
        if s['kind'] == 'pole':
            put(s['i'], (mpf(s['x']), mpf(s['y'])))
        elif s['kind'] == 'axis':
            put(s['i'], (mpf(0), vals[k])); k += 1
        elif s['kind'] == 'rim':
            from mpmath import cos, sin
            put(s['i'], (cos(vals[k]), sin(vals[k]))); k += 1
        else:
            put(s['i'], (vals[k], vals[k+1])); k += 2
    pt = lambda r: (T[r[0]] if isinstance(r, list) else P[r])
    ONE, Z = mpf(1), mpf(0)
    for j, st in enumerate(rb['steps']):
        op = st[0]
        fl = (j == flipstep)
        if op in ('s2', 's3'):
            g = (through_axis2(pt(st[2]), pt(st[3])) if op == 's2'
                 else through3(pt(st[2]), pt(st[3]), pt(st[4])))
            if g is None:
                return None
            G[st[1]] = g
            m = smir[st[1]]
            if G[m] is None:
                G[m] = g if m == st[1] else [g[0], -g[1], g[2], g[3]]
        else:
            if op == 'na':
                q = meet_line(G[st[2]], (ONE, Z, Z), (Z, mpf(st[3])), fl)
            elif op == 'nr':
                q = meet(G[st[2]], [ONE, Z, Z, -ONE], (mpf(st[3]), mpf(st[4])), fl)
            elif op == 'nx':
                q = meet(G[st[2]], G[st[3]], (mpf(st[4]), mpf(st[5])), fl)
            elif op == 't':
                q = meet(G[st[2]], [ONE, Z, Z, -R*R], (mpf(st[3]), mpf(st[4])), fl)
            else:
                q = (Z, R if st[2] > 0 else -R)
            if q is None:
                return None
            if op in ('t', 'ta'):
                T[st[1]] = q
            else:
                put(st[1], q)
    f = rb['fit']
    q = (T[f['tip']] if f['kind'] == 't' else
         pt(f['point']) if f['kind'] == 's' else P[f['node']])
    g = G[f['circle']]
    if g is None or q is None:
        return None
    return P, G, T, on(g, q)
