#!/usr/bin/env python3
"""Hunt for an arc Sri Yantra that can be drawn.

Where the question stands. One arranged coincidence is necessary - F = 0 is
UNSAT for both rings - and one is enough: F = 1 is realised, executable and
verified to 1e-15. The coincidence is never shallow. Seven square roots stand
between it and every free parameter, in both rings, in all four hundred orders
collected, and in ring 2 it is the same seven objects every time. Multiplying
the residual over all 2^7 branch signs gives the polynomial the parameter
satisfies over the rationals, and that polynomial has degree above 128.

None of which settles anything. A ruler-and-compass construction reaches
exactly the quadratic closure of the field it starts from, auxiliary circles
included - so allowing objects that are not part of the figure, which is
Chiodo's freedom, buys no number that was not already reachable. And all of
the above concerns a GENERIC choice of the other eleven parameters. We do not
need a family of arc yantras. We need one. Eleven free parameters against one
equation is a great deal of room, and a special choice may put the twelfth
number inside the closure, the way Chiodo's four heights are special.

So: give eleven parameters and the radius rational values, solve the
coincidence for the twelfth, and ask whether the answer is rational, or a
quadratic irrational, or a quartic one whose resolvent cubic has a rational
root. Any of the three and the entire figure is constructible - every seed is
then a constructible number and the sixty steps that follow are ruler and
compass throughout - and the arc yantra can be drawn.

Read this as a lottery, because that is what it is. A hit is a complete
answer. A million misses are not evidence of anything.

    python3 hunt.py --selftest                      # check the tests first
    python3 hunt.py --figure ring1 --trials 50000 -j 6
    python3 hunt.py --figure ring2 --trials 50000 -j 6 --seed 2

Everything it needs is in this file: both compiled recipes, the executor in
floating point for the search and in arbitrary precision for the verdict.
"""
import argparse, json, math, os, random, sys, time, warnings
warnings.filterwarnings('ignore')
from fractions import Fraction
from multiprocessing import Pool

from mpmath import mp, mpf, mpc, pslq, sqrt as mpsqrt

RECIPES = json.loads(r"""{"ring1":{"rb":{"R":1.2336067805243813,"seeds":[{"i":1,"kind":"free","x":0.325671002024,"y":-0.0486645285664,"a":null},{"i":0,"kind":"free","x":0.40255896897616716,"y":0.07995186209565727,"a":null},{"i":25,"kind":"free","x":0.455989849604,"y":-0.0310575207951,"a":null},{"i":3,"kind":"axis","x":0,"y":-0.4208528514168048,"a":null},{"i":27,"kind":"axis","x":0,"y":-0.6523609113718601,"a":null},{"i":55,"kind":"pole","x":0,"y":-1,"a":null},{"i":46,"kind":"axis","x":0,"y":0.458346729508,"a":null},{"i":44,"kind":"rim","x":-0.4061343130376837,"y":-0.9138133943937399,"a":-1.989016111660748},{"i":68,"kind":"pole","x":0,"y":1,"a":null},{"i":67,"kind":"axis","x":0,"y":0.908246941832509,"a":null},{"i":65,"kind":"rim","x":0.7968742534609248,"y":-0.6041452012315368,"a":-0.6486927410236439}],"tips":[{"orb":[0,7],"sides":[0,7],"h":[0.8463322499564097,0.8975006471470108]},{"orb":[0,7],"sides":[0,8],"h":[-0.8463322499563914,0.8975006471470108]},{"orb":[1,16],"sides":[1,17],"h":[0.3541252963917868,1.1816856449270774]},{"orb":[1,16],"sides":[2,16],"h":[-0.3541252963917869,1.1816856449270776]},{"orb":[3,15],"sides":[3,15],"h":[1.2333273077414075,-0.02625720700672418]},{"orb":[3,15],"sides":[3,15],"h":[-1.2333273077414073,-0.02625720700672418]},{"orb":[4,6],"sides":[4,6],"h":[1.11891390584475,0.5194394673751809]},{"orb":[4,6],"sides":[5,6],"h":[-1.11891390584475,0.5194394673751811]},{"orb":[4,22],"sides":[4,23],"h":[-0.3258721988646516,-1.1897869552858313]},{"orb":[4,22],"sides":[5,22],"h":[0.32587219886465113,-1.1897869552858313]},{"orb":[16,18],"sides":[17,18],"h":[-1.0989976818491842,-0.5603479135699916]},{"orb":[16,18],"sides":[16,18],"h":[1.0989976818491838,-0.560347913569992]},{"orb":[19,21],"sides":[20,21],"h":[-0.837481344905024,-0.9057652487790283]},{"orb":[19,21],"sides":[19,21],"h":[0.8374813449050453,-0.9057652487790274]}],"steps":[["s3",13,0,1,3],["s3",22,1,64,46],["t",8,23,-0.3258721988646516,-1.1897869552858313],["t",9,22,0.32587219886465113,-1.1897869552858313],["s2",18,27,65],["nx",45,18,23,-0.4341377901190528,-0.6380868897487062],["s3",4,59,55,[8]],["nr",52,5,-0.9860390076985804,0.16651448975028824],["t",10,18,-1.0989976818491842,-0.5603479135699916],["t",11,18,1.0989976818491838,-0.560347913569992],["s2",3,46,52],["t",4,3,1.2333273077414075,-0.02625720700672418],["t",5,3,-1.2333273077414073,-0.02625720700672418],["s2",15,3,[4]],["nx",48,4,15,0.7208237873390602,-0.2948437903487017],["s3",19,48,65,67],["nx",14,3,19,0.3946502086477745,0.4148407789195002],["nx",34,15,22,0.41566284424875966,-0.37979755316217484],["nr",47,15,0.9837805456326438,-0.17937624713082195],["s3",16,47,68,[11]],["t",2,17,0.3541252963917868,1.1816856449270774],["t",3,16,-0.3541252963917869,1.1816856449270776],["nx",57,5,17,-0.9041130350286282,0.0011055614554117826],["nx",32,3,16,0.725627448352,0.306685582764],["t",12,20,-0.837481344905024,-0.9057652487790283],["t",13,19,0.8374813449050453,-0.9057652487790274],["s2",21,44,[12]],["na",36,21,-0.9162877017129598],["s3",7,32,34,36],["nx",33,7,19,0.6127383292945495,0.0010658840105212328],["t",0,7,0.8463322499564097,0.8975006471470108],["t",1,8,-0.8463322499563914,0.8975006471470108],["s2",0,67,[0]],["nr",9,0,-0.4242460313063741,0.9055469645031069],["s3",1,12,14,[2]],["nx",7,2,23,-0.261031859642,0.0936879409509],["s2",12,0,7],["nx",56,5,21,-0.117540728686,-0.916080460435],["t",6,4,1.11891390584475,0.5194394673751809],["t",7,5,-1.11891390584475,0.5194394673751811],["nr",31,7,0.800617104794826,0.599176310871809],["s2",6,31,[6]],["na",61,6,0.6811008497583529],["s3",25,0,25,61],["na",16,1,-0.2196416970284696],["nx",20,7,25,0.522755567868,-0.189230721165],["s2",24,16,20],["nx",21,22,24,0.3776487192756541,-0.2037962070292666],["s3",10,21,25,27],["nx",24,10,19,0.5204110708748471,0.19634116281342834],["nx",15,1,25,0.3264601510301748,0.22064792181206297],["s2",9,15,24],["nx",53,3,25,0.18062161380855327,0.4493206204191347],["nx",2,13,24,0.202169610808,-0.215106292792],["nx",38,6,20,-0.20812835840107724,0.6756180254523059],["nx",10,0,17,-0.137596804234,0.907962937153],["nx",26,10,15,0.2454792833193286,-0.4066244272727091],["nx",29,9,23,-0.18065833798190625,0.23153306942326904],["nx",35,7,18,0.2315165978107828,-0.6483047101808462],["nx",13,1,6,0.43346517478669616,0.6572633684884314]],"fit":{"kind":"n","node":13,"circle":16},"mir":[6,5,4,3,2,1,0,8,7,12,11,10,9,19,18,17,16,15,14,13,23,22,21,20,28,62,49,27,24,30,29,37,43,42,41,40,36,31,39,38,35,34,33,32,64,59,46,51,50,26,48,47,54,63,52,55,60,58,57,45,56,61,25,53,44,66,65,67,68],"smir":[0,2,1,3,5,4,6,8,7,9,11,10,12,14,13,15,17,16,18,20,19,21,23,22,24,26,25]},"base":[0.325671002024,-0.0486645285664,0.40255896897616716,0.07995186209565727,0.455989849604,-0.0310575207951,-0.4208528514168048,-0.6523609113718601,0.458346729508,-1.989016111660748,0.908246941832509,-0.6486927410236439],"R":1.2336067805243813,"live":[0,1,6,7,8,9,10,11]},"ring2":{"rb":{"R":1.4533633622134141,"seeds":[{"i":2,"kind":"free","x":0.230582520975,"y":-0.246848061673,"a":null},{"i":1,"kind":"free","x":0.353808868089,"y":-0.0557960163148,"a":null},{"i":3,"kind":"axis","x":0,"y":-0.5185831526576876,"a":null},{"i":26,"kind":"free","x":0.2211871793376591,"y":-0.5024029446964909,"a":null},{"i":29,"kind":"free","x":-0.20892461747066643,"y":0.2601842882556762,"a":null},{"i":36,"kind":"axis","x":0,"y":-0.9461811945769771,"a":null},{"i":68,"kind":"pole","x":0,"y":1,"a":null},{"i":55,"kind":"pole","x":0,"y":-1,"a":null},{"i":65,"kind":"rim","x":0.7679166964008204,"y":-0.6405497228075665,"a":-0.6952139146379206}],"tips":[{"orb":[0,4],"sides":[0,4],"h":[1.1394740314647713,0.9021441094582932]},{"orb":[0,4],"sides":[0,5],"h":[-1.1394740314647713,0.9021441094582922]},{"orb":[1],"sides":[1,2],"h":[0,1.4533633622221158]},{"orb":[3,18],"sides":[3,18],"h":[1.3962235042277202,-0.4035157851282567]},{"orb":[3,18],"sides":[3,18],"h":[-1.3962235042277207,-0.4035157851282566]},{"orb":[4,19],"sides":[4,20],"h":[-0.7160635053761206,-1.2647205695037762]},{"orb":[4,19],"sides":[5,19],"h":[0.7160635053761206,-1.2647205695037764]},{"orb":[6,15],"sides":[6,15],"h":[1.3967051587413883,0.40184544566078995]},{"orb":[6,15],"sides":[6,15],"h":[-1.3967051587413886,0.40184544566078984]},{"orb":[7,16],"sides":[7,17],"h":[0.7304074160469536,1.2564911735528534]},{"orb":[7,16],"sides":[8,16],"h":[-0.7304074160469535,1.2564911735528534]},{"orb":[16,21],"sides":[17,21],"h":[-1.1294463613906023,-0.9146671412392486]},{"orb":[16,21],"sides":[16,21],"h":[1.1294463613906223,-0.9146671412392464]},{"orb":[22],"sides":[22,23],"h":[0,-1.4533633622478448]}],"steps":[["s2",15,3,26],["ta",13,-1.4533633622478448],["s3",22,1,30,[13]],["nx",34,15,22,0.4166104831331163,-0.4603723158004749],["nr",44,23,-0.3314609174080737,-0.9434689503269296],["s2",21,36,44],["t",11,21,-1.1294463613906023,-0.9146671412392486],["t",12,21,1.1294463613906223,-0.9146671412392464],["na",46,22,0.519671332585],["nr",47,15,0.9881425793728429,-0.1535390596245291],["s3",16,47,68,[12]],["t",9,17,0.7304074160469536,1.2564911735528534],["t",10,16,-0.7304074160469535,1.2564911735528534],["s3",7,34,36,[9]],["nx",32,7,16,0.735281286155,0.32962569248],["s2",3,32,46],["nr",52,3,-0.98809034732324,0.15387483721076342],["t",3,3,1.3962235042277202,-0.4035157851282567],["t",4,3,-1.3962235042277207,-0.4035157851282566],["s2",18,65,[3]],["na",27,18,-0.7376883062026446],["nx",45,18,23,-0.3945529876498614,-0.7123447152109137],["nx",35,7,18,0.21455848482598133,-0.7302154405463895],["s3",4,59,54,55],["nx",48,4,15,0.7292799414953426,-0.332196544202395],["t",5,4,-0.7160635053761206,-1.2647205695037762],["t",6,5,0.7160635053761206,-1.2647205695037764],["s3",19,48,65,[6]],["t",0,4,1.1394740314647713,0.9021441094582932],["t",1,5,-1.1394740314647713,0.9021441094582922],["na",67,19,0.934868675632344],["nx",57,5,17,-0.9230200127007564,0.0018280803653482058],["s2",0,67,[0]],["nr",31,7,0.7751404394471852,0.631788967245887],["nx",14,3,19,0.40733787195256077,0.46397588635834275],["nr",9,0,-0.363643250308201,0.9315382904128456],["t",7,15,1.3967051587413883,0.40184544566078995],["t",8,15,-1.3967051587413886,0.40184544566078984],["s2",6,31,[7]],["na",61,6,0.7287358625569401],["ta",2,1.4533633622221158],["s3",13,1,2,3],["nx",13,6,16,0.41092652587410056,0.7017942609143348],["s3",1,12,13,14],["nx",38,6,20,-0.21296574927504947,0.721522300938653],["na",16,1,-0.25503903573360764],["s2",24,2,16],["nx",21,22,24,0.39623122014127893,-0.23079216261274754],["s3",10,21,26,27],["nx",24,10,19,0.554422739102252,0.19999843291879077],["s2",9,24,29],["nx",15,1,9,0.35517635495839833,0.24153392520708566],["nx",20,7,24,0.557216988564,-0.206909476983],["nx",56,5,21,-0.0890146580208,-0.945985598466],["s3",25,15,20,61],["nx",7,2,23,-0.294341668776,0.100289409112],["nx",0,13,25,0.43271054831273176,0.09431479391980747],["nx",53,3,25,0.1815479476398729,0.5087717185075423],["nx",10,0,17,-0.108736428026,0.93457091974],["s2",12,0,7],["nx",25,10,25,0.487455868822,-0.0255055795869],["nx",33,7,19,0.6404430041172859,-0.008242875850457973]],"fit":{"kind":"s","circle":1,"point":[2]},"mir":[6,5,4,3,2,1,0,8,7,12,11,10,9,19,18,17,16,15,14,13,23,22,21,20,28,62,49,27,24,30,29,37,43,42,41,40,36,31,39,38,35,34,33,32,64,59,46,51,50,26,48,47,54,63,52,55,60,58,57,45,56,61,25,53,44,66,65,67,68],"smir":[0,2,1,3,5,4,6,8,7,9,11,10,12,14,13,15,17,16,18,20,19,21,23,22,24,26,25]},"base":[0.230582520975,-0.246848061673,0.353808868089,-0.0557960163148,-0.5185831526576876,0.2211871793376591,-0.5024029446964909,-0.20892461747066643,0.2601842882556762,-0.9461811945769771,-0.6952139146379206],"R":1.4533633622134141,"live":[2,3,4,5,6,7,8,9,10]}}""")


# ------------------------------------------------------- the executor, fast
# One pass of the compiled recipe in double precision. The search runs this a
# few hundred times per trial, so it is written out flat.

def unit(g):
    n = math.sqrt(sum(x * x for x in g))
    if n:
        g = [x / n for x in g]
    k = max(range(len(g)), key=lambda i: abs(g[i]))
    return g if g[k] > 0 else [-x for x in g]


def through3(p1, p2, p3):
    (x1, y1), (x2, y2), (x3, y3) = p1, p2, p3
    r1, r2, r3 = x1*x1 + y1*y1, x2*x2 + y2*y2, x3*x3 + y3*y3
    d3 = lambda a1,b1,c1,a2,b2,c2,a3,b3,c3: (a1*(b2*c3-b3*c2) - b1*(a2*c3-a3*c2)
                                             + c1*(a2*b3-a3*b2))
    return unit([ d3(x1,y1,1, x2,y2,1, x3,y3,1),
                 -d3(r1,y1,1, r2,y2,1, r3,y3,1),
                  d3(r1,x1,1, r2,x2,1, r3,x3,1),
                 -d3(r1,x1,y1, r2,x2,y2, r3,x3,y3)])


def through_axis2(p1, p2):
    (x1, y1), (x2, y2) = p1, p2
    r1, r2 = x1*x1 + y1*y1, x2*x2 + y2*y2
    return unit([y1 - y2, 0.0, r2 - r1, r1*y2 - r2*y1])


def _roots(a, b, c):
    if abs(a) < 1e-15:
        return [] if b == 0 else [-c / b]
    d = b * b - 4 * a * c
    if d < 0:
        return []
    s = math.sqrt(d)
    return [(-b + s) / (2 * a), (-b - s) / (2 * a)]


def meet_line(g, L, hint):
    A, B, C = L
    a, b, c, d = g
    # A = B = 0 with C non-zero is the line at infinity: the two circles are
    # concentric, there is no radical axis, and dividing by A here is how the
    # search first fell over. numpy hid it by returning inf with a warning.
    if A == 0 and B == 0:
        return None
    a, b, c, d = g
    if abs(B) > abs(A):
        p, q = -A / B, -C / B
        pts = [(x, p*x + q) for x in
               _roots(a*(1 + p*p), b + c*p + 2*a*p*q, a*q*q + c*q + d)]
    else:
        p, q = -B / A, -C / A
        pts = [(p*y + q, y) for y in
               _roots(a*(1 + p*p), c + b*p + 2*a*p*q, a*q*q + b*q + d)]
    if not pts:
        return None
    return min(pts, key=lambda P: (P[0]-hint[0])**2 + (P[1]-hint[1])**2)


def meet(g1, g2, hint):
    L = (g1[1]*g2[0] - g2[1]*g1[0], g1[2]*g2[0] - g2[2]*g1[0],
         g1[3]*g2[0] - g2[3]*g1[0])
    return None if max(abs(x) for x in L) < 1e-14 else meet_line(g1, L, hint)


def run(rb, vals, R):
    """returns (nodes, residual), or None if a step has no answer"""
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
            put(s['i'], (s['x'], s['y']))
        elif s['kind'] == 'axis':
            put(s['i'], (0.0, vals[k])); k += 1
        elif s['kind'] == 'rim':
            put(s['i'], (math.cos(vals[k]), math.sin(vals[k]))); k += 1
        else:
            put(s['i'], (vals[k], vals[k+1])); k += 2
    pt = lambda r: (T[r[0]] if isinstance(r, list) else P[r])
    for st in rb['steps']:
        op = st[0]
        if op in ('s2', 's3'):
            a = pt(st[2]); b = pt(st[3])
            if a is None or b is None:
                return None
            if op == 's2':
                g = through_axis2(a, b)
            else:
                c = pt(st[4])
                if c is None:
                    return None
                g = through3(a, b, c)
            G[st[1]] = g
            m = smir[st[1]]
            if G[m] is None:
                G[m] = g if m == st[1] else [g[0], -g[1], g[2], g[3]]
            continue
        if op == 'na':
            q = meet_line(G[st[2]], (1.0, 0.0, 0.0), (0.0, st[3]))
        elif op == 'nr':
            q = meet(G[st[2]], [1.0, 0.0, 0.0, -1.0], (st[3], st[4]))
        elif op == 'nx':
            q = meet(G[st[2]], G[st[3]], (st[4], st[5]))
        elif op == 't':
            q = meet(G[st[2]], [1.0, 0.0, 0.0, -R*R], (st[3], st[4]))
        else:
            q = (0.0, R if st[2] > 0 else -R)
        if q is None:
            return None
        if op in ('t', 'ta'):
            T[st[1]] = q
        else:
            put(st[1], q)
    f = rb['fit']
    q = pt(f['point']) if f['kind'] == 's' else P[f['node']]
    g = G[f['circle']]
    if g is None or q is None:
        return None
    x, y = q
    return P, g[0]*(x*x + y*y) + g[1]*x + g[2]*y + g[3]


def solve(rb, vals, R, idx, span=0.06):
    """the arranged coincidence, by a ladder of samples then bisection.

    Newton is not used: its steps leave the region where the sixty steps are
    all defined, and outside it there is no residual to follow.
    """
    def f(t):
        w = list(vals); w[idx] = t
        r = run(rb, w, R)
        return None if r is None else r[1]
    t0 = vals[idx]
    xs, ys, d = [], [], 1e-9
    while d <= span:
        for sgn in (1, -1):
            y = f(t0 + sgn*d)
            if y is not None:
                xs.append(t0 + sgn*d); ys.append(y)
        d *= 1.7
    y0 = f(t0)
    if y0 is not None:
        xs.append(t0); ys.append(y0)
    o = sorted(range(len(xs)), key=lambda i: xs[i])
    xs = [xs[i] for i in o]; ys = [ys[i] for i in o]
    best = None
    for i in range(1, len(xs)):
        if ys[i-1]*ys[i] < 0 and (best is None or abs(xs[i]-t0) < abs(best[1]-t0)):
            best = (xs[i-1], xs[i], ys[i-1], ys[i])
    if best is None:
        return None
    a, b, fa, fb = best
    for _ in range(120):
        m = (a + b)/2
        fm = f(m)
        if fm is None:
            return None
        if fa*fm <= 0:
            b, fb = m, fm
        else:
            a, fa = m, fm
    return (a + b)/2


# --------------------------------------------- the executor, arbitrary digits

def mp_run(rb, vals, R):
    P = [None]*69
    G = [None]*27
    T = [None]*len(rb['tips'])
    mir, smir = rb['mir'], rb['smir']
    def put(i, q):
        P[i] = q
        if mir[i] != i:
            P[mir[i]] = (-q[0], q[1])
    def th3(p1, p2, p3):
        (x1,y1),(x2,y2),(x3,y3) = p1,p2,p3
        r1,r2,r3 = x1*x1+y1*y1, x2*x2+y2*y2, x3*x3+y3*y3
        d3 = lambda a1,b1,c1,a2,b2,c2,a3,b3,c3: (a1*(b2*c3-b3*c2)-b1*(a2*c3-a3*c2)
                                                 + c1*(a2*b3-a3*b2))
        one = mpf(1)
        return [ d3(x1,y1,one,x2,y2,one,x3,y3,one), -d3(r1,y1,one,r2,y2,one,r3,y3,one),
                 d3(r1,x1,one,r2,x2,one,r3,x3,one), -d3(r1,x1,y1,r2,x2,y2,r3,x3,y3)]
    def th2(p1, p2):
        (x1,y1),(x2,y2) = p1,p2
        r1,r2 = x1*x1+y1*y1, x2*x2+y2*y2
        return [y1-y2, mpf(0), r2-r1, r1*y2-r2*y1]
    def rts(a, b, c):
        # A negative discriminant means the circles miss: no point, as in the
        # floating-point executor. Taking mpmath's square root of it instead
        # returns an mpc, and two complex points cannot be sorted by distance
        # to the hint - which is how the overnight run died.
        if a == 0:
            return [] if b == 0 else [-c/b]
        d = b*b - 4*a*c
        if d < 0:
            return []
        s = mpsqrt(d)
        return [(-b+s)/(2*a), (-b-s)/(2*a)]
    def mline(g, L, hint):
        A, B, C = L
        if A == 0 and B == 0:
            return None
        a, b, c, d = g
        if abs(B) > abs(A):
            p, q = -A/B, -C/B
            pts = [(x, p*x+q) for x in rts(a*(1+p*p), b+c*p+2*a*p*q, a*q*q+c*q+d)]
        else:
            p, q = -B/A, -C/A
            pts = [(p*y+q, y) for y in rts(a*(1+p*p), c+b*p+2*a*p*q, a*q*q+b*q+d)]
        if not pts:
            return None
        return min(pts, key=lambda P_: (P_[0]-hint[0])**2 + (P_[1]-hint[1])**2)
    def mmeet(g1, g2, hint):
        L = (g1[1]*g2[0]-g2[1]*g1[0], g1[2]*g2[0]-g2[2]*g1[0], g1[3]*g2[0]-g2[3]*g1[0])
        return None if max(abs(x) for x in L) == 0 else mline(g1, L, hint)
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
    for st in rb['steps']:
        op = st[0]
        if op in ('s2', 's3'):
            g = th2(pt(st[2]), pt(st[3])) if op == 's2' else th3(pt(st[2]), pt(st[3]), pt(st[4]))
            G[st[1]] = g
            m = smir[st[1]]
            if G[m] is None:
                G[m] = g if m == st[1] else [g[0], -g[1], g[2], g[3]]
            continue
        if op == 'na':
            q = mline(G[st[2]], (ONE, Z, Z), (Z, mpf(st[3])))
        elif op == 'nr':
            q = mmeet(G[st[2]], [ONE, Z, Z, -ONE], (mpf(st[3]), mpf(st[4])))
        elif op == 'nx':
            q = mmeet(G[st[2]], G[st[3]], (mpf(st[4]), mpf(st[5])))
        elif op == 't':
            q = mmeet(G[st[2]], [ONE, Z, Z, -R*R], (mpf(st[3]), mpf(st[4])))
        else:
            q = (Z, R if st[2] > 0 else -R)
        if q is None:
            return None
        if op in ('t', 'ta'):
            T[st[1]] = q
        else:
            put(st[1], q)
    f = rb['fit']
    q = pt(f['point']) if f['kind'] == 's' else P[f['node']]
    g = G[f['circle']]
    if g is None or q is None:
        return None
    x, y = q
    return g[0]*(x*x + y*y) + g[1]*x + g[2]*y + g[3]


def refine(rb, vals, R, slot, t0, dps):
    """the root again, to dps digits, from the double-precision bracket"""
    old = mp.dps
    mp.dps = dps + 15
    try:
        def f(t):
            v = [mpf(x) for x in vals]
            v[slot] = t
            return mp_run(rb, v, mpf(R))
        a, b = mpf(t0) - mpf(10)**-9, mpf(t0) + mpf(10)**-9
        fa, fb = f(a), f(b)
        if fa is None or fb is None or fa*fb > 0:
            return None
        for _ in range(int(dps*3.4) + 40):
            m = (a + b)/2
            fm = f(m)
            if fm is None:
                return None
            if fa*fm <= 0:
                b, fb = m, fm
            else:
                a, fa = m, fm
        return (a + b)/2
    finally:
        mp.dps = old


# ------------------------------------------------------- the number theory
# Every test asks whether t satisfies a small integer polynomial to far better
# than the accident rate. The nearest rational with denominator up to Q sits
# about 1/Q^2 away, so agreement to 1e-28 with Q = 1e9 is an accident about
# once in 1e10 trials; the quadratic and quartic tests have similar margins.

def _rel(t, deg, maxcoeff, tol):
    v = [mpf(1)]
    for _ in range(deg):
        v.append(v[-1]*t)
    r = pslq(v, tol=tol, maxcoeff=maxcoeff, maxsteps=2000)
    if r is None or all(c == 0 for c in r):
        return None
    err = abs(sum(mpf(c)*x for c, x in zip(r, v)))
    scale = max(abs(mpf(c)) for c in r)
    return r if err < tol*scale*100 else None


def resolvent_has_rational_root(c):
    """c = [c0..c4] for c4 t^4 + ... + c0. An irreducible quartic is
    constructible exactly when its resolvent cubic has a rational root: with
    one, the splitting field is a tower of quadratics and nothing else."""
    from math import gcd
    F = Fraction
    a, b, cc, d, e = (F(c[4]), F(c[3]), F(c[2]), F(c[1]), F(c[0]))
    if a == 0:
        return False
    b, cc, d, e = b/a, cc/a, d/a, e/a
    p2, p1, p0 = -cc, b*d - 4*e, -(b*b*e - 4*cc*e + d*d)
    den = 1
    for x in (p2, p1, p0):
        den = den*x.denominator // gcd(den, x.denominator)
    A, B, C, D = den, int(p2*den), int(p1*den), int(p0*den)
    def divisors(n):
        n = abs(n)
        if n == 0:
            return [0]
        out, i = [], 1
        while i*i <= n:
            if n % i == 0:
                out += [i, n//i]
            i += 1
        return sorted(set(out))
    for p in divisors(D):
        for q in divisors(A):
            for s in (1, -1):
                y = F(s*p, q)
                if A*y**3 + B*y**2 + C*y + D == 0:
                    return True
    return False


def classify(t, tol):
    r = _rel(t, 1, 10**9, tol)
    if r:
        return ('rational', r)
    r = _rel(t, 2, 10**7, tol)
    if r:
        return ('quadratic', r)
    r = _rel(t, 4, 10**5, tol)
    if r and resolvent_has_rational_root([int(x) for x in r]):
        return ('quartic, 2-group', r)
    return None


# ------------------------------------------------------------- the search

def valid(rb, vals, R, slot):
    t = solve(rb, vals, R, slot)
    if t is None:
        return None
    v = list(vals); v[slot] = t
    c = run(rb, v, R)
    if c is None or abs(c[1]) > 1e-11:
        return None
    if any(p is None or abs(p[0]) > 50 or abs(p[1]) > 50 for p in c[0]):
        return None
    return t


def walk(rb, base, R0, slot, rng, qmax):
    """Give the parameters rational values one at a time, never leaving the
    region where the figure exists.

    Rounding them all at once does not work. The region is a thousandth wide
    in some directions and unbounded in others, so any common denominator
    small enough to be interesting throws the figure away on the first narrow
    parameter - blind rounding produced a figure twice in six hundred tries.
    Taken one at a time, each parameter gets the smallest denominator that
    still leaves a figure standing, and the others shift to accommodate it
    when the coincidence is re-solved. Every trial then ends on a genuinely
    rational point of the variety.
    """
    vals, R = list(base), R0
    idx = [i for i in range(len(base)) if i != slot] + ['R']
    rng.shuffle(idx)
    dens = {}
    for i in idx:
        x = R if i == 'R' else vals[i]
        got, q = None, 1
        while q <= qmax and got is None:
            c0 = Fraction(x).limit_denominator(q)
            cands = [c0, c0 + Fraction(1, q), c0 - Fraction(1, q)]
            rng.shuffle(cands)
            for c in cands:
                v, RR = list(vals), R
                if i == 'R':
                    RR = float(c)
                else:
                    v[i] = float(c)
                if valid(rb, v, RR, slot) is not None:
                    got, dens[i] = (v, RR), c.denominator
                    break
            q *= 2
        if got is None:
            return None
        vals, R = got
    return vals, R, dens


def trial(args):
    figure, seed, n, qmax = args
    rec = RECIPES[figure]
    rb, base, R0, live = rec['rb'], rec['base'], rec['R'], rec['live']
    rng = random.Random(seed)
    hits, done, solved, worst = [], 0, 0, []
    bad = norefine = 0
    for _ in range(n):
        done += 1
        try:
            r = one_trial(rb, base, R0, live, rng, qmax)
        except Exception:
            bad += 1
            continue
        if r is None:
            continue
        kind, payload = r
        if kind == 'norefine':
            solved += 1
            worst.append(payload)
            norefine += 1
            continue
        solved += 1
        worst.append(payload[0])
        if payload[1] is not None:
            hits.append(payload[1])
    return done, solved, hits, worst, bad, norefine


def one_trial(rb, base, R0, live, rng, qmax):
    for _ in (0,):
        slot = rng.choice(live)
        w = walk(rb, base, R0, slot, rng, qmax)
        if w is None:
            return None
        vals, R, dens = w
        t = valid(rb, vals, R, slot)
        if t is None:
            return None
        worst = max(dens.values())
        # Refine every root to forty digits before testing anything. Gating on
        # the double-precision root would be cheaper and wrong: it carries
        # about twelve honest digits, the same order as the distance from a
        # random real to the nearest small quadratic, so the gate would throw
        # away exactly the hits it exists to catch. Forty digits leaves the
        # tests twenty-odd orders of margin and costs a twentieth of the walk.
        tt = refine(rb, vals, R, slot, t, 40)
        if tt is None:
            return ('norefine', worst)
        mp.dps = 40
        c = classify(tt, mpf(10)**-28)
        if not c:
            return ('ok', (worst, None))
        tt2 = refine(rb, vals, R, slot, t, 150)
        if tt2 is not None:
            mp.dps = 150
            c = classify(tt2, mpf(10)**-100) or c
            tt = tt2
        return ('ok', (worst, dict(kind=c[0], rel=[int(x) for x in c[1]], slot=slot,
                                   dens={str(k): v for k, v in dens.items()},
                                   vals=vals, R=R, t=mp.nstr(tt, 100))))


def selftest():
    mp.dps = 60
    cases = [('22/7', mpf(22)/7, 'rational'),
             ('sqrt 2', mpsqrt(mpf(2)), 'quadratic'),
             ('(1+sqrt 5)/2', (1 + mpsqrt(mpf(5)))/2, 'quadratic'),
             ('sqrt2 + sqrt3', mpsqrt(mpf(2)) + mpsqrt(mpf(3)), 'quartic, 2-group'),
             ('cbrt 2', mpf(2)**(mpf(1)/3), None),
             ('2cos(2pi/7)', 2*mp.cos(2*mp.pi/7), None),
             ('pi', mp.pi, None)]
    ok = True
    for name, x, want in cases:
        got = classify(x, mpf(10)**-45)
        g = got[0] if got else None
        ok &= (g == want)
        print('   %-16s -> %-18s expected %-18s %s'
              % (name, g, want, 'OK' if g == want else 'WRONG'))
    # and the executor, at both precisions, on a known figure
    for fig in sorted(RECIPES):
        rec = RECIPES[fig]
        rb, base, R = rec['rb'], rec['base'], rec['R']
        slot = rec['live'][0]
        t = valid(rb, base, R, slot)
        if t is None:
            print('   %-8s the base figure does not solve   WRONG' % fig); ok = False
            continue
        tt = refine(rb, base, R, slot, t, 60)
        mp.dps = 75
        v = [mpf(x) for x in base]; v[slot] = tt
        res = abs(mp_run(rb, v, mpf(R)))
        good = res < mpf(10)**-60
        ok &= good
        print('   %-8s base figure solves; residual at 60 digits %s   %s'
              % (fig, mp.nstr(res, 4), 'OK' if good else 'WRONG'))
    print('selftest %s' % ('passed' if ok else 'FAILED'))
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--figure', default='ring1', choices=sorted(RECIPES) or None)
    ap.add_argument('--trials', type=int, default=5000)
    ap.add_argument('--qmax', type=int, default=8192)
    ap.add_argument('-j', type=int, default=6)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--selftest', action='store_true')
    a = ap.parse_args()
    if a.selftest:
        sys.exit(0 if selftest() else 1)
    rec = RECIPES[a.figure]
    print('%s: %d parameters the residual sees, denominators up to %d, %d workers'
          % (a.figure, len(rec['live']), a.qmax, a.j), flush=True)
    # Small chunks, not one block per worker. A block of thirty thousand
    # trials prints nothing for eight hours and loses everything if the run is
    # interrupted; chunks report as they land and hits are written at once.
    CH = 25
    nch = max(1, a.trials // CH)
    work = [(a.figure, a.seed*1000000 + i, CH, a.qmax) for i in range(nch)]
    out = 'hits_%s.json' % a.figure
    t0 = time.time()
    D = V = BAD = NR = 0
    HITS, W = [], []
    with Pool(a.j) as pool:
        for done, solved, hits, worst, bad, nr in pool.imap_unordered(trial, work):
            D += done; V += solved; W += worst; BAD += bad; NR += nr
            if hits:
                HITS += hits
                json.dump(HITS, open(out, 'w'), indent=1)
                for h in hits:
                    print('  *** HIT: %s, slot %d -- written to %s'
                          % (h['kind'], h['slot'], out), flush=True)
            if True:
                el = time.time() - t0
                print('  %d tried, %d rational figures, %d hits, %d skipped, '
                      '%.0f s, %.1f/s' % (D, V, len(HITS), BAD + NR, el, D/max(el, 1)),
                      flush=True)
    print('\n%d trials, %d fully rational figures, %.0f s'
          % (D, V, time.time() - t0))
    if BAD or NR:
        print('%d trials threw and %d could not be refined; both are skipped'
              % (BAD, NR))
    if W:
        W.sort()
        print('largest denominator a figure needed: best %d, median %d, worst %d'
              % (W[0], W[len(W)//2], W[-1]))
    if not HITS:
        print('\nno rational, quadratic or constructible-quartic root.')
        print('That is the expected outcome and it proves nothing. Raise --trials,')
        print('change --seed, or try the other ring.')
        return
    print('\n%d HIT(S) -- written to %s' % (len(HITS), out))
    for h in HITS:
        print('   %s   slot %d   denominators %s'
              % (h['kind'], h['slot'], sorted(h['dens'].values())))
        print('      relation %s' % h['rel'])
        print('      t = %s' % h['t'][:60])
    print('\nSend me that file. If it holds up, the arc yantra is constructible')
    print('and the construction follows from it.')


if __name__ == '__main__':
    main()
