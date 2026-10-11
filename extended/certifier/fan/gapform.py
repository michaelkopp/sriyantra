"""Exact factorization of the conditions in the fan chart (zD = 0, zU = inf, y10 = 1, p_k = t_k t_{k+1} ... t_10):
   p_b - p_a = (t_b ... t_10) * g(a,b),   g(a,b) = 1 - t_a ... t_{b-1}   (a >= 1);   p_b - p_0 = t_b ... t_10.
Each term of a condition is a monomial in the t_k times a product of the g(a,b)."""
import json
from collections import Counter
D=json.load(open('tri.json')); PAIRL=[tuple(p) for p in D['pairs']]
def factor(cnt):
    E=Counter(); G=Counter()
    for i,e in cnt.items():
        a,b=PAIRL[int(i)]
        if a==12 or b==12: continue
        for j in range(b,11): E[j]+=e
        if a>0: G[(a,b)]+=e
    return E,G
def diff(x,y,s=1):
    r=Counter(x); 
    for k,v in y.items(): r[k]-=s*v
    return Counter({k:v for k,v in r.items() if v})
def logratio(name):
    """(E, G) of log(term1/term2) for a binomial; for DOWN: of log(T1/T0) and the ratio T2/T1"""
    T=D['eqterms'][name]; F=[factor(c) for s,c in T]
    if len(T)==2: return diff(F[0][0],F[1][0]),diff(F[0][1],F[1][1]),None
    # DOWN: T0 - T1 - T2 = 0 (signs +,-,-), L = log((T1+T2)/T0) = log(T1/T0) + log(1 + T2/T1)
    return diff(F[1][0],F[0][0]),diff(F[1][1],F[0][1]),(diff(F[2][0],F[1][0]),diff(F[2][1],F[1][1]))
if __name__=='__main__':
    tot_E=Counter(); tot_G=Counter()
    for name,s in [('target',1),('Q21',-1),('Q59',1),('DOWN',1)]:
        E,G,extra=logratio(name)
        print(name,'  t-exponents',dict(E),'\n      g-exponents',dict(G))
        if extra: print('      T2/T1 = t^',dict(extra[0]),' g^',dict(extra[1]))
        for k,v in E.items(): tot_E[k]+=s*v
        for k,v in G.items(): tot_G[k]+=s*v
    print('Phi: t-exponents',{k:v for k,v in tot_E.items() if v})
    print('Phi: g-exponents',{k:v for k,v in sorted(tot_G.items()) if v})
    for name in ['Q43','apexU','apexD','UP']:
        E,G,extra=logratio(name); print(name,dict(E),dict(G),extra)
