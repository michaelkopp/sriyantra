import json, os
D=json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','certifier','cones','tri.json'))); P=[tuple(p) for p in D['pairs']]
def nm(i):
    return 'D' if i==0 else ('U' if i==12 else str(i-1))
def br(a,b):
    x,y=nm(a),nm(b)
    if x=='D' or y=='U': return r'\br{%s}{%s}'%(x,y)
    return r'\br{%s}{%s}'%(x,y)
def mono(v):
    num=[];den=[]
    for (a,b),e in zip(P,v):
        if e==0: continue
        t=br(a,b)+('^{%d}'%abs(e) if abs(e)>1 else '')
        (num if e>0 else den).append(t)
    return num,den
def frac(v):
    num,den=mono(v)
    n=' '.join(num) if num else '1'; d=' '.join(den) if den else '1'
    return r'\dfrac{%s}{%s}'%(n,d)
def eq(v):  # monomial = 1 written as num = den
    num,den=mono(v)
    return '%s = %s'%(' '.join(num) or '1',' '.join(den) or '1')
def ge(v):
    num,den=mono(v)
    return '%s \\ge %s'%(' '.join(num) or '1',' '.join(den) or '1')
M=D['mono']
out=[]
out.append(r'\begin{longtable}{@{}l>{$}l<{$}@{}}')
out.append(r'\caption{The eight equations in the thirteen points, $\br{a}{b} = |p_b - p_a|$, with $a, b$ the levels $0, \dots, 10$ and $D = z_D$, $U = z_U$.}\label{tab:eqs13}\\')
out.append(r'\toprule equation & \text{form} \\ \midrule \endfirsthead \toprule equation & \text{form} \\ \midrule \endhead')
names={'Q21':r'$E_2$','Q43':r'$E_3$','Q59':r'$E_4$','target':r'$E_1$','apexU':r'$A_U$','apexD':r'$A_D$'}
for n in ['target','Q21','Q43','Q59','apexU','apexD']:
    out.append(r'%s & %s \\'%(names[n],eq(M[n])))
for n,lab in (('UP',r'$C_U$'),('DOWN',r'$C_D$')):
    out.append(r'%s & 1 = %s + %s \\'%(lab,frac(M[n+'.m1']),frac(M[n+'.m2'])))
out.append(r'\bottomrule\end{longtable}')
out.append('')
out.append(r'\begin{longtable}{@{}l>{$}l<{$}@{}}')
out.append(r'\caption{The twenty-five inequalities that define $\mathcal{L}$. A non-crossing clause $\mathrm{nc}\,AB$ reads $m \ge 1$; a corner clause reads $m_1 + m_2 \le 1$.}\label{tab:clauses13}\\')
out.append(r'\toprule clause & \text{form} \\ \midrule \endfirsthead \toprule clause & \text{form} \\ \midrule \endhead')
for n,c in D['clauses'].items():
    import re
    t=re.sub(r'([DU])(\d)',r'\1_\2',n.split(' ')[1])
    lab=('non-crossing $%s$' if n.startswith('nc') else 'corner $%s$')%t
    if c['type']=='bin': out.append(r'%s & %s \\'%(lab,ge(c['m'][0])))
    else: out.append(r'%s & %s + %s \le 1 \\'%(lab,frac(c['m'][0]),frac(c['m'][1])))
out.append(r'\bottomrule\end{longtable}')
open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','paper','tex','tables13.tex'),'w').write('\n'.join(out))
print('\n'.join(out[:14]))
