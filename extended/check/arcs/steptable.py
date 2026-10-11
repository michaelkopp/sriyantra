"""Write the step tables of the two constructions for the appendix (../../paper/tex/steps_arcs.tex): every step with the
object it draws, the objects it is drawn from, and for a node its position at the corners of Table 4."""
import os, sys, json, re
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
from construct import REC, build, values
def side_tex(name,m):
    t,kind=name.split('.'); T='%s_{%s}'%(t[0],t[1:])
    if kind=='base': return r'\beta(%s)'%T
    return (r"\lambda'(%s)" if m else r'\lambda(%s)')%T
def node_tex(R,i,m=0):
    n=re.sub(r'\D','',R['nodeNames'][i]); return "P_{%s}%s"%(n,"'" if m else '')
def rows(key):
    R=REC[key]; N,_=build(key,values(key,digits=4)); kinds={ni:k for ni,k in R['par']}
    out=[]
    for k,op in enumerate(R['ops'],1):
        kind,i=op[0],op[1]
        pos=lambda j: r'$(%.3f, %.3f)$'%(N[j][0]+0.0,N[j][1]+0.0)
        if kind=='seed':
            how={0:'placed',1:'placed on the axis',2:'placed on the rim'}[kinds[i]]
            out.append((k,'$%s$'%node_tex(R,i),how,pos(i)))
        elif kind=='pole':
            out.append((k,'$%s$'%node_tex(R,i),'pole',pos(i)))
        elif kind in ('s2','s3'):
            pts=', '.join(node_tex(R,op[a],op[a+1]) for a in range(2,len(op),2))
            out.append((k,'$%s$'%side_tex(R['sideNames'][i],0),'$%s$'%pts,''))
        else:
            a=side_tex(R['sideNames'][op[2]],op[3])
            b={'nx':lambda: side_tex(R['sideNames'][op[4]],op[5]),'na':lambda:r'\text{axis}','nr':lambda:r'\text{rim}'}[kind]()
            out.append((k,'$%s$'%node_tex(R,i),r'$%s \cap %s$'%(a,b),pos(i)))
    return out
def table(key,pat,label):
    rs=rows(key); h=(len(rs)+1)//2; L=rs[:h]; Rr=rs[h:]+[('','','','')]*(2*h-len(rs))
    s=[r'\begin{table}[p]',r'\centering\scriptsize',r'\setlength{\tabcolsep}{3.2pt}',
       r'\begin{tabular}{@{}rlll@{\qquad}rlll@{}}',r'\toprule',
       r'& drawn & from & position & & drawn & from & position \\ \midrule']
    for a,b in zip(L,Rr):
        s.append(' & '.join(str(x) for x in a)+' & '+' & '.join(str(x) for x in b)+r' \\')
    s+= [r'\bottomrule',r'\end{tabular}',
         r'\caption{The construction of pattern $%s$, step by step. $\beta(T)$ is the base of $T$, $\lambda(T)$ its right leg and $\lambda^\prime(T)$ its left leg, and $P^\prime$ is the mirror image of $P$. A curve is drawn through the nodes listed, and a node as the meeting listed; of the two points of a meeting it is the one at the position given, which is its position at the corners of Table~\ref{tab:seeds} to three decimals.}'%pat,
         r'\label{%s}'%label,r'\end{table}','']
    return '\n'.join(s)
tex=table('ARC_REC2','22','tab:steps22')+'\n'+table('ARC_REC3','36','tab:steps36')
open(os.path.join(HERE,'..','..','paper','tex','steps_arcs.tex'),'w').write(tex)
print('written', len(tex))
