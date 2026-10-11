"""Write Table 4 of the paper (../../paper/tex/seeds.tex): the placed corners of the two constructions, the number
of real parameters each carries, and its value, the viewer's default rounded to four decimals (the corners at which
certify.py works)."""
import os, sys
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
from construct import REC, values
def corners(key):
    R=REC[key]; v=values(key,digits=4); k=0; out=[]
    for ni,kind in R['par']:
        nm=R['nodeNames'][ni]
        name='ring corner %s'%nm[1:] if nm[0]=='r' else '$P_{%s}$'%nm[1:]
        if kind==0: out.append((name,'free',2,r'$(%.4f,\ %.4f)$'%(v[k],v[k+1]))); k+=2
        elif kind==1: out.append((name,'axis',1,r'$y = %.4f$'%v[k])); k+=1
        else: out.append((name,'rim',1,r'$\theta = %.4f$'%v[k])); k+=1
    return out
A,B=corners('ARC_REC2'),corners('ARC_REC3')
assert len(A)==len(B)
lines=[r'\begin{table}[h]', r'\centering\footnotesize', r'\begin{tabular}{@{}llrl@{\qquad}llrl@{}}', r'\toprule',
       r'\multicolumn{4}{c}{pattern $22$} & \multicolumn{4}{c}{pattern $36$} \\',
       r'\cmidrule(r){1-4}\cmidrule(l){5-8}',
       r'corner & placed & & default & corner & placed & & default \\ \midrule']
for a,b in zip(A,B):
    lines.append('%s & %s & %d & %s & %s & %s & %d & %s \\\\'%(a+b))
lines+=[r'\midrule', r'total & & %d & & total & & %d & \\ \bottomrule'%(sum(x[2] for x in A),sum(x[2] for x in B)),
        r'\end{tabular}',
        r'\caption{The placed corners of the two constructions, the number of real parameters each carries, and its'
        r' value at the figures of Figure~\ref{fig:two}. Node numbers are those of the accompanying files; the ring'
        r' corners are numbered from $69$, and a corner on the rim is given by its polar angle.}',
        r'\label{tab:seeds}', r'\end{table}']
out=os.path.join(HERE,'..','..','paper','tex','seeds.tex')
if len(sys.argv)>1: out=sys.argv[1]
open(out,'w').write('\n'.join(lines)+'\n')
print('written',out)
