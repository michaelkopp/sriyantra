"""The classical straight yantra built ring by ring inside nested kites. Usage: python3 drawkites.py [point.json [out_base]]; defaults classical_pt.json, ../paper/figs/kites"""
import json,sys,os,numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,os.path.join(HERE,'fig1'))
from ringpaint import classify9, kites9
from model import base
from palette import RING,PALE
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
d=json.load(open(sys.argv[1] if len(sys.argv)>1 else os.path.join(HERE,'classical_pt.json')))
out=sys.argv[2] if len(sys.argv)>2 else os.path.join(HERE,'..','paper','figs','kites')
y=d['y']; W=d['W']; kap=d.get('kappa',0.6)
pts,F,ring,pale,outer=classify9(y,W,kap)
KS=kites9(y,W,kap)
plt.rcParams.update({'font.family':'serif','font.serif':['DejaVu Serif'],'font.size':8,'mathtext.fontset':'cm','pdf.fonttype':42})
INK='#2a2a2a'; KITE='#111111'
NAMES=[('U_3','D_5'),('U_4','D_4'),('U_2','D_2'),('U_1','D_1'),('U_0','D_0')]
RINGNAME=['the central triangle','the ring of eight','the inner ring of ten','the outer ring of ten','the ring of fourteen']
# outer vertices of ring k that lie on the legs of kite k+1: the vertices of ring cells not shared with deeper cells
def outer_vertices(k):
    vs=set()
    for i,r in ring.items():
        if r!=k: continue
        for v in F[i]['verts']: vs.add(v)
    inner=set()
    for i in list(ring)+list(pale):
        if (i in ring and ring[i]<k) or (i in pale and pale[i]<k):
            for v in F[i]['verts']: inner.add(v)
    return [pts[v] for v in vs-inner]
def onkite(p,K,tol=1e-6):
    for a,b in ((K[0],K[1]),(K[1],K[2]),(K[0],K[3]),(K[3],K[2])):
        d=b-a; t=np.dot(p-a,d)/np.dot(d,d)
        if -1e-9<=t<=1+1e-9 and np.linalg.norm(a+t*d-p)<tol: return True
    return False
ext_=[]
for K in KS:
    K=np.array(K); xm=K[1][0]*1.06; h=(K[0][1]-K[2][1])*1.08; ext_.append(2*xm/h)
fig,axs=plt.subplots(1,5,figsize=(6.4,2.5),gridspec_kw={'wspace':0.08,'width_ratios':ext_})
for k,ax in enumerate(axs):
    for i,r in ring.items():
        if r<=k: ax.add_patch(Polygon(pts[F[i]['verts']],fc=RING[r],ec=INK,lw=0.35,joinstyle='miter'))
    for i,p in pale.items():
        if p<k: ax.add_patch(Polygon(pts[F[i]['verts']],fc=PALE[p],ec=INK,lw=0.35,joinstyle='miter'))
    K=np.array(KS[k])
    # a side is drawn solid when it carries every outer vertex of the ring on its half, dashed otherwise
    ov=[p for p in outer_vertices(k) if abs(p[0])>1e-9]
    up_ok=all(onkite(p,K) for p in ov if p[1]>K[1][1]); dn_ok=all(onkite(p,K) for p in ov if p[1]<K[1][1])
    dash=(0,(3.5,1.8))
    ax.plot([K[3][0],K[0][0],K[1][0]],[K[3][1],K[0][1],K[1][1]],color=KITE,lw=1.0,ls='-' if up_ok else dash,solid_joinstyle='miter')
    ax.plot([K[3][0],K[2][0],K[1][0]],[K[3][1],K[2][1],K[1][1]],color=KITE,lw=1.0,ls='-' if dn_ok else dash,solid_joinstyle='miter')
    # the outer vertices of the newest ring: on the kite (filled) or off it (open)
    for p in outer_vertices(k):
        if abs(p[0])<1e-9: continue
        on=onkite(p,K)
        ax.plot([p[0]],[p[1]],'o',ms=2.2,mfc=KITE if on else 'white',mec=KITE,mew=0.6,zorder=6)
    xm=K[1][0]*1.06; ym0=K[2][1]; ym1=K[0][1]; pad=0.04*(ym1-ym0)
    ax.set_xlim(-xm,xm); ax.set_ylim(ym0-pad,ym1+pad); ax.set_aspect('equal'); ax.axis('off')
    ax.text(0.5,-0.03,'(%s)  $%s$, $%s$'%('abcde'[k],*NAMES[k]),transform=ax.transAxes,ha='center',va='top',fontsize=7.5)
plt.savefig(out+'.pdf',bbox_inches='tight',pad_inches=0.02); plt.savefig(out+'.png',dpi=220,bbox_inches='tight',pad_inches=0.02)
