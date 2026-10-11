"""Figure 2: the straight extended figure in which one concurrence fails, its rings shaded and its five kites drawn;
(a) as constructed, with the levels and an enlargement of the miss; (b), (c) the same figure with the triangles in
the positions of the wall patterns 22 and 36 carried out beyond the fifth kite.
Usage: python3 drawstraight.py [point.json [out_base]]; defaults straight_pt.json, ../paper/figs/straight_nearmiss"""
import json,sys,os,numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,os.path.join(HERE,'fig1'))
from ringpaint import cells, kite
from walls import build, segs
from model import base
from palette import RING,PALE
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
d=json.load(open(sys.argv[1] if len(sys.argv)>1 else os.path.join(HERE,'straight_pt.json')))
out=sys.argv[2] if len(sys.argv)>2 else os.path.join(HERE,'..','paper','figs','straight_nearmiss')
kap=d.get('kappa',0.6); q=d['q']; y=q[1:12]; zD,zU=q[0],q[12]; W0=d['W']; W={k:kap*v for k,v in W0.items()}
plt.rcParams.update({'font.family':'serif','font.serif':['DejaVu Serif'],'font.size':8,'mathtext.fontset':'cm','pdf.fonttype':42})
INK='#2a2a2a'; KITE='#111111'
lev=y[4]; wD2=W['D2']*(lev-y[1])/(y[8]-y[1]); g=wD2-W['U4']
def kites(T):
    C=lambda t: (T[t][2],T[t][1])
    K=[kite(T[u][0],(W[u],y[base[u]]),T[dd][0],(W[dd],y[base[dd]])) for u,dd in (('U3','D5'),('U4','D4'),('U2','D2'),('U1','D1'))]
    # the fifth kite: the legs of U0 and D0 from their apexes to the points where they cross
    sU=T['U0'][2]/(T['U0'][1]-T['U0'][0]); sD=T['D0'][2]/(T['D0'][1]-T['D0'][0])
    K.append(kite(T['U0'][0],(sU*(0-T['U0'][0]),0.0),T['D0'][0],(sD*(0-T['D0'][0]),0.0)))
    return K
def paint(ax,T,lw,klw):
    S,_=segs(T); pts,F,depth,outer,nb=cells(S)
    ring={i:(11-depth[i])//2 for i in range(len(F)) if i!=outer and depth[i]%2==1}
    for i in range(len(F)):
        if i==outer: continue
        if i in ring: c=RING[ring[i]]
        else: c=PALE[min(ring[j] for j in nb[i] if j in ring)]
        ax.add_patch(Polygon(pts[F[i]['verts']],fc=c,ec='none'))
    for (a,b) in S: ax.plot([a[0],b[0]],[a[1],b[1]],color=INK,lw=lw,solid_capstyle='butt')
    for K in kites(T):
        K=np.array(K); ax.plot(K[[3,0,1,2,3],0],K[[3,0,1,2,3],1],color=KITE,lw=klw,solid_joinstyle='miter',zorder=4)
    ax.plot([W['U4']],[lev],'o',ms=2.2,mfc='white',mec=KITE,mew=0.6,zorder=6)
    xs=[abs(p[0]) for s in S for p in s]; ys=[p[1] for s in S for p in s]
    return max(xs),min(ys),max(ys)
fig=plt.figure(figsize=(6.4,5.5))
axm=fig.add_axes([0.0,0.43,0.66,0.57]); axi=fig.add_axes([0.70,0.52,0.29,0.40])
axb=fig.add_axes([0.02,0.02,0.46,0.36]); axc=fig.add_axes([0.52,0.02,0.46,0.36])
T0=build(y,W0,kap)
X,_,_=paint(axm,T0,0.4,0.95)
for t in ['D1','D2','D3','U1','U2','U3']: axm.plot([W[t]],[y[base[t]]],'o',ms=2.0,color=KITE,zorder=6)
bw=max(4.0*g,0.10); bh=bw*0.62; cx=(W['U4']+wD2)/2
axm.add_patch(plt.Rectangle((cx-bw/2,lev-bh/2),bw,bh,fill=False,lw=0.5,color='k',zorder=7))
xl=-X-0.12
axm.plot([xl,xl],[zD,zU],color=INK,lw=0.4)
L=[(v,r'$y_{%d}$'%k) for k,v in enumerate(y)]+[(zU,r'$z_U$'),(zD,r'$z_D$')]
L.sort(); sep=0.092; pos=[v for v,_ in L]
for _ in range(200):
    moved=False
    for i in range(len(pos)-1):
        if pos[i+1]-pos[i]<sep:
            c=(pos[i]+pos[i+1])/2; pos[i]=c-sep/2; pos[i+1]=c+sep/2; moved=True
    if not moved: break
for (v,lab),pv in zip(L,pos):
    axm.plot([xl-0.04,xl],[v,v],color=INK,lw=0.4)
    if abs(pv-v)>1e-3: axm.plot([xl-0.10,xl-0.04],[pv,v],color=INK,lw=0.3)
    axm.text(xl-0.12,pv,lab,va='center',ha='right',fontsize=6.3,color=INK)
xbU=T0['U0'][2]; xbD=T0['D0'][2]
axm.text(xbU+0.02,y[0]-0.03,r'$U_0$',ha='left',va='top',fontsize=8)
axm.text(xbD+0.02,y[10]+0.03,r'$D_0$',ha='left',va='bottom',fontsize=8)
axm.set_xlim(xl-0.47,X+0.22); axm.set_ylim(min(zD,y[0])-0.05,max(zU,y[10])+0.05)
axm.set_aspect('equal'); axm.axis('off')
axm.text(0.5,-0.01,'(a)',transform=axm.transAxes,ha='center',va='top',fontsize=8)
# inset
paint(axi,T0,0.8,1.6)
axi.plot([W['U4']],[lev],'o',ms=3.2,mfc='white',mec=KITE,mew=0.8,zorder=6); axi.plot([wD2],[lev],'o',ms=3.0,color=KITE,zorder=6)
axi.annotate('',xy=(W['U4'],lev-0.10*bh),xytext=(wD2,lev-0.10*bh),arrowprops=dict(arrowstyle='<->',lw=0.5,shrinkA=0,shrinkB=0),zorder=7)
axi.text(cx,lev-0.15*bh,'%.3f'%g,ha='center',va='top',fontsize=7,zorder=7,bbox=dict(fc='white',ec='none',pad=0.6))
axi.text(W['U4']-0.03*bw,lev-0.06*bh,'corner of $U_4$',ha='right',va='top',fontsize=6.5,zorder=7,bbox=dict(fc='white',ec='none',pad=0.6))
axi.text(wD2+0.05*bw,lev-0.06*bh,'leg of $D_2$',ha='left',va='top',fontsize=6.5,zorder=7,bbox=dict(fc='white',ec='none',pad=0.6))
axi.set_xlim(cx-bw/2,cx+bw/2); axi.set_ylim(lev-bh/2,lev+bh/2)
axi.set_aspect('equal'); axi.set_xticks([]); axi.set_yticks([])
for s in axi.spines.values(): s.set_linewidth(0.6)
# (b), (c): walls in the positions of the patterns 22 and 36
for ax,(lab,wD,wU) in ((axb,('(b)',('D2',),('U2',))),(axc,('(c)',('D3','D2'),('U3','U2')))):
    T=build(y,W0,kap,wD,wU); X2,y0_,y1_=paint(ax,T,0.35,0.8)
    ax.set_xlim(-X2*1.02,X2*1.02); ax.set_ylim(y0_-0.03,y1_+0.03); ax.set_aspect('equal'); ax.axis('off')
    ax.text(0.5,-0.01,lab,transform=ax.transAxes,ha='center',va='top',fontsize=8)
plt.savefig(out+'.pdf',bbox_inches='tight',pad_inches=0.02); plt.savefig(out+'.png',dpi=220,bbox_inches='tight',pad_inches=0.02)
print('miss %.4f relative %.4f'%(g,g/W['U4']))
