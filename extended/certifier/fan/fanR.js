'use strict';
/* Chart B in coordinates adapted to the face of the six fan diagonals.
   Fan chart zD = 0, zU = inf, y10 = 1, lambda_k = log t_k (k = 1..10).  Base coordinates
     b0 = e = lambda_3,  b2 = lambda_1,  b3 = lambda_10,
   and, for each of b1, b4..b9, either a relative coordinate (mode 'r') or the absolute lambda (mode 'a'):
     b1: d  = lambda_6 - lambda_3                 | lambda_6
     b4: r2 = lambda_2 - lambda_3                 | lambda_2
     b5: r4 = lambda_4 - lambda_6                 | lambda_4
     b6: r5 = lambda_5 + lambda_6 - lambda_3      | lambda_5
     b7: r7 = lambda_7 + lambda_6 - lambda_3      | lambda_7
     b8: r8 = lambda_8 - lambda_3                 | lambda_8
     b9: r9 = lambda_9 - lambda_3                 | lambda_9
   In the relative coordinates the linear parts of the conditions are
     Q21: r5,  Q43: r7,  Q59: r4 + r5 - r8,  apexU: r8 + r9 - r4 - r5,  apexD: r2 - r7 - r8,
     DOWN.m1: r8,  target: -r4,  UP.m2: -r7 - r8,
   so the common scale e cancels.  Each of Q21, Q43, Q59, apexU, apexD, DOWN is solved for its
   coordinate (b6, b7, b5, b9, b4, b8; coefficient 1 in either mode) by interval propagation.
   A box is refuted when propagation empties it, or when target, UP, a clause or Phi/t3 has a sign
   on it.  An upper tail of a relative coordinate is continued in the absolute mode. */
const G=require('./lg.js');
const FB=require('./fanB.js');
const FN=require('./funcs.json'), CL=require('./clauses.json');
/* caps: lambda_j <= CAP[j] (default 0, i.e. t_j <= 1); a region with smaller caps is set by setCaps */
const CAP=new Array(11).fill(0);
function setCaps(c){ for(let j=1;j<=10;j++) CAP[j]=(c&&c[j]!==undefined)? Math.min(0,c[j]) : 0; for(const C of CHARTS.values()) C.cache.clear(); }
if(process.env.FANCAP) setCaps(Object.fromEntries(process.env.FANCAP.split(',').map(s=>s.split(':').map(Number))));
const NAMES=['e','d','l1','l10','r2','r4','r5','r7','r8','r9'], ANAMES=['e','l6','l1','l10','l2','l4','l5','l7','l8','l9'];
const finm=(...xs)=>xs.reduce((s,x)=>s+(Number.isFinite(x)?Math.abs(x):0),0);
const pad=(x,m)=>{ const mm= (m===undefined||!Number.isFinite(m))? finm(x[0],x[1])+1 : m; return [G.dn(x[0],mm),G.up(x[1],mm)]; };
const meet=(x,y)=>{ const a=Math.max(x[0],y[0]), b=Math.min(x[1],y[1]); return a<=b?[a,b]:null; };
function linv(c,b,skip){ let v=[0,0]; for(let i=0;i<10;i++) if(c[i]&&i!==skip) v=G.add(v,G.scale(c[i],b[i])); return v; }
/* ---------- charts: lambda_j = sum_i M[j][i] b_i for a mode string ---------- */
const CHARTS=new Map();
function chart(mode){ let C=CHARTS.get(mode); if(C) return C;
  const z=()=>new Array(10).fill(0), M=[null]; for(let j=1;j<=10;j++) M.push(z());
  const add=(j,v,s=1)=>{ for(let i=0;i<10;i++) M[j][i]+=s*v[i]; };
  const unit=i=>{ const v=z(); v[i]=1; return v; };
  /* mode 'n' (cluster mode, coordinates b1..b9): b_i = log(-lambda_j), for lambda_j close to 0 (t_j close to 1);
     lambda_j is then not linear in b and its row of M is zero (NJ[j] = i) */
  const NJ=new Array(11).fill(-1);
  if(mode[0]==='n') throw new Error('mode n not allowed for e');
  if(mode[1]==='n'&&(mode[5]==='r'||mode[6]==='r'||mode[7]==='r')) throw new Error('mode n for b1 needs b5, b6, b7 absolute');
  if(mode[2]==='n') NJ[1]=2; else M[1][2]=1;
  if(mode[3]==='n') NJ[10]=3; else M[10][3]=1;
  M[3][0]=1;
  if(mode[1]==='r'){ add(6,M[3]); add(6,unit(1)); } else if(mode[1]==='n') NJ[6]=1; else add(6,unit(1));
  for(const [j,i] of [[2,4],[8,8],[9,9]]){ if(mode[i]==='r'){ add(j,M[3]); add(j,unit(i)); } else if(mode[i]==='n') NJ[j]=i; else add(j,unit(i)); }
  for(const [j,i] of [[5,6],[7,7]]){ if(mode[i]==='r'){ add(j,unit(i)); add(j,M[6],-1); add(j,M[3]); } else if(mode[i]==='n') NJ[j]=i; else add(j,unit(i)); }
  if(mode[5]==='r'){ add(4,unit(5)); add(4,M[6]); } else if(mode[5]==='n') NJ[4]=5; else add(4,unit(5));
  const coefOf=k=>{ const c=z(); for(let j=1;j<=10;j++) if(k[j]) for(let i=0;i<10;i++) c[i]+=k[j]*M[j][i]; return c; };
  const SC=new Map(); for(let a=1;a<=10;a++) for(let b=a+1;b<=11;b++){ const k=new Array(11).fill(0); for(let j=a;j<b;j++) k[j]=1; SC.set(a*16+b,coefOf(k)); }
  const NL=[]; for(let j=1;j<=10;j++) if(NJ[j]>=0) NL.push(j);
  C={mode,M,coefOf,SC,NJ,NL,LAMC:[null,...[1,2,3,4,5,6,7,8,9,10].map(j=>M[j])],cache:new Map()};
  CHARTS.set(mode,C); return C; }
/* forms in lambda: {k: [11], g: [[a,b,P]]} */
const Fk={}; for(const n in FN) Fk[n]=FB.form(FN[n]);
const CLS=Object.entries(CL).map(([n,c])=>[n,c.type,c.m.map(m=>FB.form(m))]);
const PHI=FB.PHI, RD=FB.RD;
class RBox{
  constructor(b,mode){ this.b=b; this.mode=mode; this.C=chart(mode); this.gc=new Map(); this.lc=[]; }
  coef(F){ let c=this.C.cache.get(F); if(!c){ c=this.C.coefOf(F.k); this.C.cache.set(F,c); } return c; }
  lamI(j){ let v=this.lc[j]; if(v) return v;
    const i=this.C.NJ[j];
    if(i>=0){ const x=this.b[i]; v=[-Math.exp(x[1])*(1+1e-14)-1e-300, x[0]===-Infinity? 0 : -Math.exp(x[0])*(1-1e-14)]; }
    else v=linv(this.C.LAMC[j],this.b);
    v=[v[0],Math.min(v[1],CAP[j])]; if(v[0]>v[1]) v=[v[1],v[1]]; this.lc[j]=v; return v; }
  linK(k,c){ // sum_j k_j lambda_j: symbolic in the base coordinates, met with the sum of the clipped lambda_j
    let v=linv(c||this.C.coefOf(k),this.b); for(const j of this.C.NL) if(k[j]) v=G.add(v,G.scale(k[j],this.lamI(j)));
    let w=[0,0]; for(let j=1;j<=10;j++) if(k[j]) w=G.add(w,G.scale(k[j],this.lamI(j)));
    const m=[Math.max(v[0],w[0]),Math.min(v[1],w[1])]; return m[0]<=m[1]? m : [m[1],m[1]]; }
  S(a,bb){ let v=linv(this.C.SC.get(a*16+bb),this.b); for(const j of this.C.NL) if(j>=a&&j<bb) v=G.add(v,this.lamI(j));
    let w=[0,0]; for(let j=a;j<bb;j++) w=G.add(w,this.lamI(j));
    const m=[Math.max(v[0],w[0]),Math.min(v[1],w[1],0)]; if(m[0]>m[1]) return [m[1],m[1]]; return m; }
  lg(a,bb){ const key=a*16+bb; let r=this.gc.get(key); if(r) return r;
    if(this.C.NL.some(j=>j>=a&&j<bb)){ r=this.lgN(a,bb); this.gc.set(key,r); return r; }
    const S=this.S(a,bb);
    if(S[0]>=0) r=[-Infinity,-Infinity];
    else { const lo=G.l1m1(S[1]), hi=G.l1m1(S[0]); r=pad([lo,hi]); if(r[1]>0) r[1]=0; }
    this.gc.set(key,r); return r; }
  /* log g(a,b) for a span with cluster-mode coordinates: -S = A + sum_j e^(b_j), A = -(sum of the other lambda_j),
     log g = log(-S) + h(S),  h(S) = log((1 - e^S)/(-S)) increasing in S, h(0-) = 0 */
  lgN(a,bb){ let v=linv(this.C.SC.get(a*16+bb),this.b); let wo=[0,0], logN=[-Infinity,-Infinity];
    for(let j=a;j<bb;j++){ const i=this.C.NJ[j]; if(i>=0) logN=G.LSE(logN,this.b[i]); else wo=G.add(wo,this.lamI(j)); }
    // A = -(sum of non-cluster lambda_j) >= 0: symbolic and clipped enclosures met
    let A=[Math.max(0,-v[1],-wo[1]),Math.min(-v[0],-wo[0])]; if(A[0]>A[1]) A=[A[1],A[1]];
    const logA=[A[0]>0? G.dn(Math.log(A[0])) : -Infinity, A[1]>0? G.up(Math.log(A[1])) : -Infinity];
    const L=G.LSE(logA,logN);                                  // log(-S)
    if(L[1]===-Infinity) return [-Infinity,-Infinity];
    const Slo=-Math.exp(L[1]), Shi=L[0]===-Infinity? 0 : -Math.exp(L[0]);
    const h=S=>{ if(S===0) return 0; if(S===-Infinity) return -Infinity; return Math.log(-Math.expm1(S)/(-S)); };
    const hr=[G.dn(h(Slo),Math.abs(h(Slo))+1e-6), Math.min(0,G.up(h(Shi),Math.abs(h(Shi))+1e-6))];
    const r=G.add(L,hr); if(r[1]>0) r[1]=0; return r; }
  gpart(F){ let v=[0,0]; for(const [a,bb,P] of F.g) v=G.add(v,G.scale(P,this.lg(a,bb))); return v; }
  f(F){ return G.add(this.linK(F.k,this.coef(F)),this.gpart(F)); }
  phiOverC(){ const items=[];
    for(const [a,bb,n] of PHI.g){ const k=new Array(11).fill(0); for(let j=a;j<bb;j++) k[j]=1; k[3]-=1;
      const lP=this.linK(k); const S=this.S(a,bb);
      const lh=[lhS(S[0]),lhS(S[1])];
      const it={s:-Math.sign(n),L:pad(G.add(lP,lh))}; for(let i=0;i<Math.abs(n);i++) items.push(it); }
    const k=RD.k.slice(); k[3]-=1; const lR=G.add(this.linK(k),this.gpart(RD));
    const r=G.add(lR,this.lamI(3));
    items.push({s:1,L:pad(G.add(lR,[lhtR(r[1]),lhtR(r[0])]))});
    return G.ssum(items); }
}
function lhS(S){ if(S===-Infinity) return 0; if(S>=0) return Infinity; return Math.log(-G.l1m1(S))-S; }
function lhtR(r){ if(r===-Infinity) return 0; if(r===Infinity) return -Infinity; return Math.log(G.lse1(0,r))-r; }
/* pivots: condition = 0 solved for one base coordinate (coefficient 1 in both modes) */
const PIV=[['Q21',6],['Q43',7],['Q59',5],['apexU',9],['apexD',4],['DOWN',8]];
function solveFor(box,F,i,rhs){ const c=box.coef(F); if(c[i]!==1) throw new Error('pivot coefficient '+c[i]);
  let rest=G.add(linv(c,box.b,i),box.gpart(F));
  for(const j of box.C.NL) if(F.k[j]) rest=G.add(rest,G.scale(F.k[j],box.lamI(j)));   // cluster-mode lambda_j
  return G.sub(rhs,rest); }
function propagate(b,mode,opt={}){
  let box=new RBox(b,mode);
  const sweeps=opt.sweeps||6;
  for(let s=0;s<sweeps;s++){
    let changed=0;
    const nb=box.b.slice();
    for(const j of box.C.NL) if(CAP[j]<0){ const i=box.C.NJ[j], lb=G.dn(Math.log(-CAP[j]),1);    // cluster mode: e^(b_i) >= -CAP[j]
      if(lb>nb[i][0]){ nb[i]=[lb,nb[i][1]]; changed++; if(!(nb[i][0]<=nb[i][1])) return null; } }
    for(let j=1;j<=10;j++){ const c=box.C.LAMC[j]; for(let i=0;i<10;i++) if(Math.abs(c[i])===1){   // lambda_j <= CAP[j] (<= 0)
        const rest=linv(c,nb,i);
        if(c[i]>0){ const ub=G.up(CAP[j]-rest[0],Math.abs(CAP[j])+Math.abs(rest[0])+1); if(ub<nb[i][1]){ nb[i]=[nb[i][0],ub]; changed++; } }
        else { const lb=G.dn(rest[0]-CAP[j],Math.abs(CAP[j])+Math.abs(rest[0])+1); if(lb>nb[i][0]){ nb[i]=[lb,nb[i][1]]; changed++; } }
        if(!(nb[i][0]<=nb[i][1])) return null; } }
    box=new RBox(nb,mode);
    for(const [n,i] of PIV){ let v; if(mode[i]==='n') continue;
      if(n==='DOWN'){ const m2=box.f(Fk['DOWN.m2']); if(m2[0]>=0) return null;
        const m1=[m2[1]>=0?-Infinity:G.l1m1(m2[1]), G.l1m1(m2[0])];
        v=solveFor(box,Fk['DOWN.m1'],i,pad(m1)); }
      else v=solveFor(box,Fk[n],i,[0,0]);
      v=pad(v);
      const m=meet(box.b[i],v); if(!m) return null;
      if(m[0]>box.b[i][0]+1e-9||m[1]<box.b[i][1]-1e-9){ changed++; const nb2=box.b.slice(); nb2[i]=m; box=new RBox(nb2,mode); } }
    if(!changed) break; }
  return box; }
function evaluate(b,mode,opt={}){
  const box=propagate(b,mode,opt); if(!box) return {status:'refuted',why:'prop'};
  const t=box.f(Fk.target); if(t[0]>0||t[1]<0) return {status:'refuted',why:'target',box};
  const up=G.LSE(box.f(Fk['UP.m1']),box.f(Fk['UP.m2'])); if(up[0]>0||up[1]<0) return {status:'refuted',why:'UP',box};
  if(opt.c0){ const u=box.f(FB.U11); if(u[0]>=Math.log(opt.c0)) return {status:'refuted',why:'chartA',box}; }
  if(!opt.nophi){ const p=box.phiOverC(); if(p.s!==0) return {status:'refuted',why:'Phi'+(p.s>0?'+':'-'),box}; }
  if(!opt.noclauses) for(const [n,type,fs] of CLS){
    if(type==='bin'){ const v=box.f(fs[0]); if(v[1]<0) return {status:'refuted',why:n,box}; }
    else { const v=G.LSE(box.f(fs[0]),box.f(fs[1])); if(v[0]>0) return {status:'refuted',why:n,box}; } }
  return {status:'open',box}; }
/* base coordinates of a configuration t[1..10] in the relative mode */
const REL='rrrrrrrrrr';
function toBase(t){ const l=t.map(Math.log); const L=j=>l[j-1];
  const e=L(3), d=L(6)-L(3);
  return [e,d,L(1),L(10),L(2)-e,L(4)-L(6),L(5)+d,L(7)+d,L(8)-e,L(9)-e]; }
module.exports={NAMES,ANAMES,RBox,propagate,evaluate,toBase,chart,Fk,PHI,RD,REL,setCaps,CAP};
