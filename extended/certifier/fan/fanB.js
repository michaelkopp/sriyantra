'use strict';
/* Chart B: the region u(y0,y1,y2 | y3..y10) < c0 left by chart A, in the fan chart
     zD = 0, zU = inf, y10 = 1,   p_k = t_k t_{k+1} ... t_10,   t_k in (0, 1)   (k = 1..10).
   Every condition and clause is  f = sum_k kappa_k log t_k + sum_(a,b) P_ab log g(a,b),
   g(a,b) = 1 - t_a ... t_(b-1)  (funcs.json, clauses.json).  A box gives, for each k, one of
     'L'  lambda_k = log t_k            in [-inf, -ln 2]
     'M'  mu_k = log(-log t_k)          in [-inf, log ln 2]      (t_k near 1: 1 - t_k ~ e^mu)
     'S'  sigma_k = log(t_k / t_3)      (fan coordinates, with t_3 in 'L')
     'F'  undetermined, t_k in (0, 1)   (split into an 'L' and an 'M' part)
   Linear forms are evaluated with the fan coordinates relative to t_3, so that the common scale
   log t_3 cancels in every balanced condition.  A box is refuted when
     - u11 >= c0 on the whole box (chart A covers it), or
     - one of the seven conditions or the target is nonzero on the whole box, or
     - one of the 25 clauses fails on the whole box, or
     - Phi / t_3 has a sign on the whole box, where
         Phi = f_target - f_Q21 + f_Q59 + log(m1 + m2)_DOWN = sum_g n_g log g + log(1 + R),
       R = m2/m1 (DOWN); Phi = 0 at every point of V8, and Phi/t_3 is evaluated as a signed sum of
       the terms -n_g (P_g/t_3) h(P_g), (R/t_3) h~(R) with h(P) = -log(1-P)/P, h~(R) = log(1+R)/R. */
const G=require('./lg.js');
const FN=require('./funcs.json'), CL=require('./clauses.json');
const LN2=Math.log(2), LLN2=Math.log(LN2), FAN=[2,5,7,8,9];
const fin=Number.isFinite;
const finm=(...xs)=>xs.reduce((s,x)=>s+(Number.isFinite(x)?Math.abs(x):0),0);
const pad=(x,m)=>{ const mm= (m===undefined||!Number.isFinite(m))? finm(x[0],x[1])+1 : m; return [G.dn(x[0],mm),G.up(x[1],mm)]; };
/* ---------- forms ---------- */
function form(f){ const k=new Array(11).fill(0); f.kappa.forEach((v,i)=>k[i+1]=v);
  return {k,g:f.terms.map(t=>[t.a,t.b,t.P])}; }
function comb(list){ const k=new Array(11).fill(0), g=new Map();
  for(const [c,f] of list){ for(let j=1;j<=10;j++) k[j]+=c*f.k[j]; for(const [a,b,P] of f.g){ const key=a*16+b; g.set(key,(g.get(key)||0)+c*P); } }
  return {k,g:[...g].filter(([key,P])=>P!==0).map(([key,P])=>[key>>4,key&15,P])}; }
const Fm={}; for(const n in FN) Fm[n]=form(FN[n]);
const EQS=[['Q21',[Fm.Q21]],['Q43',[Fm.Q43]],['Q59',[Fm.Q59]],['apexU',[Fm.apexU]],['apexD',[Fm.apexD]],
  ['UP',[Fm['UP.m1'],Fm['UP.m2']]],['DOWN',[Fm['DOWN.m1'],Fm['DOWN.m2']]],['target',[Fm.target]]];
const CLS=Object.entries(CL).map(([n,c])=>[n,c.type,c.m.map(form)]);
const PHI=comb([[1,Fm.target],[-1,Fm.Q21],[1,Fm.Q59],[1,Fm['DOWN.m1']]]);
if(PHI.k.some(x=>x!==0)) throw new Error('Phi is not balanced');
const RD=comb([[1,Fm['DOWN.m2']],[-1,Fm['DOWN.m1']]]);           // log R, R = m2/m1 of DOWN
const U11=comb([[1,{k:[0,0,0,1,0,0,0,0,0,0,0],g:[[1,3,1],[4,11,1],[1,4,-1],[3,11,-1]]}]]);
/* ---------- evaluation of a box ---------- */
function l1m1I(S){ // log(1 - e^S) for S in [S0, S1] (S1 <= 0)
  const a=G.l1m1(S[1]), b=G.l1m1(S[0]); return pad([a,b],Math.abs(a)+Math.abs(S[1])+1); }
function phiS(S){ // log((1 - e^S)/(-S)), increasing in S; -> 0 as S -> 0
  if(S===0) return 0; if(S===-Infinity) return -Infinity; return Math.log(-Math.expm1(S)/(-S)); }
class Box{
  constructor(mode,iv){ this.mode=mode; this.iv=iv; this.lam=[]; this.mu=[]; this.gc=new Map(); this.bad=false;
    for(let j=1;j<=10;j++){ const m=mode[j], x=iv[j]; let lam, mu;
      if(m==='L'){ lam=x.slice(); }
      else if(m==='M'){ mu=x.slice(); lam=[-Math.exp(x[1]),-Math.exp(x[0])]; lam=pad(lam,Math.abs(lam[0])); }
      else if(m==='F'){ lam=[-Infinity,0]; mu=[-Infinity,Infinity]; }
      else if(m==='S'){ continue; }
      this.lam[j]=lam; this.mu[j]=mu||this.muOf(lam); }
    for(let j=1;j<=10;j++) if(mode[j]==='S'){ let lam=G.add(this.lam[3],iv[j]); if(lam[0]>=0){ this.bad=true; lam=[-1e-300,0]; } if(lam[1]>0) lam[1]=0;
      this.lam[j]=lam; this.mu[j]=this.muOf(lam); } }
  muOf(lam){ // mu = log(-lambda)
    const a= lam[1]>=0? -Infinity : Math.log(-lam[1]), b= lam[0]===-Infinity? Infinity : Math.log(-lam[0]);
    return pad([a,b],Math.abs(a)+Math.abs(b)+1); }
  lin(k){ // sum_j k_j lambda_j with fan coordinates relative to lambda_3
    let c3=k[3]||0, v=[0,0];
    for(let j=1;j<=10;j++){ if(j===3||!k[j]) continue;
      if(this.mode[j]==='S'){ c3+=k[j]; v=G.add(v,G.scale(k[j],this.iv[j])); } else v=G.add(v,G.scale(k[j],this.lam[j])); }
    if(c3) v=G.add(v,G.scale(c3,this.lam[3]));
    return v; }
  S(a,b){ const k=new Array(11).fill(0); for(let j=a;j<b;j++) k[j]=1; const v=this.lin(k); if(v[1]>0) v[1]=0; return v; }
  lg(a,b){ const key=a*16+b; let r=this.gc.get(key); if(r) return r;
    const S=this.S(a,b);
    r=l1m1I(S);
    // near S = 0: log g = log(-S) + phi(S), log(-S) = LSE of mu_j
    let lms=[-Infinity,-Infinity]; for(let j=a;j<b;j++) lms=G.LSE(lms,this.mu[j]);
    const r2=pad(G.add(lms,[phiS(S[0]),phiS(S[1])]));
    const m=G.meet(r,r2); if(m) r=m;
    this.gc.set(key,r); return r; }
  f(F){ let v=this.lin(F.k); for(const [a,b,P] of F.g) v=G.add(v,G.scale(P,this.lg(a,b))); return v; }
  /* Phi / t_3 as a signed sum */
  phiOverC(){ const items=[];
    for(const [a,b,n] of PHI.g){ const k=new Array(11).fill(0); for(let j=a;j<b;j++) k[j]=1; k[3]-=1;
      const lP=this.lin(k);                         // log(P/t3)
      const S=this.S(a,b);
      const lh=[lhS(S[0]),lhS(S[1])];               // log h(P), increasing in S
      items.push({s:-Math.sign(n),L:pad(G.add(lP,lh)),n:Math.abs(n)}); }
    const k=RD.k.slice(); k[3]-=1; let lR=this.lin(k); for(const [a,b,P] of RD.g) lR=G.add(lR,G.scale(P,this.lg(a,b)));
    const r=G.add(lR,this.lam[3]);                  // log R
    items.push({s:1,L:pad(G.add(lR,[lhtR(r[1]),lhtR(r[0])])),n:1});
    const ex=[]; for(const it of items){ for(let i=0;i<it.n;i++) ex.push(it); }
    return G.ssum(ex); }
}
function lhS(S){ if(S===-Infinity) return 0; if(S>=0) return Infinity; return Math.log(-G.l1m1(S))-S; }
function lhtR(r){ if(r===-Infinity) return 0; if(r===Infinity) return -Infinity; return Math.log(G.lse1(0,r))-r; }
function evaluate(box,opt={}){
  if(box.bad) return {status:'refuted',why:'range'};
  const c0=opt.c0||0.05;
  const u=box.f(U11); if(u[0]>=Math.log(c0)) return {status:'refuted',why:'chartA'};
  const vals={};
  for(const [n,fs] of EQS){ const v= fs.length===1? box.f(fs[0]) : G.LSE(box.f(fs[0]),box.f(fs[1]));
    vals[n]=v; if(v[0]>0||v[1]<0) return {status:'refuted',why:n}; }
  if(!opt.noclauses) for(const [n,type,fs] of CLS){
    if(type==='bin'){ const v=box.f(fs[0]); if(v[1]<0) return {status:'refuted',why:n}; }
    else { const v=G.LSE(box.f(fs[0]),box.f(fs[1])); if(v[0]>0) return {status:'refuted',why:n}; } }
  if(!opt.nophi){ const p=box.phiOverC(); if(p.s!==0) return {status:'refuted',why:'Phi'+(p.s>0?'+':'-')}; }
  return {status:'open',vals}; }
module.exports={Box,evaluate,EQS,CLS,PHI,RD,U11,comb,form,LN2,LLN2,FAN};
