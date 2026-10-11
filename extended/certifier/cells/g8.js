'use strict';
/* Chart A in gap coordinates, gauge y7 = inf.
   On the line the twelve finite points are, in order,
       y8 < y9 < y10 < zU < zD < y0 < y1 < y2 < y3 < y4 < y5 < y6.
   Known: y8, y0, y2, y3, y4, y5, through the logarithms of their consecutive gaps relative to |y0 - y8|:
       eta = (log|y2-y0|, log|y3-y2|, log|y4-y3|, log|y5-y4|) - log|y0-y8|,  eta in [-inf, inf]^4.
   Solved in turn, each from a condition affine in the new gap: y9 (Q43), zU (apexU), y1 (DOWN),
   y6 (Q21), y10 (Q59), zD (apexD).  Checked: N9 = P9 and UP, the 25 clauses (as monomials in the
   dihedral coordinates of y0..y10), and u(F2) against c0.
   Every quantity is the logarithm of a positive number, enclosed in an interval; differences are
   taken only where their sign is known.  Point labels: 0 zD, 1..11 y0..y10, 12 zU. */
const G=require('./lg.js');
const D=require('./tri.json'), U11=require('./u11.json');
const INF=8, ORDER=[9,10,11,12,0,1,2,3,4,5,6,7];
const POS=new Array(13).fill(-1); ORDER.forEach((p,i)=>{ POS[p]=i; });
const PAIRS=D.pairs;
const SCH=D['sched8_'+(process.env.SCHED||'134569_r')];
const STEPS=SCH.order.map(([n,x])=>({name:n,x,terms:D.eqterms[n].map(([s,c])=>{
  let other=null; const F=[];
  for(const k in c){ const i=+k, e=c[k], [a,b]=PAIRS[i];
    if(a===x||b===x) other= a===x? b : a;
    else if(a!==INF&&b!==INF) F.push([a,b,e]); }
  return {s,other,F}; })}));
const ZERO=[-Infinity,-Infinity], ONE0=[0,0];
const SP=x=>G.LSE(ONE0,x);
function L1mI(y){ const hi=Math.min(y[1],0), lo=y[0];
  const a=G.l1m1(hi), b=G.l1m1(lo); return [G.dn(a,Math.abs(a)+1), Math.min(0,G.up(b,Math.abs(b)+1))]; }
/* known points (sorted by position on the line) and the logs of their consecutive gaps */
class Chain{
  constructor(pts,gaps,parent){ this.pts=pts; this.gaps=gaps; this.cache=new Map(); this.parent=parent||null; }
  idx(p){ return this.pts.indexOf(p); }
  /* a bracket between points already known earlier is taken from the earlier (tighter) chain */
  br(a,b){ if(a===b) return ZERO;
    if(this.parent&&this.parent.idx(a)>=0&&this.parent.idx(b)>=0) return this.parent.br(a,b);
    let ia=this.idx(a), ib=this.idx(b); if(ia<0||ib<0) throw new Error('unknown point '+a+' '+b);
    if(ia>ib){ const t=ia; ia=ib; ib=t; }
    const key=ia*16+ib; let r=this.cache.get(key); if(r) return r;
    r=this.gaps[ia]; for(let k=ia+1;k<ib;k++) r=G.LSE(r,this.gaps[k]);
    this.cache.set(key,r); return r; }
  insert(x,p,z,zq){ const pts=this.pts.slice(), gaps=this.gaps.slice();
    if(zq){ const ip=this.idx(p); pts.splice(ip+1,0,x); gaps.splice(ip,1,z,zq); }
    else if(POS[x]>POS[this.pts[this.pts.length-1]]){ pts.push(x); gaps.push(z); }
    else { pts.unshift(x); gaps.unshift(z); }
    return new Chain(pts,gaps,this); }
}
function ratioLog(F,F0,ch){ let v=[0,0];
  for(const [a,b,e] of F) v=G.add(v,G.scale(e,ch.br(a,b)));
  for(const [a,b,e] of F0) v=G.add(v,G.scale(-e,ch.br(a,b)));
  return v; }
/* B(Z) at z = e^Z, signed log */
function Bat(t,Z){
  if(t.type==='one') return {s:1,L:ONE0};
  if(Z===-Infinity) return {s:1,L:t.C};
  if(t.type==='K') return {s:1,L:G.LSE(t.C,[Z,Z])};
  if(Z<t.C[0]) return {s:1,L:G.LDF(t.C,[Z,Z])};
  if(Z>t.C[1]) return {s:-1,L:G.LDF([Z,Z],t.C)};
  return {s:0}; }
/* One solve step.  With p, q the known neighbours of x and z = |x - p|, the condition reads
     E(z) = sum_i s_i g_i B_i(z),  B_i = e^{C_i} + z (a_i before x), e^{C_i} - z (after), 1 (a_i = inf),
   g_i = F_i/F_0.  E is affine in z and multi-affine in the g's.  (i) The root -E(0)/E' is a Mobius
   function of each g, so where E' keeps its sign its range is the hull of the corner values.
   (ii) If E has one sign at both ends of the admissible range for every g, the box is refuted. */
function solve(st,ch,zc){
  const x=st.x; let p=null, q=null;
  for(const r of ch.pts){ if(POS[r]<POS[x]) p=r; else if(POS[r]>POS[x]&&q===null) q=r; }
  if(p===null) throw new Error('no left neighbour');
  const T=st.terms.map((t,i)=>{
    const lg= i===0? ONE0 : ratioLog(t.F,st.terms[0].F,ch);
    const a=t.other; let type, C=null;
    if(a===INF) type='one';
    else if(POS[a]<POS[x]){ type='K'; C=ch.br(a,p); }
    else { type='D'; C=ch.br(p,a); }
    return {s:t.s,lg,type,C}; });
  const pq= q!==null? ch.br(p,q) : null;
  let Zlo=-Infinity, Zhi= pq? pq[1] : Infinity;
  if(zc){ Zlo=Math.max(Zlo,zc[0]); Zhi=Math.min(Zhi,zc[1]); }
  if(Zlo>Zhi) return {refuted:'order '+st.name};
  const lowB=T.map(t=>Bat(t,Zlo));
  const upExact= pq&&Zhi>=pq[1];
  const upB= upExact? T.map(t=> t.type==='one'? {s:1,L:ONE0} : t.type==='K'? {s:1,L:G.LSE(t.C,pq)}
                                  : (()=>{ const d=G.LDF(t.C,pq); return d? {s:1,L:d} : {s:1,L:ZERO}; })())
           : Zhi===Infinity? null : T.map(t=>Bat(t,Zhi));
  const dB=T.map(t=> t.type==='K'? 1 : t.type==='D'? -1 : 0);
  const free=[]; for(let i=1;i<T.length;i++) free.push(i);
  let sLo=null, sHi=null, sD=null, zlo=Infinity, zhi=-Infinity, rootOK=true, anyPos=false;
  const merge=(cur,s)=> cur===null? s : (cur===s? cur : 0);
  const Esign=(Bv,L)=>{ if(Bv.some(b=>b.s===0)) return 0;
    return G.ssum(T.map((t,i)=>({s:t.s*Bv[i].s, L:G.add([L[i],L[i]],Bv[i].L)}))).s; };
  for(let mask=0;mask<(1<<free.length);mask++){
    const L=T.map(t=>t.lg[0]); free.forEach((j,k)=>{ L[j]= (mask>>k)&1? T[j].lg[1] : T[j].lg[0]; });
    const infs=free.filter(j=>L[j]===Infinity);
    if(infs.length>1){ sLo=0; sHi=0; rootOK=false; continue; }
    let e0, e1, dsg;
    if(infs.length===1){ const j=infs[0], t=T[j];
      e0= lowB[j].s===0||lowB[j].L[1]===-Infinity? 0 : t.s*lowB[j].s;
      e1= upB? (upB[j].s===0||upB[j].L[1]===-Infinity? 0 : t.s*upB[j].s) : (dB[j]===0? 0 : t.s*dB[j]);
      dsg= dB[j]===0? 0 : t.s*dB[j];
      if(t.type==='D'){ zlo=Math.min(zlo,t.C[0]); zhi=Math.max(zhi,t.C[1]); anyPos=true; }
      else if(t.type==='K') zlo=-Infinity;
      else rootOK=false;
    } else {
      e0=Esign(lowB,L);
      const dn=G.ssum(T.map((t,i)=> dB[i]===0? null : {s:t.s*dB[i], L:[L[i],L[i]]}).filter(Boolean));
      dsg= dn.zero? 0 : dn.s;
      e1= upB? Esign(upB,L) : dsg;
      const E0=G.ssum(T.map((t,i)=>({s:t.s, L:G.add([L[i],L[i]], t.type==='one'? ONE0 : t.C)})));
      if(E0.s===0||dsg===0) rootOK=false;
      else if(E0.s!==dsg){ const lz=G.sub(E0.L,dn.L); zlo=Math.min(zlo,lz[0]); zhi=Math.max(zhi,lz[1]); anyPos=true; }
      else zlo=-Infinity;
    }
    sLo=merge(sLo,e0); sHi=merge(sHi,e1); sD=merge(sD,dsg);
  }
  if(!sD) rootOK=false;
  if(sLo&&sHi&&sLo===sHi) return {refuted:'order '+st.name};
  const signDecided= !!(sLo&&sHi&&sLo!==sHi);
  let z, undecided;
  if(rootOK){ if(!anyPos) return {refuted:'order '+st.name};
    z=[Math.max(zlo,Zlo),Math.min(zhi,Zhi)]; undecided= !signDecided&&!(zlo>-Infinity); }
  else { z=[Zlo,Zhi]; undecided=!signDecided; }
  if(z[0]>z[1]) return {refuted:'order '+st.name};
  let zq=null;
  if(pq){ zq=G.LDF(pq,z); if(!zq) return {refuted:'order '+st.name}; }
  return {chain:ch.insert(x,p,z,zq),undecided,z}; }
/* log u of four points a, b, c, d in cyclic order (labels 0..12), from the arcs between them */
function arcHasInf(x,y){ for(let k=(x+1)%13;k!==y;k=(k+1)%13) if(k===INF) return true; return false; }
function logU4(ch,a,b,c,d){
  if(a===INF) return G.neg(SP(G.sub(ch.br(c,d),ch.br(b,c))));
  if(b===INF) return G.neg(SP(G.sub(ch.br(c,d),ch.br(d,a))));
  if(c===INF) return G.neg(SP(G.sub(ch.br(a,b),ch.br(d,a))));
  if(d===INF) return G.neg(SP(G.sub(ch.br(a,b),ch.br(b,c))));
  const al=()=>ch.br(a,b), be=()=>ch.br(b,c), ga=()=>ch.br(c,d), de=()=>ch.br(d,a);
  if(arcHasInf(d,a)) return L1mI(G.neg(G.add(SP(G.sub(be(),al())),SP(G.sub(be(),ga())))));
  if(arcHasInf(b,c)) return L1mI(G.neg(G.add(SP(G.sub(de(),ga())),SP(G.sub(de(),al())))));
  if(arcHasInf(a,b)) return G.neg(G.add(SP(G.sub(ga(),de())),SP(G.sub(ga(),be()))));
  return G.neg(G.add(SP(G.sub(al(),be())),SP(G.sub(al(),de())))); }            // inf in arc (c,d)
const CH11=U11.chords;
const Q11=CH11.map(([i,j])=>[i+1,i+2,j+1,((j+1)%11)+1]);        // 13-labels of q_i, q_i+1, q_j, q_j+1
const cm11=v=>v.map((e,i)=>[i,e]).filter(([i,e])=>e!==0);
const LEFT={target:['eq',['target']],UP:['tri=',['UP.m1','UP.m2']],Q59:['eq',['Q59']],DOWN:['tri=',['DOWN.m1','DOWN.m2']],Q21:['eq',['Q21']],Q43:['eq',['Q43']]};
const COND=SCH.left.map(n=>[n,LEFT[n][0],LEFT[n][1].map(k=>cm11(U11.mono[k]))]);
for(const n in U11.clauses){ const c=U11.clauses[n]; COND.push([n, c.type==='bin'?'bin':'tri<', c.m.map(cm11)]); }
function mono(m,lu){ let s=[0,0]; for(const [c,e] of m) s=G.add(s,G.scale(e,lu(c))); return s; }
const KNOWN=SCH.known.slice().sort((a,b)=>POS[a]-POS[b]);
function evalBox(E,opt={},zc){
  // the largest known gap (|y2 - y9| for the default schedule) is the unit
  // opt.chart = j: the j-th known gap is the largest and is the unit; the box gives the other four (<= 0)
  const REF= opt.chart!==undefined? opt.chart : SCHED_REF; let g=[E[0],E[1],E[2],E[3]]; g.splice(REF,0,ONE0);
  if(opt.chart===undefined&&g.some(x=>x[1]>30)){ let r=0; for(let i=1;i<5;i++) if(g[i][0]>g[r][0]) r=i;
    const gr=g[r]; g=g.map((x,i)=> i===r? ONE0 : G.sub(x,gr)); }
  let ch=new Chain(KNOWN.slice(),g);
  let und=-1, urange=null; const unds=[], zs=[];
  for(let k=0;k<STEPS.length;k++){ const res=solve(STEPS[k],ch,zc&&zc[k]);
    if(res.refuted) return {status:'refuted',why:res.refuted};
    if(res.undecided&&und<0&&k<STEPS.length-2){ und=k; urange=res.z; }
    if(res.undecided) unds.push({k,z:res.z}); zs.push({k,z:res.z,und:!!res.undecided});
    ch=res.chain; }
  // T* region: the dihedral coordinate of y0, y1, y2 among the eleven heights (equal to u(F2) near T*)
  const lc=logU4(ch,3,4,11,1);
  if(opt.c0&&lc[1]<Math.log(opt.c0)) return {status:'delegated',why:'T*',lc};
  const cache=new Map(); const lu=c=>{ let v=cache.get(c); if(!v){ const [a,b,cc,d]=Q11[c]; v=logU4(ch,a,b,cc,d); cache.set(c,v); } return v; };
  const out={status:'open',why: und>=0? 'undecided '+STEPS[und].name : '',und,urange,unds,zs,ch,lc,cond:{}};
  for(const [n,kind,ms] of COND){
    if(n==='target'&&opt.notarget) continue;
    if((kind==='bin'||kind==='tri<')&&opt.noclauses) continue;
    const v= kind==='eq'||kind==='bin'? mono(ms[0],lu) : G.LSE(mono(ms[0],lu),mono(ms[1],lu));
    out.cond[n]=v; if(opt.full) continue;
    if(kind==='eq'||kind==='tri='){ if(v[0]>0||v[1]<0) return {status:'refuted',why:n}; }
    else if(kind==='bin'){ if(v[1]<0) return {status:'refuted',why:n}; }
    else if(v[0]>0) return {status:'refuted',why:n}; }
  return out; }
const SCHED_REF= KNOWN_REF();
function KNOWN_REF(){ return (process.env.SCHED||'134569_r').startsWith('134569')? 0 : 1; }
module.exports={evalBox,solve,Chain,STEPS,COND,logU4,POS,ORDER,INF,KNOWN};
