'use strict';
/* enclosure test for cluster mode: random configurations (some t_j within 1e-9..1e-2 of 1), boxes in
   modes with 'n' coordinates around them; every form enclosure must contain the point value. */
const R=require('./fanR.js'), FB=require('./fanB.js');
let seed=5; const rnd=()=>{ seed=(seed*1103515245+12345)%2147483648; return seed/2147483648; };
const forms=Object.values(R.Fk).concat(FB.CLS.flatMap(c=>c[2]),[FB.PHI,FB.RD,FB.U11]);
function val(F,lam){ let v=0; for(let j=1;j<=10;j++) v+=F.k[j]*lam[j];
  for(const [a,b,P] of F.g){ let S=0; for(let j=a;j<b;j++) S+=lam[j]; v+=P*Math.log(-Math.expm1(S)); } return v; }
// base coordinates (absolute where needed) from lambda for a given mode
function base(lam,mode){ const e=lam[3]; const b=new Array(10);
  b[0]=e; b[2]=mode[2]==='n'?Math.log(-lam[1]):lam[1]; b[3]=mode[3]==='n'?Math.log(-lam[10]):lam[10];
  const L6=lam[6];
  b[1]= mode[1]==='r'? L6-e : mode[1]==='n'? Math.log(-L6) : L6;
  for(const [j,i] of [[2,4],[8,8],[9,9]]) b[i]= mode[i]==='r'? lam[j]-e : mode[i]==='n'? Math.log(-lam[j]) : lam[j];
  for(const [j,i] of [[5,6],[7,7]]) b[i]= mode[i]==='r'? lam[j]+L6-e : mode[i]==='n'? Math.log(-lam[j]) : lam[j];
  b[5]= mode[5]==='r'? lam[4]-L6 : mode[5]==='n'? Math.log(-lam[4]) : lam[4];
  return b; }
let bad=0, n=0;
for(let it=0;it<3000;it++){
  const lam=[null]; const near=[];
  for(let j=1;j<=10;j++){ const u=rnd(); if(u<0.4){ const d=Math.exp(-2-20*rnd()); lam.push(-d); near.push(j); } else lam.push(-Math.exp(-3+4*rnd())); }
  // mode: n for near-1 coordinates where allowed, else a (or r for a few)
  const JI={6:1,1:2,10:3,2:4,4:5,5:6,7:7,8:8,9:9};
  let mode=['r','a','r','r','a','a','a','a','a','a'];
  for(const j of near) if(JI[j]!==undefined) mode[JI[j]]='n';
  if(mode[1]==='n'){ for(const i of [5,6,7]) if(mode[i]==='r') mode[i]='a'; }
  if(rnd()<0.3&&mode[8]!=='n') mode[8]='r';
  mode=mode.join('');
  const b=base(lam,mode); if(!b.every(Number.isFinite)) continue;
  const w=1e-9*(1+3*rnd()); const box=new R.RBox(b.map(x=>[x-w*(1+Math.abs(x)),x+w*(1+Math.abs(x))]),mode);
  for(const F of forms){ const v=val(F,lam), e=box.f(F); n++;
    if(!(e[0]<=v+1e-9*(1+Math.abs(v))&&v-1e-9*(1+Math.abs(v))<=e[1])&&Number.isFinite(v)){ bad++; if(bad<5) console.log('BAD',mode,v,e); } } }
console.log({checks:n,bad});
// enclosure widths: cluster mode vs absolute mode for configurations with t_5..t_9 within 1e-6 of 1
{ const lam=[null,-0.5,-2,-3,-0.3,-1e-6,-2e-6,-5e-7,-3e-6,-1e-6,-1];
  for(const mode of ['raraaaaaaa','raraaannnn','rnraannnnn']){ let b=base(lam,mode); const box=new R.RBox(b.map(x=>[x-1e-6*(1+Math.abs(x)),x+1e-6*(1+Math.abs(x))]),mode);
    const ws=forms.map(F=>{ const e=box.f(F); return e[1]-e[0]; }); ws.sort((a,b)=>a-b);
    console.log(mode,'median width',ws[ws.length>>1].toExponential(2),'max',ws[ws.length-1].toExponential(2)); } }
