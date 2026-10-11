'use strict';
/* Chart B certification of the region D = { t3 <= c, t_k <= theta (k in K) } outside the boxes already
   certified (R1 = {t3 <= 0.1, d >= 2, |r| <= 2}, R2 = {t3 <= 0.1, 0.5 <= d <= 2.2, |r| <= 0.4}).
     node fanD.js C THETA BUDGET [K=2,5,7,8,9] [OUTFILE]
   Every point of D lies in R1, R2, or one of the pieces
     A: d <= 0.5;   B(i,+-): 0.5 <= d <= 2, +-r_i >= 0.4;   C(i,+-): d >= 2, +-r_i >= 2,
   and the run refutes every piece (interval tests of fanR.js with the caps of D). */
const c=+(process.argv[2]||0.1), theta=+(process.argv[3]||0.5), budget=+(process.argv[4]||2e6);
const K=(process.argv[5]||'2,5,7,8,9').split(',').map(Number);
const caps={3:Math.log(c)}; for(const k of K) caps[k]=Math.log(theta);
const R=require('./fanR.js'); R.setCaps(caps);
const BB=require('./fanRBB.js');
const I=[-Infinity,Infinity], E=Math.log(c);
const base=()=>[[-Infinity,E],I,[-Infinity,0],[-Infinity,0],I,I,I,I,I,I];
const roots=[]; let b=base(); b[1]=[-Infinity,0.5]; roots.push({b,mode:R.REL,tag:'A'});
for(let i=4;i<10;i++) for(const [s,rg] of [['+',[0.4,Infinity]],['-',[-Infinity,-0.4]]]){ b=base(); b[1]=[0.5,2]; b[i]=rg; roots.push({b,mode:R.REL,tag:'B'+R.NAMES[i]+s}); }
for(let i=4;i<10;i++) for(const [s,rg] of [['+',[2,Infinity]],['-',[-Infinity,-2]]]){ b=base(); b[1]=[2,Infinity]; b[i]=rg; roots.push({b,mode:R.REL,tag:'C'+R.NAMES[i]+s}); }
const only=process.env.ONLY? new Set(process.env.ONLY.split(',')) : null;
const fmt=x=>'['+x.map(v=>Number.isFinite(v)?v.toFixed(3):String(v)).join(',')+']';
let allClosed=true; const t0=Date.now(); const openAll=[];
for(const r of roots){ if(only&&!only.has(r.tag)) continue;
  const o=BB.run([{b:r.b,mode:r.mode}],{maxBoxes:budget,progress:+(process.env.PROGRESS||0)});
  if(o.open.length) allClosed=false;
  console.log(JSON.stringify({piece:r.tag,boxes:o.boxes,open:o.open.length,why:o.why,maxDepth:o.maxDepth,seconds:o.seconds}));
  for(const u of o.open.slice(0,3)) console.log('   ',u.why,u.d,u.b.map((x,i)=>(u.mode[i]==='a'?R.ANAMES[i]:R.NAMES[i])+fmt(x)).join(' '));
  for(const u of o.open) openAll.push({piece:r.tag,b:u.b.map(x=>x.map(v=>Number.isFinite(v)?v:(v>0?'inf':'-inf'))),mode:u.mode,why:u.why}); }
console.log(JSON.stringify({c,theta,K,closed:allClosed,seconds:(Date.now()-t0)/1000}));
if(process.argv[6]) require('fs').writeFileSync(process.argv[6],JSON.stringify(openAll));
