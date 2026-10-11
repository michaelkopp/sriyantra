// accurate conversion of fan-chart l (l_k = log t_k) to the gap coordinates eta
function logX(l){ const L=new Array(12).fill(0); for(let k=10;k>=1;k--) L[k]=l[k-1]+L[k+1]; return L; } // log p_k, k=1..11
function logGap(l,a,b){ // log(p_b - p_a), 1<=a<b<=11 : log p_b + log(1 - p_a/p_b)
  const L=logX(l); let S=0; for(let j=a;j<b;j++) S+=l[j-1]; return L[b]+Math.log(-Math.expm1(S)); }
function toEta(l){ const s=logGap(l,3,5); return [logGap(l,1,2)-s,logGap(l,2,3)-s,logGap(l,5,6)-s,logGap(l,6,8)-s]; }
module.exports={toEta,logGap};
// gauge y7 = inf: log|b - a|' = log|b - a| - log|a - y7| - log|b - y7|  (fan labels: y_k = p_{k+1})
function lg8(l,a,b){ const c=8, L=(u,v)=>logGap(l,Math.min(u,v),Math.max(u,v)); return L(a,b)-L(a,c)-L(b,c); }
function toEta8(l){ const s=lg8(l,9,1); return [lg8(l,1,3)-s,lg8(l,3,4)-s,lg8(l,4,5)-s,lg8(l,5,6)-s]; }
module.exports.toEta8=toEta8; module.exports.lg8=lg8;
function toEtaS(l,known,ref){ // gaps between consecutive known points (chain order), minus the reference gap
  const POS=[4,5,6,7,8,9,10,11,-1,0,1,2,3]; const k=known.slice().sort((a,b)=>POS[a]-POS[b]);
  const g=[]; for(let i=0;i<k.length-1;i++) g.push(lg8(l,k[i],k[i+1]));
  const r=g[ref]; return g.filter((_,i)=>i!==ref).map(x=>x-r); }
module.exports.toEtaS=toEtaS;
