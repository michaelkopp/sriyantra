'use strict';
/* Points of V6 (Q21, Q43, Q59, apexU, apexD, DOWN) by Newton in the six pivot coordinates; propagation
   with the caps of D must not refute tiny boxes around them. */
const c=+(process.argv[2]||0.1), theta=+(process.argv[3]||0.5), K=(process.argv[4]||'2,5,7,8,9').split(',').map(Number), N=+(process.argv[5]||3000);
const caps={3:Math.log(c)}; for(const k of K) caps[k]=Math.log(theta);
const R=require('./fanR.js'), G=require('./lg.js');
const mid=x=>0.5*(x[0]+x[1]);
function res(b){ const box=new R.RBox(b.map(v=>[v,v]),R.REL);
  const F=n=>mid(box.f(R.Fk[n]));
  const dn=mid(G.LSE(box.f(R.Fk['DOWN.m1']),box.f(R.Fk['DOWN.m2'])));
  return [F('Q21'),F('Q43'),F('Q59'),F('apexU'),F('apexD'),dn]; }
const PIV=[6,7,5,9,4,8];
function solve(b){ b=b.slice();
  for(let it=0;it<60;it++){ const f=res(b); if(!f.every(Number.isFinite)) return null; const nf=Math.hypot(...f); if(nf<1e-11) return b;
    const J=[]; for(let i=0;i<6;i++){ const h=1e-7, bb=b.slice(); bb[PIV[i]]+=h; const g=res(bb); J.push(g.map((x,k)=>(x-f[k])/h)); }
    // solve J^T dx = -f  (J[i][k] = d f_k / d b_PIV[i])
    const A=[...Array(6)].map((_,k)=>[...Array(6)].map((_,i)=>J[i][k])), y=f.map(x=>-x);
    for(let col=0;col<6;col++){ let p=col; for(let r=col+1;r<6;r++) if(Math.abs(A[r][col])>Math.abs(A[p][col])) p=r;
      [A[col],A[p]]=[A[p],A[col]]; [y[col],y[p]]=[y[p],y[col]]; if(Math.abs(A[col][col])<1e-14) return null;
      for(let r=0;r<6;r++) if(r!==col){ const q=A[r][col]/A[col][col]; for(let k=col;k<6;k++) A[r][k]-=q*A[col][k]; y[r]-=q*y[col]; } }
    let st=1; const dx=y.map((v,i)=>v/A[i][i]); const mx=Math.max(...dx.map(Math.abs)); if(mx>1) st=1/mx;
    for(let i=0;i<6;i++) b[PIV[i]]+=st*dx[i]; }
  return null; }
let seed=11; const rnd=()=>{ seed=(seed*1103515245+12345)%2147483648; return seed/2147483648; };
const st={found:0,valid:0,inD:0,inPieces:0,refuted:0}; const pieceCount={};
const lamOf=b=>{ const C=R.chart(R.REL); return [null,...[1,2,3,4,5,6,7,8,9,10].map(j=>C.LAMC[j].reduce((s,x,i)=>s+x*b[i],0))]; };
const inR12=b=>{ const [e,d]=b; const rs=b.slice(4); return (e<=Math.log(0.1)&&d>=2&&rs.every(r=>Math.abs(r)<=2))||(e<=Math.log(0.1)&&d>=0.5&&d<=2.2&&rs.every(r=>Math.abs(r)<=0.4)); };
for(let n=0;n<N;n++){
  const b=[Math.log(c)-8*rnd(), -4+8*rnd(), -6*rnd(), -6*rnd(), ...[0,0,0,0,0,0].map(()=>-3+6*rnd())];
  const s=solve(b); if(!s) continue; st.found++;
  const lam=lamOf(s); if(!lam.slice(1).every(x=>x<=1e-12)) continue; st.valid++;
  if(!(lam[3]<=caps[3]&&K.every(k=>lam[k]<=caps[k]))) continue; st.inD++;
  if(inR12(s)) continue; st.inPieces++;
  const piece= s[1]<=0.5? 'A' : s[1]<=2? 'B' : 'C'; pieceCount[piece]=(pieceCount[piece]||0)+1;
  for(const w of [1e-8,1e-5]){ const box=R.propagate(s.map(v=>[v-w,v+w]),R.REL,{}); if(!box){ st.refuted++; if(st.refuted<5) console.log('REFUTED',w,s.map(v=>v.toFixed(3)).join(' ')); } } }
console.log(st,pieceCount);
