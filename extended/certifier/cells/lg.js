'use strict';
/* Arithmetic on logarithms of positive quantities, with outward padding.
   An interval [lo, hi] encloses log x; lo may be -inf (x may be 0), hi may be +inf.
   Library functions (exp, log, log1p, expm1) are assumed accurate to a few ulps; every
   result is padded by EPS relative to the magnitudes involved. */
const EPS=1e-14, TINY=1e-300;
const fin=Number.isFinite;
const dn=(x,m=Math.abs(x))=> fin(x)? x-m*EPS-TINY : x;
const up=(x,m=Math.abs(x))=> fin(x)? x+m*EPS+TINY : x;
/* log(e^a + e^b) */
function lse1(a,b){ if(a===-Infinity) return b; if(b===-Infinity) return a;
  if(a===Infinity||b===Infinity) return Infinity;
  const m=Math.max(a,b); return m+Math.log1p(Math.exp(Math.min(a,b)-m)); }
/* log(1 - e^s), s <= 0 */
function l1m1(s){ if(s===-Infinity) return 0; if(s>=0) return -Infinity;
  return s>-0.6931? Math.log(-Math.expm1(s)) : Math.log1p(-Math.exp(s)); }
/* log(e^a - e^b), a > b */
function ldf1(a,b){ if(b===-Infinity) return a; if(a===Infinity) return Infinity; return a+l1m1(b-a); }
function LSE(x,y){ const a=lse1(x[0],y[0]), b=lse1(x[1],y[1]); return [dn(a,Math.abs(a)+1), up(b,Math.abs(b)+1)]; }
function LSEn(xs){ let r=[-Infinity,-Infinity]; for(const x of xs) r=LSE(r,x); return r; }
/* log(e^A - e^B); null if A <= B everywhere.  Lower end -inf if the difference can vanish. */
function LDF(A,B){ if(!(A[1]>B[0])) return null;
  let lo=-Infinity; if(A[0]>B[1]){ const v=ldf1(A[0],B[1]); lo=dn(v,Math.abs(A[0])+Math.abs(v-A[0])+1); }
  const v=ldf1(A[1],B[0]); const hi=up(v,Math.abs(A[1])+Math.abs(v-A[1])+1);
  return [lo,hi]; }
function add(x,y){ const a=x[0]+y[0], b=x[1]+y[1];
  return [fin(a)? a-(Math.abs(x[0])+Math.abs(y[0]))*EPS-TINY : (Number.isNaN(a)?-Infinity:a),
          fin(b)? b+(Math.abs(x[1])+Math.abs(y[1]))*EPS+TINY : (Number.isNaN(b)?Infinity:b)]; }
const neg=x=>[-x[1],-x[0]];
const sub=(x,y)=>add(x,neg(y));
function scale(e,x){ if(e===0) return [0,0]; const a=e*x[0], b=e*x[1];
  return e>0? [dn(a),up(b)] : [dn(b),up(a)]; }
const hull=(x,y)=>[Math.min(x[0],y[0]),Math.max(x[1],y[1])];
const meet=(x,y)=>{ const a=Math.max(x[0],y[0]), b=Math.min(x[1],y[1]); return a<=b?[a,b]:null; };
/* signed sums: items {s: +1|-1, L: log-interval}; returns {s, L} with s in {+1,-1} or s = 0 (undetermined) */
function ssum(items){
  const P=[], N=[]; for(const it of items){ if(it.L[1]===-Infinity) continue; (it.s>0?P:N).push(it.L); }
  const p=P.length?LSEn(P):null, n=N.length?LSEn(N):null;
  if(!n) return p? {s:1,L:p} : {s:0,L:[-Infinity,-Infinity],zero:true};
  if(!p) return {s:-1,L:n};
  if(p[0]>n[1]) return {s:1,L:LDF(p,n)};
  if(n[0]>p[1]) return {s:-1,L:LDF(n,p)};
  return {s:0,L:[-Infinity,Math.max(p[1],n[1])]}; }
module.exports={EPS,dn,up,lse1,l1m1,ldf1,LSE,LSEn,LDF,add,sub,neg,scale,hull,meet,ssum};
