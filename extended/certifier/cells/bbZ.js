'use strict';
/* Branch and bound in the gap chart.  A node is a box of eta together with optional ranges zc[k]
   for the gap solved at step k.  Where a solve step cannot decide whether its root lies in the
   admissible range, the node is split on that gap; otherwise on the widest eta. */
const GC=require(process.env.CHART||'./g8.js');
const TAIL=40;
function width(x){ const [a,b]=x;
  if(a===-Infinity&&b===Infinity) return 1e9;
  if(a===-Infinity) return b<-TAIL? 0 : 8;
  if(b===Infinity) return a>TAIL? 0 : 8;
  return b-a; }
function split(x){ const [a,b]=x;
  if(a===-Infinity&&b===Infinity) return 0;
  if(a===-Infinity) return b-Math.max(2,Math.abs(b));
  if(b===Infinity) return a+Math.max(2,Math.abs(a));
  return 0.5*(a+b); }
function zsplit(r){ const [a,b]=r;
  if(a===-Infinity&&b===Infinity) return 0;
  if(a===-Infinity) return b-6;
  if(b===Infinity) return a+6;
  if(b-a<1e-3) return null;
  return 0.5*(a+b); }
function run(B0,opt={}){
  const maxBoxes=opt.maxBoxes||1e7, minW=opt.minW||1e-3, ZW=opt.zw||+(process.env.ZW||0.05), ZDEEP=opt.zdeep||+(process.env.ZDEEP||-30), EZS=+(process.env.EZS||0.1);
  const stack=[{B:B0,zc:null,d:0,w:1}], out={mass:0,boxes:0,refuted:0,delegated:0,why:{},open:[],maxDepth:0,openWhy:{},zsplits:0};
  const t0=Date.now(); let last=t0;
  while(stack.length){
    const node=stack.pop(), {B,zc,d}=node, wm=node.w||0; out.boxes++; if(d>out.maxDepth) out.maxDepth=d;
    if(opt.progress&&Date.now()-last>opt.progress){ last=Date.now();
      console.error(`  ${((last-t0)/1000).toFixed(0)}s done ${(100*out.mass).toFixed(3)}% boxes ${out.boxes} refuted ${out.refuted} delegated ${out.delegated} open ${out.open.length} zsplits ${out.zsplits} stack ${stack.length}`); }
    if(out.boxes>maxBoxes){ out.open.push({B,zc,d,why:'budget'}); out.mass+=wm; continue; }
    const r=GC.evalBox(B,opt,zc);
    if(opt.onEval) opt.onEval(B,zc,d,r);
    if(r.status!=='open') out.mass+=wm;
    if(r.status==='refuted'){ out.refuted++; out.why[r.why]=(out.why[r.why]||0)+1; continue; }
    if(r.status==='delegated'){ out.delegated++; continue; }
    // split on an undecided gap?
    let ew=0; for(let i=0;i<4;i++) ew=Math.max(ew,width(B[i]));
    if(r.und>=0&&ew<=EZS){ const rg=r.urange, w=rg[1]-rg[0];
      if(w>ZW&&rg[1]>ZDEEP&&rg[0]<-ZDEEP){ const s=zsplit(rg);
        if(s!==null){ out.zsplits++;
          const z1=(zc||[]).slice(), z2=(zc||[]).slice(); z1[r.und]=[rg[0],s]; z2[r.und]=[s,rg[1]];
          stack.push({B,zc:z1,d:d+1,w:wm/2},{B,zc:z2,d:d+1,w:wm/2}); continue; } } }
    let k=-1, w=0; for(let i=0;i<4;i++){ const wi=width(B[i]); if(wi>w){ w=wi; k=i; } }
    if(k<0||w<minW){ out.mass+=wm; out.open.push({B,zc,d,why:r.why}); out.openWhy[r.why]=(out.openWhy[r.why]||0)+1; continue; }
    const s=split(B[k]);
    const L=B.map(x=>x.slice()), R=B.map(x=>x.slice()); L[k][1]=s; R[k][0]=s;
    stack.push({B:L,zc,d:d+1,w:wm/2},{B:R,zc,d:d+1,w:wm/2});
  }
  out.seconds=(Date.now()-t0)/1000;
  return out; }
module.exports={run,width,split};
if(require.main===module){
  const a=process.argv.slice(2);
  const B0=a[0]? JSON.parse(a[0].replace(/inf/g,'1e999')) : [[-Infinity,Infinity],[-Infinity,Infinity],[-Infinity,Infinity],[-Infinity,Infinity]];
  const o=run(B0,{c0:parseFloat(a[1]||'0.05'),minW:parseFloat(a[2]||'1e-3'),maxBoxes:parseFloat(a[3]||'1e7'),progress:20000});
  console.log(JSON.stringify({boxes:o.boxes,refuted:o.refuted,delegated:o.delegated,open:o.open.length,zsplits:o.zsplits,maxDepth:o.maxDepth,seconds:o.seconds,why:o.why,openWhy:o.openWhy}));
  const fmt=x=>'['+x.map(v=>Number.isFinite(v)?v.toFixed(3):String(v)).join(',')+']';
  for(const u of o.open.slice(0,10)) console.log('open',u.d,u.why,u.B.map(fmt).join(' '),u.zc?JSON.stringify(u.zc):'');
  require('fs').writeFileSync(a[4]||'bbZ_open.json',JSON.stringify(o.open.map(u=>({B:u.B,zc:u.zc}))));
}
