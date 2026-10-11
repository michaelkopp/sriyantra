'use strict';
/* Branch and bound for chart B in the adapted coordinates of fanR.js. */
const R=require('./fanR.js');
const TAIL=+(process.env.TAIL||40);
const WT=(process.env.WT||'1,1,1,1,1,1,1,1,1,1').split(',').map(Number);   // split weights per coordinate
function widthT(x){ const [a,b]=x; return Math.exp(b)-Math.exp(a); }     // t1, t10 enter only through g's: width in t
function width(x){ const [a,b]=x;
  if(a===-Infinity&&b===Infinity) return 1e9;
  if(a===-Infinity) return b<-TAIL? 0 : 8;
  if(b===Infinity) return a>TAIL? 0 : 8;
  return b-a; }
function splitAt(x){ const [a,b]=x;
  if(a===-Infinity&&b===Infinity) return 0;
  if(a===-Infinity) return b-Math.max(2,Math.abs(b));
  if(b===Infinity) return a+Math.max(2,Math.abs(a));
  return 0.5*(a+b); }
const SWITCH=new Set([1,4,5,6,7,8,9]), LAMJ={1:6,4:2,5:4,6:5,7:7,8:8,9:9};
function run(roots,opt={}){
  const maxBoxes=opt.maxBoxes||1e6, minW=opt.minW||1e-4;
  const stack=roots.map(r=>({b:r.b,mode:r.mode,w:1/roots.length,d:0}));
  const out={boxes:0,refuted:0,why:{},open:[],closedMass:0,maxDepth:0};
  const t0=Date.now(); let last=t0;
  while(stack.length){
    const node=stack.pop(); out.boxes++; if(node.d>out.maxDepth) out.maxDepth=node.d;
    if(opt.progress&&Date.now()-last>opt.progress){ last=Date.now();
      console.error(`  ${((last-t0)/1000).toFixed(0)}s closed ${(100*out.closedMass).toFixed(4)}% boxes ${out.boxes} open ${out.open.length} stack ${stack.length} depth ${node.d}`); }
    if(out.boxes>maxBoxes){ out.open.push({...node,why:'budget'}); continue; }
    const r=R.evaluate(node.b,node.mode,opt);
    if(r.status==='refuted'){ out.refuted++; out.why[r.why]=(out.why[r.why]||0)+1; out.closedMass+=node.w; continue; }
    const b=r.box.b;
    let k=-1, w=0; for(let i=0;i<10;i++){ const wi=WT[i]*((i===2||i===3)? 4*widthT(b[i]) : width(b[i])); if(wi>w){ w=wi; k=i; } }
    if(k<0||w<minW){ out.open.push({b,mode:node.mode,w:node.w,d:node.d,why:'small'}); continue; }
    const s=splitAt(b[k]); const L=b.slice(), Rr=b.slice(); L[k]=[b[k][0],s]; Rr[k]=[s,b[k][1]];
    let mR=node.mode;
    if(b[k][1]===Infinity&&node.mode[k]==='r'&&SWITCH.has(k)){      // upper tail of a relative coordinate: absolute mode
      const C=R.chart(node.mode); const j=LAMJ[k]; let rest=[0,0];
      for(let i=0;i<10;i++) if(i!==k&&C.M[j][i]) rest=[rest[0]+Math.min(C.M[j][i]*b[i][0],C.M[j][i]*b[i][1]),0];
      const lo=Number.isFinite(rest[0])? s+rest[0]-1e-9*(Math.abs(s+rest[0])+1) : -Infinity;
      Rr[k]=[Math.min(lo,0),0]; mR=node.mode.slice(0,k)+'a'+node.mode.slice(k+1); out.switches=(out.switches||0)+1; }
    stack.push({b:Rr,mode:mR,w:node.w/2,d:node.d+1},{b:L,mode:node.mode,w:node.w/2,d:node.d+1});
  }
  out.seconds=(Date.now()-t0)/1000; return out; }
module.exports={run,width,splitAt};
if(require.main===module){
  const a=process.argv.slice(2);
  const c0=+(a[0]||0.05), budget=+(a[1]||2e5);
  const I=[-Infinity,Infinity];
  const root=[[-Infinity,Math.log(c0)],I,[-Infinity,0],[-Infinity,0],I,I,I,I,I,I];
  if(a[2]) JSON.parse(a[2].replace(/inf/g,'1e999')).forEach((x,i)=>{ if(x) root[i]=x; });
  const o=run([{b:root,mode:R.REL}],{maxBoxes:budget,progress:10000,noclauses:process.env.NOCL==='1',nophi:process.env.NOPHI==='1'});
  console.log(JSON.stringify({switches:o.switches,boxes:o.boxes,refuted:o.refuted,open:o.open.length,closed:o.closedMass,why:o.why,maxDepth:o.maxDepth,seconds:o.seconds}));
  const fmt=x=>'['+x.map(v=>Number.isFinite(v)?v.toFixed(3):String(v)).join(',')+']';
  for(const u of o.open.slice(0,10)) console.log(u.why,u.d,u.b.map((x,i)=>(u.mode[i]==='a'?R.ANAMES[i]:R.NAMES[i])+fmt(x)).join(' '));
  if(a[3]) require('fs').writeFileSync(a[3],JSON.stringify(o.open.map(u=>({b:u.b,mode:u.mode}))));
}
