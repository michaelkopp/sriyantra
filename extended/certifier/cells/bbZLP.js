'use strict';
/* Branch and bound in the gap chart.  A node is a box of eta together with optional ranges zc[k]
   for the gap solved at step k.  Where a solve step cannot decide whether its root lies in the
   admissible range, the node is split on that gap; otherwise on the widest eta. */
const GC=require(process.env.CHART||'./g8.js');
const cp=require('child_process'), path=require('path');
const ZALL=process.env.ZALL!=='0', ZWD=+(process.env.ZWD||0.05), ZWR=+(process.env.ZWR||30);
const LPEVERY=+(process.env.LPEVERY||1), LPMIN=+(process.env.LPMIN||0), HANDOFF=+(process.env.HANDOFF||0.25), USELP=process.env.NOLP!=='1';
const PTS=GC.ORDER.filter(p=>p!==GC.INF);
function lpServer(){
  const py=cp.spawn(process.env.PYTHON||'python3',[path.join(__dirname,'lpserverA.py')],{stdio:['pipe','pipe','inherit']});
  py.on('exit',code=>{ if(code){ console.error('LP server ('+(process.env.PYTHON||'python3')+') exited with code '+code+': it needs numpy and scipy'); process.exit(2); } });
  let buf='', id=0; const pend=new Map();
  py.stdout.on('data',d=>{ buf+=d; let k; while((k=buf.indexOf('\n'))>=0){ const line=buf.slice(0,k); buf=buf.slice(k+1);
    const m=JSON.parse(line); const f=pend.get(m.id); pend.delete(m.id); f(m.refuted?m.how:null); } });
  return { call:(ch,opt)=>new Promise(res=>{ const i=id++; pend.set(i,res); const br={};
             for(let a=0;a<PTS.length;a++) for(let b=a+1;b<PTS.length;b++){ const v=ch.br(PTS[a],PTS[b]); br[PTS[a]+','+PTS[b]]=v.map(x=>Number.isFinite(x)?x:null); }
             py.stdin.write(JSON.stringify({id:i,br,opt})+'\n'); }),
           close:()=>py.stdin.end() }; }
const NT=process.env.NEWTON==='1'? require('./newton.js') : null;
const NW=+(process.env.NW||0.25);
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
async function run(B0,opt={}){
  const LP= USELP? lpServer() : null;
  const maxBoxes=opt.maxBoxes||1e7, minW=opt.minW||1e-3, ZW=opt.zw||+(process.env.ZW||0.5), ZDEEP=opt.zdeep||+(process.env.ZDEEP||-30), EZS=+(process.env.EZS||0.3);
  const stack=[{B:B0,zc:opt.zc0||null,d:0,w:1}], out={mass:0,closedMass:0,boxes:0,refuted:0,delegated:0,why:{},open:[],maxDepth:0,openWhy:{},zsplits:0,lpCalls:0,lpRefuted:0};
  const t0=Date.now(); let last=t0;
  while(stack.length){
    const node=stack.pop(), {B,zc,d}=node, wm=node.w||0; out.boxes++; if(d>out.maxDepth) out.maxDepth=d;
    if(opt.progress&&Date.now()-last>opt.progress){ last=Date.now(); /* progressMass */
      console.error(`  ${((last-t0)/1000).toFixed(0)}s closed ${(100*out.closedMass).toFixed(3)}% boxes ${out.boxes} lp ${out.lpRefuted}/${out.lpCalls} refuted ${out.refuted} delegated ${out.delegated} open ${out.open.length} zsplits ${out.zsplits} stack ${stack.length}`); }
    if(out.boxes>maxBoxes){ out.open.push({B,zc,d,why:'budget'}); out.mass+=wm; if(opt.onOpen) opt.onOpen(out.open[out.open.length-1],null); continue; }
    const r=GC.evalBox(B,opt,zc);
    if(opt.onEval) opt.onEval(B,zc,d,r);
    if(r.status!=='open'){ out.mass+=wm; out.closedMass+=wm; }
    if(r.status==='refuted'){ out.refuted++; out.why[r.why]=(out.why[r.why]||0)+1; continue; }
    if(r.status==='delegated'){ out.delegated++; continue; }
    // LP certificate on the brackets of the solved chain
    let ew0=0; for(let i=0;i<4;i++) ew0=Math.max(ew0,width(B[i]));
    const stuck0= ew0<minW;
    if(LP&&r.ch&&(stuck0||(d>=LPMIN&&d%LPEVERY===0))){ out.lpCalls++;
      const ho= HANDOFF>0&&r.lc&&r.lc[1]<Math.log(HANDOFF);
      const how=await LP.call(r.ch,{noclauses:!!opt.noclauses,handoff:ho});
      if(how){ out.lpRefuted++; out.refuted++; out.why[how]=(out.why[how]||0)+1; out.mass+=wm; out.closedMass+=wm; continue; } }
    // interval Newton on the two leftover equations: contract the box towards V8, or refute it
    if(NT&&r.und<0&&!zc&&B.every(x=>Number.isFinite(x[0])&&Number.isFinite(x[1]))&&Math.max(...B.map(x=>x[1]-x[0]))<NW&&!node.nt){
      const nb=NT.contract(B,opt.chart);
      if(nb==='empty'){ out.refuted++; out.why.newton=(out.why.newton||0)+1; out.mass+=wm; out.closedMass+=wm; continue; }
      if(nb){ const vr=nb.reduce((s,x,i)=>s*(x[1]-x[0])/(B[i][1]-B[i][0]),1);
        if(vr<0.8){ out.closedMass+=wm*(1-vr); stack.push({B:nb,zc,d,w:wm*vr,nt:false}); out.newton=(out.newton||0)+1; out.mass+=wm*(1-vr); continue; } }
      node.nt=true; }
    // split on an undecided gap?
    let ew=0; for(let i=0;i<4;i++) ew=Math.max(ew,width(B[i]));
    // split the range of an undecided solve step: the first one (ZALL=0), or the first one wide enough (ZALL=1)
    if(ew<=EZS){ const list= ZALL? (r.unds||[]) : (r.und>=0? [{k:r.und,z:r.urange}] : []);
      let done=false;
      for(const u of list){ const rg=u.z, w=rg[1]-rg[0];
        if(w>ZW&&rg[1]>ZDEEP&&rg[0]<-ZDEEP){ const s=zsplit(rg);
          if(s!==null){ out.zsplits++;
            const z1=(zc||[]).slice(), z2=(zc||[]).slice(); z1[u.k]=[rg[0],s]; z2[u.k]=[s,rg[1]];
            stack.push({B,zc:z1,d:d+1,w:wm/2},{B,zc:z2,d:d+1,w:wm/2}); done=true; break; } }
        if(!ZALL) break; }
      // a decided step whose root enclosure is still wide: bisect it (ZWD > 0)
      if(!done&&ZWD>0){ let best=null;
        for(const u of (r.zs||[])){ if(u.und) continue; const rg=u.z, w=rg[1]-rg[0];
          if(Number.isFinite(w)&&w>Math.max(ZWD,ZWR*ew)&&(!best||w>best.w)) best={k:u.k,rg,w}; }
        if(best){ const s=0.5*(best.rg[0]+best.rg[1]); out.zsplits++; out.zdsplits=(out.zdsplits||0)+1;
          const z1=(zc||[]).slice(), z2=(zc||[]).slice(); z1[best.k]=[best.rg[0],s]; z2[best.k]=[s,best.rg[1]];
          stack.push({B,zc:z1,d:d+1,w:wm/2},{B,zc:z2,d:d+1,w:wm/2}); done=true; } }
      if(done) continue; }
    let k=-1, w=0; for(let i=0;i<4;i++){ const wi=width(B[i]); if(wi>w){ w=wi; k=i; } }
    if(k<0||w<minW){ out.mass+=wm; out.open.push({B,zc,d,why:r.why}); out.openWhy[r.why]=(out.openWhy[r.why]||0)+1; if(opt.onOpen) opt.onOpen(out.open[out.open.length-1],r); continue; }
    const s=split(B[k]);
    const L=B.map(x=>x.slice()), R=B.map(x=>x.slice()); L[k][1]=s; R[k][0]=s;
    stack.push({B:L,zc,d:d+1,w:wm/2},{B:R,zc,d:d+1,w:wm/2});
  }
  if(LP) LP.close();
  out.seconds=(Date.now()-t0)/1000;
  return out; }
module.exports={run,width,split};
if(require.main===module){
  const a=process.argv.slice(2);
  const B0=a[0]? JSON.parse(a[0].replace(/inf/g,'1e999')) : [[-Infinity,Infinity],[-Infinity,Infinity],[-Infinity,Infinity],[-Infinity,Infinity]];
  run(B0,{c0:parseFloat(a[1]||'0.05'),minW:parseFloat(a[2]||'1e-3'),maxBoxes:parseFloat(a[3]||'1e7'),progress:20000,chart:a[5]!==undefined?+a[5]:undefined}).then(o=>{
  console.log(JSON.stringify({lp:[o.lpRefuted,o.lpCalls],boxes:o.boxes,refuted:o.refuted,delegated:o.delegated,open:o.open.length,zsplits:o.zsplits,maxDepth:o.maxDepth,seconds:o.seconds,why:o.why,openWhy:o.openWhy}));
  const fmt=x=>'['+x.map(v=>Number.isFinite(v)?v.toFixed(3):String(v)).join(',')+']';
  for(const u of o.open.slice(0,10)) console.log('open',u.d,u.why,u.B.map(fmt).join(' '),u.zc?JSON.stringify(u.zc):'');
  require('fs').writeFileSync(a[4]||'bbZ_open.json',JSON.stringify(o.open.map(u=>({B:u.B,zc:u.zc}))));
  });
}
