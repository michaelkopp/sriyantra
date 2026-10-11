'use strict';
/* Chart A over a grid of cells.  Each cell is closed by branch and bound with a budget; the result of
   every cell is appended to results.jsonl, and cells already there are skipped when the run is resumed.
     node runcells.js --jobs 10 --budget 2e5 [--charts 0,1,2,3,4] [--edges -inf,-8,...,0] [--out results2.jsonl]
   A cell is 'closed' when every box in it is refuted or lies in u(F2) < c0 (left to the chart of T*). */
const fs=require('fs'), path=require('path'), cp=require('child_process'), os=require('os');
const args=process.argv.slice(2), opt={};
for(let i=0;i<args.length;i++) if(args[i].startsWith('--')) opt[args[i].slice(2)]=args[i+1];
const OUT=opt.out||'results2.jsonl', JOBS=+(opt.jobs||Math.max(1,os.cpus().length-1)), BUDGET=+(opt.budget||2e5);
const C0=+(opt.c0||0.05), MINW=+(opt.minw||1e-3);
if(opt.worker){                                   // worker: read cells from stdin lines, write results
  process.env.SCHED=process.env.SCHED||'134569_r';
  const bb=require('./bbZ.js');
  process.on('message',m=>{
    if(m==='exit') process.exit(0);
    const t=Date.now();
    const o=bb.run(m.B,{c0:C0,minW:MINW,maxBoxes:BUDGET,chart:m.chart});
    const openWhy={}; for(const u of o.open) openWhy[u.why]=(openWhy[u.why]||0)+1;
    process.send({id:m.id,chart:m.chart,B:m.B,status:o.open.length?'open':'closed',boxes:o.boxes,refuted:o.refuted,
      delegated:o.delegated,open:o.open.length,openWhy,mass:o.mass,why:o.why,seconds:(Date.now()-t)/1000});
  });
  return;
}
/* five charts: in chart j the j-th of the five known gaps is the largest and is the unit, and the
   box gives the logarithms of the other four relative to it, each in [-inf, 0] */
const edges=(opt.edges||'-inf,-8,-6,-4,-3,-2,-1.5,-1,-0.5,0').split(',').map(x=> x==='-inf'?-Infinity: x==='inf'?Infinity:+x);
const charts=(opt.charts||'0,1,2,3,4').split(',').map(Number);
const cells=[]; const n=edges.length-1;
for(const j of charts) for(let a=0;a<n;a++)for(let b=0;b<n;b++)for(let c=0;c<n;c++)for(let d=0;d<n;d++)
  cells.push({id:[j,a,b,c,d].join('.'),chart:j,B:[a,b,c,d].map(i=>[edges[i],edges[i+1]])});
const done=new Set();
if(fs.existsSync(OUT)) for(const line of fs.readFileSync(OUT,'utf8').split('\n')) if(line.trim()){ try{ done.add(JSON.parse(line).id); }catch(e){} }
const todo=cells.filter(c=>!done.has(c.id));
console.log(`${cells.length} cells, ${done.size} done, ${todo.length} to run on ${JOBS} workers, budget ${BUDGET} boxes per cell`);
let next=0, running=0, closed=0, open=0; const t0=Date.now();
const enc=B=>B.map(x=>x.map(v=> v===Infinity?'inf': v===-Infinity?'-inf': v));
function start(){
  for(let j=0;j<JOBS&&j<todo.length;j++){
    const w=cp.fork(__filename,['--worker','1','--budget',String(BUDGET),'--c0',String(C0),'--minw',String(MINW)],{serialization:'advanced'});
    running++;
    const feed=()=>{ if(next<todo.length){ const c=todo[next++]; w.send({id:c.id,chart:c.chart,B:c.B}); } else { w.send('exit'); } };
    w.on('message',r=>{ r.B=enc(r.B); fs.appendFileSync(OUT,JSON.stringify(r)+'\n');
      if(r.status==='closed') closed++; else open++; feed(); });
    w.on('exit',()=>{ running--; if(running===0) finish(); });
    feed(); } }
const hb=setInterval(()=>{ const s=(Date.now()-t0)/1000;
  console.log(`  ${s.toFixed(0)}s  cells finished ${closed+open}/${todo.length}  closed ${closed}  open ${open}`); },15000);
function finish(){ clearInterval(hb); console.log(`finished: ${closed} closed, ${open} open, ${((Date.now()-t0)/1000).toFixed(0)} s`); }
if(todo.length) start(); else finish();
