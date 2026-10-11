'use strict';
/* Chart A over a grid of cells, with LP certificates (bbZLP.js, lpserverA.py).  Each cell is closed by branch and bound with a budget; the result of
   every cell is appended to results.jsonl, and cells already there are skipped when the run is resumed.
     node runcells.js --jobs 10 --budget 2e5 [--charts 0,1,2,3,4] [--edges -inf,-8,...,0] [--out results2.jsonl]
     node runcells.js --jobs 10 --budget 5e6 --from results2.jsonl --out results3.jsonl   (open cells only)
   A cell is 'closed' when every box in it is refuted: by the interval tests of g8.js, by an LP certificate
   (bracketlp.py), or by an LP certificate that every point of the box on the relaxed V7 lies in the
   region certified by chart B ('handoff').  With --c0 0 nothing is delegated by u(F2).
     node runcellsLP.js --jobs 10 --budget 2e5 --c0 0 --from results2.jsonl --out results3.jsonl
   (--from takes the cells that FILE records as open, and the closed cells that delegated boxes.)
     --box 'CHART:[[a,b],...]' runs one explicit cell; --sub k splits the selected cells; --dumpopen FILE appends
     the boxes left open (one JSON line each: cell id, B, zc, depth, reason); --redo 1 runs again the cells
     that OUT records as open (the last record of a cell counts). */
const fs=require('fs'), path=require('path'), cp=require('child_process'), os=require('os');
const args=process.argv.slice(2), opt={};
for(let i=0;i<args.length;i++) if(args[i].startsWith('--')) opt[args[i].slice(2)]=args[i+1];
const OUT=opt.out||'results3.jsonl', JOBS=+(opt.jobs||Math.max(1,os.cpus().length-1)), BUDGET=+(opt.budget||2e5);
const C0=+(opt.c0||0.05), MINW=+(opt.minw||1e-3);
if(opt.worker){                                   // worker: read cells from stdin lines, write results
  process.env.SCHED=process.env.SCHED||'134569_r';
  const bb=require('./bbZLP.js');
  process.on('message',async m=>{
    if(m==='exit') process.exit(0);
    const t=Date.now();
    const o=await bb.run(m.B,{c0:C0,minW:MINW,maxBoxes:BUDGET,chart:m.chart});
    const openWhy={}; for(const u of o.open) openWhy[u.why]=(openWhy[u.why]||0)+1;
    process.send({id:m.id,chart:m.chart,B:m.B,status:o.open.length?'open':'closed',boxes:o.boxes,refuted:o.refuted,
      delegated:o.delegated,open:o.open.length,openWhy,mass:o.mass,why:o.why,lp:[o.lpRefuted,o.lpCalls],seconds:(Date.now()-t)/1000,
      openBoxes: m.dump? o.open.map(u=>({B:u.B,zc:u.zc||null,d:u.d,why:u.why})) : undefined});
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
// --from FILE: run again only the cells that FILE records as open (typically with a larger budget)
let only=null;
if(opt.from){ only=new Set(); for(const line of fs.readFileSync(opt.from,'utf8').split('\n')) if(line.trim()){ const r=JSON.parse(line); if(r.status!=='closed'||(r.delegated>0&&!opt.keepdelegated)) only.add(r.id); } }
if(opt.ids){ only=new Set(opt.ids.split(',')); }
// --box 'CHART:[[a,b],[a,b],[a,b],[a,b]]': one explicit cell (id --boxid, default 'box')
if(opt.box){ const k=opt.box.indexOf(':'); cells.length=0;
  cells.push({id:opt.boxid||'box',chart:+opt.box.slice(0,k),B:JSON.parse(opt.box.slice(k+1).replace(/-inf/g,'-1e999'))}); }
if(only&&!opt.box) for(let i=cells.length-1;i>=0;i--) if(!only.has(cells[i].id)) cells.splice(i,1);
// --sub k: split every selected cell into k^4 sub-cells (ids CELL/i.j.k.l)
if(opt.sub){ const k=+opt.sub, out=[];
  const parts=x=>{ const [a,b]=x, r=[];
    if(Number.isFinite(a)) for(let i=0;i<k;i++) r.push([a+(b-a)*i/k, i===k-1? b : a+(b-a)*(i+1)/k]);
    else { r.push([a,b-(k-1)]); for(let i=k-1;i>=1;i--) r.push([b-i,b-i+1]); }
    return r; };
  for(const c of cells){ const P=c.B.map(parts);
    for(let i=0;i<k;i++)for(let j=0;j<k;j++)for(let l=0;l<k;l++)for(let m=0;m<k;m++)
      out.push({id:c.id+'/'+[i,j,l,m].join('.'),chart:c.chart,B:[P[0][i],P[1][j],P[2][l],P[3][m]]}); }
  cells.length=0; cells.push(...out); }
const done=new Set();
// cells already in OUT are skipped (with --redo 1: only the closed ones, so open cells run again)
if(fs.existsSync(OUT)) for(const line of fs.readFileSync(OUT,'utf8').split('\n')) if(line.trim()){ try{ const r=JSON.parse(line);
  if(opt.redo&&r.status!=='closed') done.delete(r.id); else done.add(r.id); }catch(e){} }
const todo=cells.filter(c=>!done.has(c.id));
// the LP server needs a python3 with numpy and scipy; the workers inherit PYTHON
{ const cands=[process.env.PYTHON,'python3','/opt/homebrew/bin/python3','/usr/local/bin/python3',
    '/Library/Frameworks/Python.framework/Versions/Current/bin/python3','/usr/bin/python3'].filter(Boolean);
  let ok=null; for(const c of cands){ const t=cp.spawnSync(c,['-c','import numpy, scipy.optimize'],{stdio:'ignore'}); if(t.status===0){ ok=c; break; } }
  if(!ok){ console.error('No python3 with numpy and scipy found. Run with PYTHON=/path/to/python3 (the one where python3 test_lp.py passes).'); process.exit(2); }
  process.env.PYTHON=ok; console.log('LP server python: '+ok); }
console.log(`${cells.length} cells, ${done.size} done, ${todo.length} to run on ${JOBS} workers, budget ${BUDGET} boxes per cell`);
let next=0, running=0, closed=0, open=0; const t0=Date.now();
const enc=B=>B.map(x=>x.map(v=> v===Infinity?'inf': v===-Infinity?'-inf': v));
function start(){
  for(let j=0;j<JOBS&&j<todo.length;j++){
    const w=cp.fork(__filename,['--worker','1','--budget',String(BUDGET),'--c0',String(C0),'--minw',String(MINW)],{serialization:'advanced'});
    running++;
    const feed=()=>{ if(next<todo.length){ const c=todo[next++]; w.send({id:c.id,chart:c.chart,B:c.B,dump:!!opt.dumpopen}); } else { w.send('exit'); } };
    w.on('message',r=>{ r.B=enc(r.B);
      if(opt.dumpopen&&r.openBoxes){ const lines=r.openBoxes.map(u=>JSON.stringify({id:r.id,chart:r.chart,B:enc(u.B),zc:u.zc?u.zc.map(z=>z?enc([z])[0]:null):null,d:u.d,why:u.why})+'\n');
        if(lines.length) fs.appendFileSync(opt.dumpopen,lines.join('')); }
      delete r.openBoxes; fs.appendFileSync(OUT,JSON.stringify(r)+'\n');
      if(r.status==='closed') closed++; else open++; feed(); });
    w.on('exit',()=>{ running--; if(running===0) finish(); });
    feed(); } }
const hb=setInterval(()=>{ const s=(Date.now()-t0)/1000;
  console.log(`  ${s.toFixed(0)}s  cells finished ${closed+open}/${todo.length}  closed ${closed}  open ${open}`); },15000);
function finish(){ clearInterval(hb); console.log(`finished: ${closed} closed, ${open} open, ${((Date.now()-t0)/1000).toFixed(0)} s`); }
if(todo.length) start(); else finish();
