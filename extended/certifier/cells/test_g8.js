const T=require('./g8.js'), fs=require('fs'), {toEtaS}=require('./etaconv.js');
const KN=T.KNOWN, REF=(process.env.SCHED||'134569_r').startsWith('134569')?0:1; const toEta8=l=>toEtaS(l,KN,REF);
const box=(e,h)=>e.map(x=>[x-h,x+h]);
function run(file,opt,label,h=1e-12){ const pts=JSON.parse(fs.readFileSync(require('path').join(__dirname,file),'utf8')); const cnt={};
  for(const l of pts){ const r=T.evalBox(box(toEta8(l),h),opt); const k=r.status+' '+r.why; cnt[k]=(cnt[k]||0)+1; }
  console.log(label,JSON.stringify(cnt)); }
run('./v8_l.json',{noclauses:true},'V8 no clauses (all open):');
run('./v8_l.json',{},'V8 with clauses (all refuted by clauses):');
run('./lam_l.json',{},'laminar (all refuted by target):');
run('./lam_l.json',{notarget:true},'laminar without target (open, except U4U3 near T*):');
const E=JSON.parse(fs.readFileSync(require('path').join(__dirname,'v8_l.json'),'utf8')).map(toEta8);
const q=i=>{ const v=E.map(e=>e[i]).sort((a,b)=>a-b); return [v[0],v[v.length>>1],v[v.length-1]].map(x=>x.toFixed(2)).join(' '); };
console.log('V8 eta ranges (min med max):',[0,1,2,3].map(q).join(' | '));
const L=JSON.parse(fs.readFileSync(require('path').join(__dirname,'lam_l.json'),'utf8')).map(toEta8);
const q2=i=>{ const v=L.map(e=>e[i]).sort((a,b)=>a-b); return [v[0],v[v.length>>1],v[v.length-1]].map(x=>x.toFixed(2)).join(' '); };
console.log('laminar eta ranges:',[0,1,2,3].map(q2).join(' | '));
