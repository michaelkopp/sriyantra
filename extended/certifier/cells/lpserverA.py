"""JSON-lines server for chart A: {"id":..,"br":{"a,b":[lo,hi],...},"opt":{...}} -> {"id":..,"refuted":bool,"how":..}
   opt.handoff: also try to show that every point of the box on the relaxed V7 lies in the certified chart B region."""
import sys, json, os
sys.path.insert(0,__file__.rsplit('/',1)[0])
import bracketlp as B
BASEOPT=dict(ntan=int(os.environ.get('NTAN','1')),maxrows=os.environ.get('MAXROWS','0')=='1')
for line in sys.stdin:
    line=line.strip()
    if not line: continue
    m=json.loads(line); opt=dict(BASEOPT,**{k:v for k,v in m.get('opt',{}).items() if k!='handoff'})
    how=None
    try:
        bd=B.bounds_from_dump(m)
        if B.refute_bounds(bd,opt): how='LP'
        elif m.get('opt',{}).get('handoff') and B.handoff(bd,opt): how='handoff'
    except Exception as ex: sys.stderr.write('lp error %r\n'%(ex,))
    sys.stdout.write(json.dumps({'id':m['id'],'refuted':how is not None,'how':how})+'\n'); sys.stdout.flush()
