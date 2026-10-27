from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='4'
import subprocess,json,time,datetime,traceback,concurrent.futures,threading,csv,signal
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase88_expansion_identity_20261014'
PY='/mnt/BioDisk/wangmingyuan/miniconda3/envs/pdac_spatial/bin/python'
LOCK=threading.Lock();jobs=[]
def stamp():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def dump(p,x):
    tmp=p.with_suffix(p.suffix+'.tmp');tmp.write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');tmp.replace(p)
def mem():
    for line in Path('/proc/meminfo').read_text(encoding='ascii').splitlines():
        if line.startswith('MemAvailable:'):return int(line.split()[1])/1024**2
def run(cmd,log,timeout=16*3600):
    with log.open('w',encoding='utf-8') as f:
        p=subprocess.Popen(cmd,stdout=f,stderr=subprocess.STDOUT,start_new_session=True)
        try:return p.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            os.killpg(p.pid,signal.SIGTERM)
            try:p.wait(timeout=30)
            except subprocess.TimeoutExpired:os.killpg(p.pid,signal.SIGKILL);p.wait()
            return 124
def lane(method,rows):
    result=[]
    for r in rows:
        sm=r['sample'];eligible=r['ambient_eligible']=='True' if method=='decontx' else r['cnv_eligible']=='True'
        row={'method':method,'sample':sm,'started_utc':stamp()}
        if not eligible:
            row.update(status='skipped_input_gate',reason='fewer_than_two_lineages' if method=='decontx' else r['reason'])
        else:
            while mem()<96:
                dump(O/'resources'/('MEMORY_WAIT_'+method+'.json'),{'utc':stamp(),'sample':sm,'available_GiB':mem(),'required_GiB':96,'check_interval_seconds':900})
                time.sleep(900)
            row.update(status='running')
            with LOCK:jobs.append(row.copy());dump(O/'JOB_PROGRESS.json',jobs)
            log=O/'logs'/(method+'__'+sm+'.log')
            rc=run(['bash',str(O/'code/run_tool.sh'),method,sm],log)
            (O/'logs'/(method+'__'+sm+'.exit_code')).write_text(str(rc)+'\n',encoding='utf-8')
            complete=(O/method/sm/'COMPLETE').exists()
            row.update(exit_code=rc,completion_marker=complete,status='complete' if rc==0 and complete else 'failed',ended_utc=stamp())
            print('TOOL_FINISHED',json.dumps(row),flush=True)
        result.append(row)
        with LOCK:
            jobs[:]=[j for j in jobs if (j['method'],j['sample'])!=(method,sm)];jobs.append(row);dump(O/'JOB_PROGRESS.json',jobs)
    return result
def main():
    dump(O/'STATUS.json',{'status':'preparing_and_reviewing_lineages','utc':stamp(),'pid':os.getpid(),'completion_file':str(O/'ACCEPTANCE.json')})
    for mode in ['cnv','ambient']:assert json.loads((O/'resources'/('PREFLIGHT_'+mode+'.json')).read_text(encoding='utf-8'))['status']=='PASS'
    rc=run([PY,'-u',str(O/'code/prepare.py')],O/'logs/prepare.log',4*3600)
    (O/'logs/prepare.exit_code').write_text(str(rc)+'\n',encoding='utf-8')
    assert rc==0 and (O/'INPUTS_READY').exists(),'Preparation failed; see logs/prepare.log'
    with (O/'review/input_summary.tsv').open(encoding='utf-8') as f:rows=list(csv.DictReader(f,delimiter='\t'))
    assert len(rows)==12
    rows.sort(key=lambda r:(-int(r['FibLike']),-int(r['MesAmb'])))
    dump(O/'STATUS.json',{'status':'CNV_and_ambient_running','utc':stamp(),'pid':os.getpid(),'libraries':12,'within_source_patient_or_donor_labels':10,'CNV_eligible_samples':sum(r['cnv_eligible']=='True' for r in rows),'method_lanes':['infercnv','copykat','decontx'],'threads_per_lane':4,'baseline_changed':False})
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:
        futures=[ex.submit(lane,m,rows) for m in ['infercnv','copykat','decontx']]
        alljobs=[r for f in futures for r in f.result()]
    dump(O/'TOOL_RESULTS.json',alljobs)
    dump(O/'STATUS.json',{'status':'integrating_and_validating','utc':stamp(),'pid':os.getpid()})
    rc=run([PY,'-u',str(O/'code/finalize.py')],O/'logs/finalize.log',4*3600)
    (O/'logs/finalize.exit_code').write_text(str(rc)+'\n',encoding='utf-8')
    assert rc==0 and (O/'ACCEPTANCE.json').exists(),'Integration/acceptance failed; see logs/finalize.log'
    a=json.loads((O/'ACCEPTANCE.json').read_text(encoding='utf-8'))
    dump(O/'STATUS.json',{'status':a['status'],'utc':stamp(),'report':str(O/'delivery/OVERNIGHT_SUMMARY_ZH.txt'),'acceptance':str(O/'ACCEPTANCE.json'),'baseline_changed':False,'manual_review_pending':True})
    print('OVERNIGHT_COMPUTATION_FINISHED',json.dumps(a),flush=True)
if __name__=='__main__':
    try:main()
    except Exception:
        err={'status':'failed','utc':stamp(),'traceback':traceback.format_exc()}
        dump(O/'FAILED.json',err);dump(O/'STATUS.json',err);traceback.print_exc();raise
