from pathlib import Path
import os,sys,json,datetime,fcntl,select,subprocess,traceback,csv,hashlib
BASE=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase26_client_reanalysis_20260922')
OUT=BASE/'retry_once_20260924'
PY='/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python'
R='/mnt/BioDisk/wangmingyuan/project26210_envs/phase11_infercnv_conda_v2/bin/Rscript'
def now(): return datetime.datetime.now().astimezone().isoformat()
def read(p): return json.loads(p.read_text())
def write(p,x):
    t=p.with_name(p.name+'.tmp');t.write_text(json.dumps(x,ensure_ascii=False,indent=2));t.replace(p)
def success(j):
    p=Path(j['path'])
    return j.get('exit_code')==0 and str(j.get('status','')).startswith('complete') and (p/'COMPLETE').is_file() and (p/'cell_results.tsv').is_file() and (p/'cell_results.tsv').stat().st_size>0
def targets(rows):
    return [(r,j) for r in rows for j in r.get('jobs',[]) if j.get('method') in ('copykat','infercnv') and not success(j)]
def pid_fields(pid): return (Path('/proc')/str(pid)/'stat').read_text().rsplit(')',1)[1].split()
def pidfd(pid):
    if hasattr(os,'pidfd_open'): return os.pidfd_open(pid)
    import ctypes
    lib=ctypes.CDLL(None,use_errno=True);fd=lib.syscall(434,int(pid),0)
    if fd<0: raise OSError(ctypes.get_errno(),'pidfd_open')
    return fd
def validate_cells(p,sm):
    with (BASE/'inputs'/sm/'cells.tsv').open() as f: expected=[r['cell_id'] for r in csv.DictReader(f,delimiter='\t')]
    with p.open() as f:
        reader=csv.DictReader(f,delimiter='\t')
        assert {'cell_id','call'} <= set(reader.fieldnames or []),'missing result columns'
        rows=list(reader)
    ids=[r['cell_id'] for r in rows]
    assert len(ids)==len(set(ids))==len(expected) and set(ids)==set(expected),'cell ID integrity failure'
    assert all(r['call'] for r in rows),'empty calls'
    return len(rows)
def wait_first_pass(plan):
    try:
        fields=pid_fields(plan['runner_pid'])
        assert fields[19]==plan['runner_start_ticks'],'runner identity changed'
        fd=pidfd(plan['runner_pid'])
        try: select.select([fd],[],[])
        finally: os.close(fd)
    except FileNotFoundError: pass
    locks=[]
    for name in ['runner_v6.lock','runner.lock']:
        h=(BASE/name).open('a');fcntl.flock(h,fcntl.LOCK_EX|fcntl.LOCK_NB);locks.append(h)
    assert (BASE/'RUN_FINISHED').exists() and (BASE/'CNV_STAGE_FINISHED').exists(),'First pass did not finish; inspect manually'
    first_rc=int((BASE/'exit_code').read_text().strip());qa=read(BASE/'delivery/CNV_STAGE_QA.json')
    rows=read(BASE/'tables/job_outcomes.json')
    assert len(rows)==107 and len({r['sample'] for r in rows})==107 and qa['recorded_outcomes']==107
    assert first_rc in (0,1) and read(BASE/'STATUS.json')['status']!='RUNNER_FAILED'
    write(OUT/'FIRST_PASS_SNAPSHOT.json',{'at':now(),'exit_code':first_rc,'QA':qa,'outcomes':rows})
    return rows,first_rc,locks

def reservation_for_pending(pending,budgets,allocations,pss_fn):
    reserved=0.0
    for sm in pending:
        used=sum(pss_fn(v['pid']) for key,v in allocations.items() if key.split(':',1)[-1]==sm and v.get('pid'))
        reserved+=max(0,budgets[sm]-used)
    return reserved

def early_memory_available():
    from memory_gate import available_gib,tree_pss,estimate_gib
    raw=available_gib()
    rows=read(BASE/'tables/job_outcomes.json')
    done={r['sample'] for r in rows}
    with (BASE/'tables/analysis_source_manifest.tsv').open() as f:
        pending={r['sample'] for r in csv.DictReader(f,delimiter='\t')}-done
    allocations=read(BASE/'MEMORY_GUARD_STATE.json').get('allocations',{})
    budgets={sm:estimate_gib(read(BASE/'inputs'/sm/'input_QA.json')) for sm in pending}
    reserve=reservation_for_pending(pending,budgets,allocations,tree_pss)
    write(OUT/'FIRST_PASS_MEMORY_RESERVATION.json',{'at':now(),'raw_available_GiB':raw,'first_pass_future_reserved_GiB':reserve,'effective_available_for_retry_GiB':raw-reserve,'pending_samples':sorted(pending),'pending_budgets_GiB':budgets,'includes_not_yet_started_first_pass_samples':True})
    return raw-reserve

def main():
    lock=(OUT/'retry.lock').open('a')
    try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError:return
    if (OUT/'exit_code').exists():return
    plan=read(OUT/'PLAN.json')
    state={'status':'WAITING_FOR_FIRST_PASS','at':now(),'pid':os.getpid(),'first_pass_pid':plan['runner_pid'],'retry_limit_per_failed_method':1}
    write(OUT/'STATUS.json',state);(OUT/'runner.pid').write_text(str(os.getpid())+'\n')
    rows=read(BASE/'tables/job_outcomes.json')
    assert len(rows)==len({r['sample'] for r in rows})
    write(OUT/'EARLY_FIRST_PASS_SNAPSHOT.json',{'at':now(),'outcomes':rows,'first_pass_complete':False})
    first_rc=None
    locks=[]
    todo=targets(rows)
    if not todo:
        rows,first_rc,locks=wait_first_pass(plan)
        todo=targets(rows)
    write(OUT/'retry_manifest.json',[{'sample':r['sample'],'method':j['method'],'first_pass':j} for r,j in todo])
    sys.path.insert(0,str(BASE/'code'))
    from memory_gate import MemoryGate,estimate_gib
    gate=MemoryGate(OUT/'MEMORY_GUARD_STATE.json',reserve=50,available_fn=early_memory_available)
    # Reuse the existing evidence aggregation function, directing new evidence to this version.
    ns={};src=(BASE/'code/run_background.py').read_text().split('\nlock=')[0]
    exec(compile(src,str(BASE/'code/run_background.py'),'exec'),ns);ns['O']=OUT
    manifest={r['sample']:r for r in ns['read'](BASE/'tables/analysis_source_manifest.tsv').fillna('').to_dict('records')}
    original=rows
    updated={r['sample']:r for r in json.loads(json.dumps(rows))}
    outcomes=read(OUT/'retry_outcomes.json') if (OUT/'retry_outcomes.json').exists() else []
    completed_keys={(x['sample'],x['method']) for x in outcomes}
    assert len(completed_keys)==len(outcomes),'Duplicate completed retry outcomes'
    affected={x['sample'] for x in outcomes}
    for x in outcomes:
        result=x['retry'];d=Path(result['path'])
        assert (d/'ATTEMPT_STARTED').exists() and int((d/'exit_code').read_text())==result['exit_code']
        rr=updated[x['sample']]
        rr['jobs']=[result if j['method']==x['method'] else j for j in rr['jobs']]
        rr['status']='CNV_EVIDENCE_READY' if all(success(j) for j in rr['jobs']) else 'CNV_EVIDENCE_WITH_SOFTWARE_FAILURE'
    scheduled={(r['sample'],j['method']) for r,j in todo}|completed_keys
    todo=[(r,j) for r,j in todo if (r['sample'],j['method']) not in completed_keys]
    state.update(status='RETRYING',scheduled_methods=len(scheduled),finished_methods=len(outcomes),first_pass_exit_code=first_rc,address_space_limit_GiB=None)
    write(OUT/'STATUS.json',state)
    for row,old in todo:
        sm=row['sample'];method=old['method'];assert (sm,method)==('SCPCS000525','copykat'),'Override is limited to SCPCS000525 CopyKAT';affected.add(sm)
        d=OUT/'cnv'/method/sm
        if d.exists():assert not any(d.iterdir()),'Refuse to repeat an attempted retry'
        else:d.mkdir(parents=True,exist_ok=False)
        budget=estimate_gib(read(BASE/'inputs'/sm/'input_QA.json'))
        key=method+':'+sm
        state.update(active_sample=sm,active_method=method,at=now());write(OUT/'STATUS.json',state)
        state.update(status='WAITING_FOR_SHARED_MEMORY');write(OUT/'STATUS.json',state)
        gate.admit(key,budget)
        rc=None;validation_error='';valid_cells=None;resource_abort=''
        (d/'ATTEMPT_STARTED').write_text(now()+'\n')
        try:
            env=os.environ.copy()
            for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:env[k]='1'
            env['CAF_RETRY_ROOT']=str(OUT)
            with (OUT/'logs'/f'{method}_{sm}.log').open('x') as log:
                p=subprocess.Popen([R,'--vanilla',str(OUT/'code/native_cnv_retry.R'),method,sm],stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,env=env,start_new_session=True)
                gate.attach(key,p.pid)
                state.update(status='RETRYING',native_pid=p.pid,at=now());write(OUT/'STATUS.json',state)
                fd=pidfd(p.pid)
                try:
                    while not select.select([fd],[],[],30)[0]:
                        from memory_gate import available_gib
                        available=available_gib()
                        write(OUT/'ACTIVE_SAMPLE_MEMORY_GUARD.json',{'at':now(),'sample':sm,'native_pid':p.pid,'available_GiB':available,'minimum_available_GiB':50,'action':'terminate_current_sample' if available<50 else 'continue'})
                        if available<50:
                            import signal
                            resource_abort='Stopped by current-sample resource guard: available memory fell below 50 GiB'
                            os.killpg(p.pid,signal.SIGTERM)
                            if not select.select([fd],[],[],15)[0]:os.killpg(p.pid,signal.SIGKILL)
                            break
                        with gate.cond:gate.status()
                        state.update(at=now(),native_pid=p.pid);write(OUT/'STATUS.json',state)
                finally:os.close(fd)
                rc=p.wait()
            (d/'exit_code').write_text(str(rc)+'\n')
            ok=rc==0 and (d/'COMPLETE').exists() and (d/'cell_results.tsv').exists()
            if ok:
                try:valid_cells=validate_cells(d/'cell_results.tsv',sm)
                except Exception as e:ok=False;validation_error=str(e)
        finally:gate.release(key)
        failure=resource_abort or ((d/'FAILED').read_text().strip() if (d/'FAILED').exists() else validation_error)
        result={'method':method,'status':'complete' if ok else 'software_failed','exit_code':rc,'path':str(d),'failure':failure,'validated_cells':valid_cells,'native_threads':8,'blas_threads_per_process':1,'retry_attempt':1,'address_space_limit_GiB':None,'resource_interrupted_predecessor':str(OUT/'resource_restart_20260925')}
        outcomes.append({'sample':sm,'method':method,'first_failure':old.get('failure',''),'retry':result,'interpretation':'recovered_on_same_parameter_retry' if ok else 'failure_reproduced' if failure==old.get('failure') else 'failed_again_different_or_unrecorded_error'})
        rr=updated[sm];rr['jobs']=[result if j['method']==method else j for j in rr['jobs']]
        rr['status']='CNV_EVIDENCE_READY' if all(success(j) for j in rr['jobs']) else 'CNV_EVIDENCE_WITH_SOFTWARE_FAILURE'
        rr['retry_finished_at']=now()
        write(OUT/'retry_outcomes.json',outcomes)
        state.update(finished_methods=len(outcomes),native_pid=None,at=now());write(OUT/'STATUS.json',state)
        if len(outcomes)==len(scheduled):
            state.update(status='WAITING_FOR_FIRST_PASS_TAIL',active_sample=None,active_method=None,at=now());write(OUT/'STATUS.json',state)
            if not locks:all_rows,first_rc,locks=wait_first_pass(plan)
            original=all_rows
            for r in all_rows:
                if r['sample'] not in affected:updated[r['sample']]=r
            extra=[(r,j) for r,j in targets(all_rows) if (r['sample'],j['method']) not in scheduled]
            for r,j in extra:scheduled.add((r['sample'],j['method']))
            todo.extend(extra)
            write(OUT/'retry_manifest.json',[{'sample':r['sample'],'method':j['method'],'first_pass':j} for r,j in todo])
            state.update(scheduled_methods=len(scheduled),first_pass_exit_code=first_rc,at=now());write(OUT/'STATUS.json',state)
    for sm in affected:
        updated[sm]['evidence_counts']=ns['combine'](manifest[sm],updated[sm]['jobs'])
    write(OUT/'tables/job_outcomes_after_retry.json',list(updated.values()))
    unresolved=[r['sample'] for r in updated.values() if r.get('status')!='CNV_EVIDENCE_READY']
    final={'at':now(),'first_pass_outcomes':len(original),'retried_methods':len(outcomes),'recovered_methods':sum(x['retry']['status']=='complete' for x in outcomes),'failed_again_methods':sum(x['retry']['status']!='complete' for x in outcomes),'unresolved_samples':unresolved,'all_retry_process_exit_codes_recorded':all(x['retry']['exit_code'] is not None for x in outcomes),'formal_CAF_clustering_started':False,'original_results_preserved':True,'same_seed_and_method_parameters':True,'current_retry_address_space_limit_GiB':None,'resource_restart':'SCPCS000525 incomplete capped attempt preserved under resource_restart_20260925; completed method retries not repeated','interpretation':'Recovery is compatible with a transient execution issue; repeat failure is not automatically a biological exclusion.','evidence_location':'Retried sample evidence is in this tables directory; unchanged samples use first-pass evidence.'}
    write(OUT/'delivery/RETRY_QA.json',final)
    (OUT/'RETRY_STAGE_FINISHED').write_text('ALL_SCHEDULED_METHODS_RETRIED_ONCE\n')
    state.update(status='RETRY_COMPLETE_REVIEW_REQUIRED',at=now(),active_sample=None,active_method=None);write(OUT/'STATUS.json',state)
    (OUT/'exit_code').write_text('0\n')
if __name__=='__main__':
    try: main()
    except Exception as e:
        traceback.print_exc()
        write(OUT/'ERROR.json',{'at':now(),'error':str(e)})
        write(OUT/'STATUS.json',{'status':'RETRY_CONTROLLER_FAILED','at':now(),'error':str(e)})
        (OUT/'exit_code').write_text('1\n')
        sys.exit(1)
