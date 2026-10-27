from pathlib import Path
import subprocess,json,datetime,traceback
O=Path(__file__).resolve().parents[1]
def stamp():return datetime.datetime.now(datetime.timezone.utc).isoformat()
rc=1
try:
 with (O/'logs/pipeline.log').open('w',encoding='utf-8') as f:
  p=subprocess.Popen(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/pdac_spatial/bin/python','-u',str(O/'code/runner.py')],stdout=f,stderr=subprocess.STDOUT)
  (O/'runner.pid').write_text(str(p.pid)+'\n',encoding='utf-8')
  rc=p.wait()
except Exception:
 (O/'logs/supervisor_error.txt').write_text(traceback.format_exc(),encoding='utf-8')
(O/'logs/pipeline.exit_code').write_text(str(rc)+'\n',encoding='utf-8')
has=(O/'ACCEPTANCE.json').exists()
result={'exit_code':rc,'acceptance_present':has,'orchestration_success':rc==0 and has,'utc':stamp()}
if has:result['analysis_acceptance']=json.loads((O/'ACCEPTANCE.json').read_text(encoding='utf-8'))
(O/'JOB_FINISHED.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
if rc!=0 and not (O/'delivery/OVERNIGHT_SUMMARY_ZH.txt').exists():
 (O/'delivery/OVERNIGHT_SUMMARY_ZH.txt').write_text('后台任务未通过自动完成验收。\n退出码：'+str(rc)+'\n请查看 FAILED.json、STATUS.json、JOB_FINISHED.json 及 logs。已完成的独立结果保留；不得将进程结束当作分析成功。\n',encoding='utf-8')
pointer=O.parent/'CAF_EXPANSION_IDENTITY_CURRENT.json'
if pointer.exists():
 d=json.loads(pointer.read_text(encoding='utf-8'))
 if d.get('task_dir')==str(O):
  d.update(last_observed_utc=stamp(),background_finished=True,job_finished=result)
  pointer.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
