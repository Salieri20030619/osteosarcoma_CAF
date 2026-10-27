import pathlib,subprocess,os,json,datetime,sys
p=pathlib.Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase7_review_states_20260917')
root=p/'download_expansion_20260917'
os.chdir(p)
os.environ.update(http_proxy='http://127.0.0.1:7897',https_proxy='http://127.0.0.1:7897',HTTP_PROXY='http://127.0.0.1:7897',HTTPS_PROXY='http://127.0.0.1:7897')
with (root/'logs/geo_runner_proxy.log').open('a') as log:
 rc=subprocess.run([sys.executable,'-u',str(root/'code/download_geo.py')],stdout=log,stderr=subprocess.STDOUT).returncode
(root/'logs/geo_runner_proxy.exit_code').write_text(str(rc)+'\n')
f=root/'geo_download_completion.json'
completion=json.loads(f.read_text()) if f.exists() else {}
ok=rc==0 and completion.get('total_files')==44 and completion.get('completed_files')==44
s=json.loads((root/'STATUS.json').read_text())
s['state']='PAUSED_AFTER_DOWNLOAD_ATTEMPTS'
s['GEO_completion']={'process_exit_code':rc,'all_files_complete':ok,'total_files':completion.get('total_files'),'completed_files':completion.get('completed_files'),'completion_file':str(f) if f.exists() else None}
s['finished_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
(root/'STATUS.json').write_text(json.dumps(s,indent=2))
z=json.loads((p.parent/'PROJECT_STATE.json').read_text())
z['current_status']=s['state'];z['download_expansion']=s
(p.parent/'PROJECT_STATE.json').write_text(json.dumps(z,ensure_ascii=False,indent=2))
(root/'logs/geo_proxy_completion.json').write_text(json.dumps(s['GEO_completion'],indent=2))
print('DOWNLOAD_ATTEMPT_FINISHED',s['GEO_completion'],flush=True)
sys.exit(0 if ok else 2)
