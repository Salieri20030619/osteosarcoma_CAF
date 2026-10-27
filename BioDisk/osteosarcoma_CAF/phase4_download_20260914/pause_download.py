"""SYNTAX_CHECK: python -m py_compile pause_download.py
RUN_COMMAND: python pause_download.py --root /mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase4_download_20260914
EXPECTED_SUCCESS_MARKER: DOWNLOAD_PAUSED_AND_SAVED
"""
from pathlib import Path
import os,signal,time,json,argparse,zipfile
from datetime import datetime,timezone
ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);a=ap.parse_args();root=a.root.resolve()
def processes():
 out={}
 for p in Path('/proc').iterdir():
  if not p.name.isdigit():continue
  try:
   cmd=(p/'cmdline').read_bytes().replace(b'\0',b' ').decode(errors='replace')
   st=(p/'stat').read_text();rest=st[st.rfind(')')+2:].split()
   out[int(p.name)]={'pid':int(p.name),'ppid':int(rest[1]),'state':rest[0],'cmd':cmd,'stat':st,'cwd':str((p/'cwd').resolve())}
  except (OSError,ValueError):pass
 return out
before=processes();ancestors={os.getpid()};pid=os.getppid()
while pid in before and pid not in ancestors:ancestors.add(pid);pid=before[pid]['ppid']
names=['run_http1.sh','download_validate_http1.py','wait_http1.sh','run_geo.sh','download_validate_geo.py','wait_completion.sh','finalize_download_review.py','run_download.sh','download_validate.py']
targets={pid for pid,r in before.items() if pid not in ancestors and (str(root) in r['cmd'] or r['cwd']==str(root)) and any(n in r['cmd'] for n in names)}
while True:
 extra={pid for pid,r in before.items() if r['ppid'] in targets and pid not in ancestors}
 if extra<=targets:break
 targets|=extra
for pid in targets:
 try:os.kill(pid,signal.SIGSTOP)
 except ProcessLookupError:pass
(root/'pause_process_snapshot.json').write_text(json.dumps([before[p] for p in sorted(targets)],indent=2))
for pid in targets:
 try:os.kill(pid,signal.SIGTERM);os.kill(pid,signal.SIGCONT)
 except ProcessLookupError:pass
time.sleep(2)
after=processes()
for pid in targets:
 if pid in after and after[pid]['state']!='Z':
  try:os.kill(pid,signal.SIGKILL)
  except ProcessLookupError:pass
time.sleep(1)
after=processes();live=[pid for pid in targets if pid in after and after[pid]['state']!='Z']
assert not live,live
manifest=json.loads((root/'download_manifest.json').read_text());files=[]
for s in manifest:
 for url in s['urls']:
  f=root/'downloads'/s['dataset']/url.rsplit('/',1)[-1];meta=f.with_name(f.name+'.verified.json')
  r={'dataset':s['dataset'],'sample':s['sample'],'file':f.name,'verified_complete':f.exists() and meta.exists(),'completed_bytes':f.stat().st_size if f.exists() else 0,'partials':[]}
  if meta.exists():r['verified_record']=json.loads(meta.read_text())
  for part in f.parent.glob(f.name+'*.part'):r['partials'].append({'path':str(part),'bytes':part.stat().st_size})
  files.append(r)
summary={'status':'PAUSED_BY_USER','timestamp_utc':datetime.now(timezone.utc).isoformat(),'files_total':len(files),'files_verified':sum(r['verified_complete'] for r in files),'verified_bytes':sum(r['completed_bytes'] for r in files if r['verified_complete']),'stopped_pids':sorted(targets),'remaining_live_target_pids':live,'matrix_validation_complete':(root/'COMPLETE').exists(),'no_new_cohort_QC_or_CNV_started':True,'files':files}
(root/'pause_state.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
(root/'PAUSED_BY_USER').write_text(summary['timestamp_utc']+'\nDownloads, retries, completion waiting and automatic finalization stopped.\n')
report=f'''# 下载工作已保存并暂停

暂停时间（UTC）：{summary['timestamp_utc']}。

已完整下载并通过压缩读取校验：{summary['files_verified']}/{summary['files_total']}个文件，{summary['verified_bytes']:,}字节。
其余文件和部分下载全部保留在downloads/。逐文件状态、部分文件大小、已完成文件SHA256见pause_state.json。

本项目下载、重试、完成等待和自动打包链均已停止；停止进程及资源快照见pause_process_snapshot.json。新队列QC、聚类、CAF注释和CNV未启动。全批次矩阵读取校验尚未完成，不得把本次暂停记作校验成功。

恢复时先读取本说明、pause_state.json及download_http1.log。已验证完整文件可复用；部分.geo.part文件不应当作完整输入。当前GEO备用入口的脚本对失败大文件整文件重试，不保证部分文件可直接续传；下一次先解决反复断连和续传策略，再启动剩余下载。不要直接运行旧wait脚本或假定COMPLETE存在。

非PDX来源审核已完成：GSE299023六个患者组织及GSE152048十一个患者组织。来源与临床原始值分别保存在download_manifest.json、clinical_source_audit.tsv及论文证据文件中；通过来源审核不等于已确认CAF。
'''
(root/'暂停记录.md').write_text(report,encoding='utf-8')
include=['暂停记录.md','pause_state.json','pause_process_snapshot.json','PAUSED_BY_USER','download_manifest.json','clinical_source_audit.tsv','download_status.json','download_http1.log','download_geo.log','download_validate.log','resource_usage_http1.txt','resource_usage_geo.txt','resource_usage.txt','运行与监控.md']
with zipfile.ZipFile(root/'download_pause_record.zip','w',zipfile.ZIP_DEFLATED) as z:
 for name in include:
  if (root/name).exists():z.write(root/name,name)
 for f in root.glob('*.py'):z.write(f,f.name)
 for f in root.glob('*.sh'):z.write(f,f.name)
print(json.dumps({k:v for k,v in summary.items() if k!='files'},indent=2))
print('DOWNLOAD_PAUSED_AND_SAVED')
