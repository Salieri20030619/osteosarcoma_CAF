from pathlib import Path
import json,datetime,shutil,pandas as pd
from range_zip import ROOT,TARGETS
qa=json.loads((ROOT/'DOWNLOAD_QA.json').read_text());assert qa['all_checks_passed']
s=json.loads((ROOT/'STATUS.json').read_text());assert s['status']=='DOWNLOAD_COMPLETE_VERIFIED'
t=pd.read_csv(ROOT/'sample_download_acceptance.tsv',sep='\t')
lines=['# ScPCA配对样本表达数据下载完成','',
 '本轮已通过官方API获取预生成AnnData数据包的授权链接，使用HTTP Range提取指定样本。表达文件下载、ZIP CRC32、SHA256记录、对象读取及官方元数据计数核对均已完成。','',
 '下载范围为SCPCP000017的5个样本和SCPCP000023的4个样本，对应4组患者骨肺样本关系。它们是否满足同步临床配对及CAF数量要求，需在后续分析中分别判断。','',
 '| 项目 | 样本 | 患者 | 部位 | filtered细胞 | processed细胞 |','|---|---|---|---|---:|---:|']
for _,r in t.iterrows():lines.append(f"| {r.scpca_project_id} | {r.scpca_sample_id} | {r.participant_id} | {r.tissue_location} | {r.filtered_cells} | {r.processed_cells} |")
lines+=['',f"共{qa['samples']}个样本、{qa['expression_objects']}个h5ad对象、{qa['selected_files']}个选定数据/报告文件；数据文件合计{s['total_bytes']:,}字节。每个样本均包含unfiltered、filtered、processed三个表达对象及QC/细胞类型报告。",'',
 '## 下载路线','',
 '最初下载进程没有代理环境变量，使用直连。对同一1 MiB片段的对照测试中，直连在约6秒后超时未收到数据；127.0.0.1:7897代理约4.3秒完成。随后下载程序显式设置该代理，采用8个传输线程、1 MiB分块，并以分块SHA256记录支持续传。测速为单次诊断，不能当作整个任务的稳定吞吐率。','',
 '使用项目已有登记邮箱取得官方API令牌；令牌及签名下载链接仅保存在权限受限文件中，不纳入交付包。未向作者发送申请邮件，未建立自定义打包任务。','',
 '## 验收范围','',
 '所有下载文件通过官方ZIP目录提供的CRC32及长度校验，并记录本地SHA256。27个AnnData对象可打开；每个对象抽查至多50行counts/X为有限非负数。所有processed条码包含于对应filtered对象，filtered条码包含于unfiltered对象；基因集合对应关系也通过核对。filtered和processed细胞数均与本轮官方元数据一致，九个样本均非PDX。','',
 '这是下载及文件完整性验收，尚未运行新的细胞质控、双细胞分析、CNV或CAF鉴定。filtered和processed细胞数不是CAF数量；SCPCS000794和SCPCS000812的processed数量较少，不能按旧表较大的filtered数量估计有效CAF样本量。','',
 '## 文件位置','',str(ROOT/'files'),'',
 'sample_download_acceptance.tsv记录样本验收；verified_files.json记录文件路径、长度、CRC32及SHA256；matrix_inventory.json记录对象结构；DOWNLOAD_QA.json记录验收结论。','',
 '官方说明：https://scpca.readthedocs.io/en/stable/download_files.html',
 '官方API示例：https://github.com/AlexsLemonade/scpca-portal/tree/dev/api-examples']
(ROOT/'下载完成记录.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
p=ROOT.parent;backup=p/'PROJECT_STATE.before_phase15_20260920.json'
if not backup.exists():shutil.copy2(p/'PROJECT_STATE.json',backup)
state=json.loads((p/'PROJECT_STATE.json').read_text());state.update(current_stage=ROOT.name,status='SCPCA_PAIRED_EXPRESSION_DOWNLOADED_VERIFIED',current_status='SCPCA_PAIRED_EXPRESSION_DOWNLOADED_VERIFIED',dynamic_status_file=str(ROOT/'STATUS.json'),latest_report=str(ROOT/'下载完成记录.md'),updated_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),phase15_scpca_download_20260920={'status':'DOWNLOAD_COMPLETE_VERIFIED','samples':9,'h5ad':27,'total_bytes':s['total_bytes'],'files_path':str(ROOT/'files'),'QC_CNV_started':False},next_action='Accept ScPCA paired source counts, timing and treatment evidence; then QC/CNV and customer myCAF panels; do not count processed cells as CAF or assert simultaneous pairs')
(p/'PROJECT_STATE.json').write_text(json.dumps(state,ensure_ascii=False,indent=2));print('PROJECT_DOWNLOAD_STATE_UPDATED')
