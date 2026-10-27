from pathlib import Path
import zipfile,hashlib,json,datetime
B=Path('C:/Users/31937/CAF_revision_20261002');D=Path('C:/Users/31937/Downloads')
out=D/'CAF修订与模型评估证据包_20261002.zip'
items={}
exclude={'REVISION_INTERIM_FACTS.json','summarize_interim.py','download_verified.mjs','ssh_bridge.py','early_revision_inputs.zip','final_revision_inputs.zip','surface_first_results.txt','ARTIFACT_OPERATION_MARKED.json','build_revised_package.py'}
for p in B.iterdir():
 if p.is_file() and p.name not in exclude and p.suffix in ['.py','.R','.json','.txt','.tsv','.diff','.log']:items['revision/'+p.name]=p
for directory in ['final_QA','inputs']:
 for p in (B/directory).rglob('*'):
  if p.is_file() and p.name!='PATHWAY_PROGRESS.json':items['revision/'+str(p.relative_to(B)).replace('\\','/')]=p
for name in ['CAF全流程修订与模型工作评估报告_20261002.pdf','CAF重启后审计证据包_20261002.zip','CAF重启后流程与结果审计_20261002.txt']:
 p=D/name;assert p.is_file(),name;items['reports_and_original_audit/'+name]=p
manifest=[]
readme="""本证据包供评价本次6.1sol实施与astra审核及修订。先读评估PDF，再读REVISION_CLOSURE.json及FINAL_REVISION_FACTS.json。original audit ZIP为修订前审计；revision/inputs为修订复算与比较证据，final_QA为最终验收。模型用量为报告生成时主会话快照，非账单。报告和统计大包另交付；统计大包SHA256见PACKAGE_INTEGRITY.json。未加入完整私有会话文本、凭据或隐藏推理。原服务器审计与阶段49目录继续保留。
"""
with zipfile.ZipFile(out,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 z.writestr('README_阅读顺序.txt',readme.encode('utf-8'))
 for name,p in sorted(items.items()):
  data=p.read_bytes();z.writestr(name,data);manifest.append({'path':name,'bytes':len(data),'SHA256':hashlib.sha256(data).hexdigest()})
 z.writestr('EVIDENCE_MANIFEST.json',json.dumps(manifest,ensure_ascii=False,indent=2).encode('utf-8'))
with zipfile.ZipFile(out) as z:
 assert z.testzip() is None
 for row in manifest:
  data=z.read(row['path']);assert len(data)==row['bytes'] and hashlib.sha256(data).hexdigest()==row['SHA256']
receipt={'at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'path':str(out),'bytes':out.stat().st_size,'SHA256':hashlib.sha256(out.read_bytes()).hexdigest(),'manifest_files':len(manifest),'CRC_and_all_entries_verified':True,'manifest_scope':'named evidence files; manifest and reading guide are container metadata'}
(B/'EVALUATION_EVIDENCE_INTEGRITY.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(receipt,ensure_ascii=False))
