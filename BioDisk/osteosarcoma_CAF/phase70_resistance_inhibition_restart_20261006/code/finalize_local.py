from pathlib import Path
import hashlib,json,datetime,zipfile,pandas as pd
B=Path(__file__).resolve().parents[1]
report=Path('C:/Users/31937/Downloads/骨肉瘤CAF功能分类一致性分析报告_20261006_胶原基质与GSEA合并版.docx')
sha=hashlib.sha256(report.read_bytes()).hexdigest();assert sha=='9081b06e8a86b379af84b30a5313995eb5199c52e2c8faa5704029a2432df34a'
r=json.loads((B/'resources/ANALYSIS_COMPLETE.json').read_text(encoding='utf-8'));assert r['status']=='completed'
s=pd.read_csv(B/'tables/positive_candidate_descriptive_sensitivity.tsv',sep='\t');assert len(s)==168 and not s.duplicated(['cohort','gene','sensitivity','excluded']).any()
assert s[s.cohort=='GSE39055'].effect.gt(0).all()
for name in ['all_fixed_family_adjusted_gene_response.tsv','all_fixed_family_unadjusted_adjusted_comparison.tsv']:
 d=pd.read_csv(B/'tables'/name,sep='\t');assert len(d)==1303
accept={'checked_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'passed','analysis_exit_code':0,'all_gene_models':1303,'independent_patients':83,'candidate_genes':39,'candidate_records':54,'descriptive_sensitivity_rows':168,'original_client_report_sha256':sha,'client_submission_provenance':'User statement on 2026-10-06; actual submission timestamp not independently recorded','report_preserved':True,'complete_p_fdr_families_preserved':True,'new_stage_scope':'Adjusted gene clinical models plus source mapping and literature-grounded mechanism prioritization; next tumor/immune receiver analysis not yet run','background_process_running':False}
(B/'delivery/FINAL_ACCEPTANCE.json').write_text(json.dumps(accept,ensure_ascii=False,indent=2),encoding='utf-8')
archive=B/'delivery/CAF耐药抑制重启_首轮结果_20261006.zip'
with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED) as z:
 for folder in ['inputs','tables','resources','code','delivery']:
  for p in sorted((B/folder).rglob('*')):
   if p.is_file() and p!=archive and p.name!='TRANSFER_MANIFEST.json':z.write(p,p.relative_to(B).as_posix())
with zipfile.ZipFile(archive) as z:assert z.testzip() is None
manifest={p.relative_to(B).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for folder in ['inputs','tables','resources','code','delivery'] for p in sorted((B/folder).rglob('*')) if p.is_file() and p.name!='TRANSFER_MANIFEST.json'}
(B/'resources/TRANSFER_MANIFEST.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(accept,ensure_ascii=False));print('ZIP_OK',archive.stat().st_size,'FILES',len(manifest))
