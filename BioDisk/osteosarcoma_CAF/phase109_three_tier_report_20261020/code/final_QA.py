import pathlib,json,hashlib,subprocess,shutil
import pandas as pd,numpy as np
R=pathlib.Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');P=R/'phase109_three_tier_report_20261020'
c=pd.read_csv(P/'tables/tier_candidate_registry.tsv',sep='\t')
assert [int(c[k].sum()) for k in ['tier1','tier2','tier3']]==[14,8,6]
assert not (c.tier3&~c.tier2).any()
d=pd.read_csv(P/'tables/final_six_complete_statistics.tsv',sep='\t').set_index('gene')
src=pd.read_csv(R/'phase108_each_source_significant_20261020/tables/final_six_genes.tsv',sep='\t').set_index('gene')
cols=[k for k in src.columns if k!='feature']
assert np.allclose(d.loc[src.index,cols],src[cols],rtol=0,atol=1e-15)
for name in ['effect_pooled','effect_NCI','source17_effect','source23_effect','mean_lung_minus_bone','third_effect']:assert d[name].gt(0).all()
assert (d.FDR_pooled.lt(.05)&d.FDR_NCI.lt(.05)).sum()==5
assert d.loc['PNKD','P_pooled']<.05 and d.loc['PNKD','P_NCI']<.05 and d.loc['PNKD','FDR_NCI']>.05
bio=json.loads((P/'resources/previous_six_gene_biology.json').read_text())
assert len(bio['genes'])==6 and len(bio['references'])==10
for item in json.loads((P/'resources/INPUT_MANIFEST.json').read_text()):
 assert hashlib.sha256(pathlib.Path(item['path']).read_bytes()).hexdigest()==item['sha256']
plots=json.loads((P/'resources/PLOT_QA.json').read_text());assert len(plots)==6 and all(x['label_collisions']==0 and x['labels']==6 for x in plots)
text=(P/'render/report_text.txt').read_text(encoding='utf-8')
for k in ['PNKD','TOR1AIP1','GRID2','TAF10','BEST1','RGS2-AS1','0.0715','6,633','65,536','95%','阴性','GACAT2','MR1','p53']:assert k in text,k
pdf=P/'report/CAF_three_tier_six_genes_biology_20261020.pdf'
info=subprocess.run(['pdfinfo',str(pdf)],capture_output=True,text=True,check=True).stdout
assert 'Pages:           8' in info
q=dict(status='PASS',pages=8,figures=6,all_pages_visually_inspected=True,revised_page2_reinspected=True,source_statistics_unchanged=True,all_sources_hash_unchanged=True,tier_counts=[14,8,6],both_background_FDR_genes=5,PNKD_nominal_only_for_NCI_background=True,biology_source_preserved=True,key_negative_evidence_preserved=True,plot_labels_no_overlap=True,PDF_text_checked=True)
(P/'resources/FINAL_QA.json').write_text(json.dumps(q,indent=2),encoding='utf-8')
status=dict(status='complete',phase=str(P),report=str(pdf),local_report=r'C:\Users\31937\Downloads\CAF_three_tier_20261020\CAF_three_tier_six_genes_biology_20261020.pdf',sha256=hashlib.sha256(pdf.read_bytes()).hexdigest(),tiers=[14,8,6],final_genes=d.index.tolist(),pages=8,volcanoes=6,literature_references=10,QA='PASS',original_reports_preserved=True)
(P/'FINAL_STATUS.json').write_text(json.dumps(status,ensure_ascii=False,indent=2),encoding='utf-8')
snap=P/'resources/previous_current';snap.mkdir(exist_ok=True)
for name in ['CAF_RESUME_STATUS_CURRENT.json','CAF_PATIENT_VALIDATION_CURRENT.json']:
 f=R/name
 if f.exists() and not (snap/name).exists():shutil.copyfile(f,snap/name)
(R/'CAF_THREE_TIER_REPORT_CURRENT.json').write_text(json.dumps(status,ensure_ascii=False,indent=2),encoding='utf-8')
f=R/'CAF_RESUME_STATUS_CURRENT.json';old=json.loads(f.read_text());old.update(authoritative_current_file=str(R/'CAF_THREE_TIER_REPORT_CURRENT.json'),current_phase=str(P),current_report=str(pdf),latest_analysis_report=str(pdf),status='three_tier_PDF_biology_report_complete',active_task='Standalone8-page report14to8to6 with6 whole-gene volcanoes and6 previous biology explanations completed')
f.write_text(json.dumps(old,ensure_ascii=False,indent=2),encoding='utf-8')
(R/'CAF_PATIENT_VALIDATION_CURRENT.json').write_text(json.dumps(status,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(status,ensure_ascii=False))
