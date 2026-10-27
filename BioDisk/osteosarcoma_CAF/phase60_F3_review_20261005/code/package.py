from pathlib import Path
import pandas as pd,json,hashlib,zipfile,datetime
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase60_F3_review_20261005'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
assert (O/'ANALYSIS_COMPLETE').exists() and (O/'analysis.exit_code').read_text().strip()=='0' and (O/'QA_COMPLETE').exists();qa=json.loads((O/'delivery/QA_ACCEPTANCE.json').read_text());assert qa['background_MT_canonical_mapping_corrected']
c=pd.read_csv(O/'tables/F3_cells_review_QC.tsv.gz',sep='\t');assert c.lineage_stratum.value_counts().to_dict()=={'Fibro':166,'Osteo':69,'Other_unresolved':16}
f=pd.read_csv(O/'tables/F4_two_main_strata_with_residual.tsv.gz',sep='\t');assert f.analysis_group.value_counts().to_dict()=={'F4_Osteo_like':3300,'F4_other_unresolved':1095,'F4_Fibro_like':909}
t=pd.read_csv(O/'tables/F3_threshold_sensitivity_summary.tsv',sep='\t');raw=t[(t.matrix=='raw')&(t.threshold_rule=='at_least_1')];assert raw.MHCII_fibro_co.sum()==12 and raw.cycle_ge2.sum()==57
report=O/'delivery/CAF_F3深入核查与F4表达分层报告_20261005.txt';assert report.stat().st_size>5000;assert len(list((O/'figures').glob('*.png')))==5
acc={'at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'phase':O.name,'status':'complete','F3_cells':251,'F4_cells':5304,'F4_strata':f.analysis_group.value_counts().to_dict(),'original_labels_preserved':True,'new_clustering':False,'new_CNV':False,'all_five_figures_visually_checked':True,'QA':qa,'report':str(report),'findings':['F4 has two main expression strata plus1095 residual cells; not newly validated subclusters','242/251 F3 cells from two patients; original matched reference only2/3 cells','MYC/OXPHOS advantage depends on sparse historical comparator; not retained against UMI-matched broader Fibro','MHCII-fibro codetection12/251; local cycling57/251; no whole-F3 apCAF or proliferating-CAF identity','Elevated MT fractions and source-specific DecontX effects retained; F3 not blanket deleted as low quality']};(O/'delivery/FINAL_ACCEPTANCE.json').write_text(json.dumps(acc,ensure_ascii=False,indent=2));a=json.loads((O/'delivery/ANALYSIS_ACCEPTANCE.json').read_text());a['QA']='passed';(O/'delivery/ANALYSIS_ACCEPTANCE.json').write_text(json.dumps(a,indent=2))
(O/'delivery/WORKFLOW_REPAIRS.json').write_text(json.dumps({'background_QC':'Initial MT diagnostic used absent symbol column; corrected using frozen canonical Ensembl mapping in sensitivity.py. Final background_QC.tsv and SUMMARY.json use corrected data.','fractional_counts':'Nonzero and0.5 detection thresholds added after reviewing fractional corrected counts; original>=1 threshold retained. No inference P values changed.','SSH':'Transient connection closures resolved using a persistent task-owned SSH session and ordinary SCP transfers.'},indent=2))
files=[p for folder in ['code','resources','tables','objects','QA','figures','delivery'] for p in (O/folder).rglob('*') if p.is_file() and p.suffix!='.zip' and p.name not in ['FILES_SHA256.tsv','ZIP_ACCEPTANCE.json','LOCAL_DELIVERY.json'] and '__pycache__' not in str(p)];m=pd.DataFrame([{'path':str(p.relative_to(O)),'bytes':p.stat().st_size,'SHA256':sha(p)} for p in sorted(files)]);m.to_csv(O/'delivery/FILES_SHA256.tsv',sep='\t',index=False);files.append(O/'delivery/FILES_SHA256.tsv');archive=O/'delivery/CAF_F3_review_F4_strata_20261005.zip'
with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for p in files:z.write(p,str(p.relative_to(O)))
with zipfile.ZipFile(archive) as z:
 assert z.testzip() is None
 for r in m.itertuples(index=False):
  d=z.read(r.path);assert len(d)==r.bytes and hashlib.sha256(d).hexdigest()==r.SHA256
receipt={'archive':str(archive),'bytes':archive.stat().st_size,'SHA256':sha(archive),'files':len(files),'CRC_pass':True,'all_entries_SHA256_verified':True};(O/'delivery/ZIP_ACCEPTANCE.json').write_text(json.dumps(receipt,indent=2));print(json.dumps(receipt),flush=True)
