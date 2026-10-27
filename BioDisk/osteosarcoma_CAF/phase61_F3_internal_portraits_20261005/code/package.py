from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,zipfile,datetime
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase61_F3_internal_portraits_20261005'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
assert (O/'ANALYSIS_COMPLETE').exists() and (O/'QA_COMPLETE').exists();assert (O/'prepare.exit_code').read_text().strip()=='0' and (O/'gsea.exit_code').read_text().strip()=='0';report=O/'delivery/F3内部亚群表达与功能画像_20261005.txt';assert report.stat().st_size>5000;labels=pd.read_csv(O/'tables/F3_internal_working_labels.tsv.gz',sep='\t');assert len(labels)==251 and labels.cell_id.is_unique and labels.internal_working_description.notna().all();assert labels.state.value_counts().to_dict()=={'C0':119,'C2':90,'C1':42};assert len(list((O/'figures').glob('*.png')))==4
cons=[]
for name in ['C0_vs_other_F3','C1_vs_other_F3','C2_vs_other_F3','C1_vs_C2_same_patient']:
 a=pd.read_csv(O/'tables/genes'/f'raw__{name}.tsv.gz',sep='\t');b=pd.read_csv(O/'tables/genes'/f'corrected__{name}.tsv.gz',sep='\t');assert a.symbol.equals(b.symbol);valid=a.any_expression&a.mean_delta.ne(0);cons.append({'contrast':name,'raw_expressed_nonzero_effect_genes':int(valid.sum()),'direction_preserved':int((np.sign(a.loc[valid,'mean_delta'])==np.sign(b.loc[valid,'mean_delta'])).sum())})
pd.DataFrame(cons).to_csv(O/'tables/raw_corrected_direction_summary.tsv',sep='\t',index=False)
acc={'phase':O.name,'status':'complete','at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'user_priority':'Coarse CAF working labels; internal subgroup expression and functional portraits are the mainline; other-F specificity and cross-patient validation are not prerequisites for exploratory internal analysis','original_F3_cells':251,'main_source_cells':245,'states':labels.state.value_counts().to_dict(),'original_state_ids_preserved':True,'new_clustering':False,'new_CNV':False,'contrasts':12,'GSEA_jobs':36,'GSEA_rows':21971,'figures':4,'all_figures_visually_checked':True,'QA':json.loads((O/'delivery/QA_ACCEPTANCE.json').read_text()),'report':str(report),'working_descriptions':{'C0':'FBLN2/COL4-associated inflammatory-response state','C1':'ISLR/CLU/LOX matrix-component state','C2':'MFAP2/ANXA2 MYC-translation state'},'key_result':'C1/C2 share matrix background; within-patient differences support state characterization; program-model P values kept separate from patient inference'}
(O/'delivery/FINAL_ACCEPTANCE.json').write_text(json.dumps(acc,ensure_ascii=False,indent=2));(O/'delivery/WORKFLOW_REPAIRS.json').write_text(json.dumps({'NPZ_string_arrays':'Converted task-generated object strings to Unicode arrays; cache now readable without pickle; expression values unchanged.','approval_timeout':'Read-only result query retried once successfully after automatic review timeout.'},indent=2))
files=[p for folder in ['code','resources','objects','tables','QA','figures','delivery'] for p in (O/folder).rglob('*') if p.is_file() and p.suffix!='.zip' and p.name not in ['FILES_SHA256.tsv','ZIP_ACCEPTANCE.json','LOCAL_DELIVERY.json'] and '__pycache__' not in str(p)];m=pd.DataFrame([{'path':str(p.relative_to(O)),'bytes':p.stat().st_size,'SHA256':sha(p)} for p in sorted(files)]);m.to_csv(O/'delivery/FILES_SHA256.tsv',sep='\t',index=False);files.append(O/'delivery/FILES_SHA256.tsv');archive=O/'delivery/CAF_F3_internal_portraits_20261005.zip'
with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for p in files:z.write(p,str(p.relative_to(O)))
with zipfile.ZipFile(archive) as z:
 assert z.testzip() is None
 for r in m.itertuples(index=False):
  d=z.read(r.path);assert len(d)==r.bytes and hashlib.sha256(d).hexdigest()==r.SHA256
receipt={'archive':str(archive),'bytes':archive.stat().st_size,'SHA256':sha(archive),'files':len(files),'CRC_pass':True,'all_entries_SHA256_verified':True};(O/'delivery/ZIP_ACCEPTANCE.json').write_text(json.dumps(receipt,indent=2));print(json.dumps(receipt),flush=True)
