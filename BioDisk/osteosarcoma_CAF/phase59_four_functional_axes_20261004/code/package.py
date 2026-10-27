from pathlib import Path
import json,hashlib,zipfile,datetime,pandas as pd
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase59_four_functional_axes_20261004'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
assert (O/'ANALYSIS_COMPLETE').exists() and (O/'analysis.exit_code').read_text().strip()=='0' and (O/'QA_COMPLETE').exists()
qa=json.loads((O/'delivery/QA_ACCEPTANCE.json').read_text());report=O/'delivery/CAF四功能轴与F2F4信号来源分析报告_20261004.txt';assert report.stat().st_size>5000
figs=sorted((O/'figures').glob('*.png'));assert len(figs)==7
acc={'phase':'phase59_four_functional_axes_20261004','status':'complete','at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'authorized_scope':'F1-F5 four overlapping axes; F2 contractile/matrix; F4 osteogenic source attribution','cells':8174,'sources':9,'programs':27,'union_genes':1300,'all_labels_preserved':True,'new_CNV':False,'new_clustering':False,'spatial_eggshell_geometry_performed':False,'treatment_resistance_inference_performed':False,'statistics_QA':qa,'figures':7,'all_figures_visually_reviewed':True,'key_findings':['F2 contractile marker increase retained; broad contraction and core matrisome not nominally elevated in main source','F4 bone-mineralization transcriptional set FDR significant in full historical class; not confirmed physical mineralization','F4 osteogenic transcripts mainly in expression-reviewed Osteo stratum; full osteogenic/mineralization programs not nominally elevated within Fibro-only comparison','OMD/OGN in F4 Fibro-only nominal signal, not FDR significant; no exclusive CAF marker claim','Cross-source disagreement and low-N contexts retained'], 'report':str(report)}
(O/'delivery/FINAL_ACCEPTANCE.json').write_text(json.dumps(acc,ensure_ascii=False,indent=2))
# Keep failed QA summary attempt as documented repair, without claiming the analysis was re-run.
(O/'delivery/WORKFLOW_REPAIRS.json').write_text(json.dumps({'repairs':[{'component':'QA read performance','action':'cache NPZ arrays once; initial QA process terminated and restarted before result acceptance'},{'component':'summary export','action':'rename count-matrix index to count_matrix to distinguish matrix functional axis','statistics_changed':False,'independent_R_pass_preserved':True},{'component':'co-detection','action':'record per-panel coverage eligibility; no missing measurements called negative','patient_inference_changed':False}],'analysis_exit_code':0},indent=2))
files=[p for folder in ['code','resources','tables','objects','QA','figures','delivery'] for p in (O/folder).rglob('*') if p.is_file() and p.suffix!='.zip' and p.name not in ['FILES_SHA256.tsv','ZIP_ACCEPTANCE.json','LOCAL_DELIVERY.json'] and '__pycache__' not in str(p)]
manifest=pd.DataFrame([{'path':str(p.relative_to(O)),'bytes':p.stat().st_size,'SHA256':sha(p)} for p in sorted(files)]);manifest.to_csv(O/'delivery/FILES_SHA256.tsv',sep='\t',index=False);files.append(O/'delivery/FILES_SHA256.tsv')
archive=O/'delivery/CAF_four_axes_F2F4_results_20261004.zip'
with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for p in files:z.write(p,str(p.relative_to(O)))
with zipfile.ZipFile(archive) as z:
 assert z.testzip() is None
 for r in manifest.itertuples(index=False):
  data=z.read(r.path);assert len(data)==r.bytes and hashlib.sha256(data).hexdigest()==r.SHA256
receipt={'archive':str(archive),'bytes':archive.stat().st_size,'SHA256':sha(archive),'files':len(files),'CRC_pass':True,'all_manifest_entries_SHA256_verified':True};(O/'delivery/ZIP_ACCEPTANCE.json').write_text(json.dumps(receipt,indent=2));print(json.dumps(receipt),flush=True)
