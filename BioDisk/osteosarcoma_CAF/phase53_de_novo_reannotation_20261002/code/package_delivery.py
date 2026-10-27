from pathlib import Path
import json,hashlib,re,subprocess,zipfile,datetime
R=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase53_de_novo_reannotation_20261002');B=R/'blind_review';D=R/'delivery';U=R/'unblinding'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for v in iter(lambda:f.read(8*1024*1024),b''):h.update(v)
 return h.hexdigest()
def dump(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2))
maps=json.loads((U/'identifier_maps.json').read_text())
pat=re.compile('|'.join(re.escape(x) for x in sorted(list(maps['sources'])+list(maps['samples'])+list(maps['patients']),key=len,reverse=True)))
pdfs=sorted(B.glob('D*/*.pdf'));assert len(pdfs)==18
for p in pdfs:
 s=subprocess.run(['pdftotext',str(p),'-'],capture_output=True,text=True,check=True).stdout
 assert not pat.search(s),(p,pat.search(s).group() if pat.search(s) else '')
 assert not any(x in s for x in ['proposed_lineage','HPCA_label','BlueprintEncode_label','supported_class_id'])
dump(D/'FINAL_FIGURE_LEAKAGE_QA.json',{'at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'final_PDF_count':18,'passed':True,'files':[{'path':str(p.relative_to(B)),'SHA256':sha(p)} for p in pdfs]})
pngs=sorted(B.glob('D*/*.png'));assert len(pngs)==18
dump(D/'FIGURE_VISUAL_QA.json',{'final_PNG_count':18,'inspection':'All final PNGs viewed by analysis assistant; patient legends, cluster labels, chromosome labels, clipping and missing-value display inspected. Not an independent biological review.','passed':True,'files':[{'path':str(p.relative_to(B)),'SHA256':sha(p)} for p in pngs]})
stages={'references':'REFERENCES_COMPLETE','structure_pipeline':'STRUCTURE_PIPELINE_COMPLETE','annotation_pipeline':'ANNOTATION_PIPELINE_COMPLETE','copykat_export':'COPYKAT_PROFILES_COMPLETE','infercnv_region_export':'INFER_CNV_REGIONS_COMPLETE','raw_blinding':'RAW_BLINDING_COMPLETE','integration':'ANNOTATION_INTEGRATION_COMPLETE','blind_view_completion':'BLIND_VIEWS_COMPLETE','blind_final_render':'BLIND_FIGURES_FINAL_COMPLETE','final_data_QA':'FINAL_DATA_QA_COMPLETE'}
checks={}
for stage,marker in stages.items():
 code=int((R/(stage+'.exit_code')).read_text().strip());assert code==0 and (R/marker).exists()
 checks[stage]={'exit_code':code,'completion_marker':marker}
for p in (R/'logs').glob('*.exit_code'):assert int(p.read_text().strip())==0,p
assert sha(R/'resources/CURATION_LOCKED.json')=='7f1730f5b351a4ec3a7ca1228b4449e92eec3c4b16ac7383662c8db2de20674a'
dump(D/'FINAL_ACCEPTANCE.json',{'state':'reannotation_v1_ready_for_independent_blind_review','cells':876648,'samples':107,'sources':9,'primary_clusters':307,'stages':checks,'historical_failure_retained':'blind_views.exit_code=1; repaired by blind_view_completion.exit_code=0; see BLIND_EXPORT_REPAIR.json','independent_blind_review_performed':False,'raw_counts_in_compact_zip':False,'raw_counts_server_directory':str(B/'raw_counts'),'visual_and_final_PDF_QA_passed':True})
# Separate manifests include only files admitted to that package.
blind=sorted(p for p in B.rglob('*') if p.is_file() and 'raw_counts' not in p.relative_to(B).parts)
assert not any(pat.search(str(p.relative_to(B))) for p in blind)
private=[]
for folder in ['unblinding','code']:
 private.extend(p for p in (R/folder).rglob('*') if p.is_file() and '__pycache__' not in p.parts)
private.extend(p for p in (R/'resources').iterdir() if p.is_file() and p.suffix!='.rds')
private.extend(p for p in (R/'tables').iterdir() if p.is_file())
for d in (R/'tables').iterdir():
 if d.is_dir():
  private.extend(p for p in d.iterdir() if p.is_file() and p.name in ['reference_annotations.tsv','reference_profile_metadata.tsv','cluster_evidence_summary.tsv','STRUCTURE_QA.json','primary_clusters.tsv','structure_features.tsv'])
private.extend(p for p in (R/'logs').iterdir() if p.is_file() and p.suffix in ['.log','.exit_code'])
private.extend(p for p in R.iterdir() if p.is_file() and (p.name.endswith('.exit_code') or p.name.endswith('_COMPLETE')))
refs=[{'server_path':str(p),'bytes':p.stat().st_size,'SHA256':sha(p)} for p in (R/'resources').glob('*.rds')]
dump(D/'REFERENCE_OBJECT_MANIFEST.json',refs)
private.extend(p for p in D.iterdir() if p.is_file() and p.suffix in ['.json','.tsv','.txt'] and not p.name.startswith('PACKAGE_'))
def pack(files,base,out):
 files=sorted(set(files));entries=[]
 for p in files:entries.append({'path':str(p.relative_to(base)),'bytes':p.stat().st_size,'SHA256':sha(p)})
 manifest='\n'.join(['relative_path\tbytes\tSHA256']+[x['path']+'\t'+str(x['bytes'])+'\t'+x['SHA256'] for x in entries])+'\n'
 part=out.with_suffix('.zip.partial')
 with zipfile.ZipFile(part,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=4,allowZip64=True) as z:
  for p,e in zip(files,entries):z.write(p,e['path'])
  z.writestr('PACKAGE_FILE_MANIFEST.tsv',manifest)
 with zipfile.ZipFile(part) as z:
  assert z.testzip() is None
  for e in entries:
   assert z.getinfo(e['path']).file_size==e['bytes']
   h=hashlib.sha256()
   with z.open(e['path']) as f:
    for v in iter(lambda:f.read(8*1024*1024),b''):h.update(v)
   assert h.hexdigest()==e['SHA256']
 part.replace(out)
 return {'archive':str(out),'bytes':out.stat().st_size,'SHA256':sha(out),'files':len(entries),'CRC_and_all_entry_SHA256_verified':True}
result={'blind':pack(blind,B,D/'CAF重注释_盲审材料_20261002.zip'),'private':pack(private,R,D/'CAF重注释_新注释与揭盲对应表_20261002.zip')}
dump(D/'PACKAGE_RECEIPTS.json',result);(R/'PACKAGES_COMPLETE').write_text('SUCCESS')
print(json.dumps(result,ensure_ascii=False),flush=True)
