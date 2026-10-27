import pathlib,pandas as pd,json,hashlib
root=pathlib.Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');old=root/'phase97_normal_background_20261020';p=root/'phase98_reference_thresholds_20261020'
files=[old/'tables/bone_reference_pseudobulk.tsv',old/'tables/lung_reference_pseudobulk.tsv',old/'tables/normal_background_all_results.tsv']
(p/'resources/INPUT_MANIFEST.json').write_text(json.dumps([dict(path=str(f),sha256=hashlib.sha256(f.read_bytes()).hexdigest()) for f in files],indent=2),encoding='utf-8')
coverage=[]
for tissue,name,pool,outname in [('bone','bone_reference_pseudobulk.tsv','Fibro-MSC','bone_reference_with_coverage_sensitivity.tsv'),('lung','lung_reference_pseudobulk.tsv','all_fibroblasts','lung_reference_pseudobulk.tsv')]:
 d=pd.read_csv(old/'tables'/name,sep='\t',encoding='utf-8',dtype={'donor':str});d=d[d.specification.eq(pool)].copy()
 assert not d.duplicated(['donor','gene']).any()
 frames=[]
 for n in [1,20,50,100]:
  z=d[d.cells.ge(n)].copy();z['specification']='min'+str(n);frames.append(z)
  q=z.drop_duplicates('donor')
  coverage.append(dict(tissue=tissue,min_cells=n,donors=len(q),selected_cells=int(q.cells.sum()),donor_ids=';'.join(q.donor),donor_cells=';'.join(q.cells.astype(str)),excluded_donors=';'.join(sorted(set(d.donor)-set(q.donor)))))
 pd.concat(frames).to_csv(p/'tables'/outname,sep='\t',encoding='utf-8',index=False,na_rep='NA')
pd.DataFrame(coverage).to_csv(p/'tables/reference_threshold_coverage.tsv',sep='\t',encoding='utf-8',index=False)
print(pd.DataFrame(coverage).to_string(index=False))
