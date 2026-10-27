from pathlib import Path
import json,hashlib,shutil
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
R=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF')
O=R/'phase111_reference_figure_adaptation_20261021'
for s in ['code','tables','resources','figures','logs']:(O/s).mkdir(parents=True,exist_ok=True)
genes=['PNKD','TOR1AIP1','GRID2','TAF10','BEST1','RGS2-AS1']
registry=json.loads((R/'phase96_all14_literature_20261016/project_registry.json').read_text())
features=[next(g['feature_ids'] for g in registry if g['gene']==s) for s in genes]
inputs=[]
def read(rel,**kw):
 p=R/rel;inputs.append(p);return pd.read_csv(p,sep='\t',**kw)
m=read('phase103_pooled_NCI_sensitivities_20261020/tables/old_cell_membership.tsv.gz');m=m[m.cellset.eq('expanded')].copy()
meta=read('phase66_doublet_repair_DSC1_20261005/tables/cells_repaired_with_confidence.tsv.gz')
m=m.merge(meta[['cell_id','expanded_class_id','class_name','review_lineage_code','doublet_positive_any','full_library_counts']],on='cell_id',validate='one_to_one')
assert len(m)==5103 and m.expanded_class_id.notna().all()
p=R/'phase41_full_transcriptome_setup_20261001/objects/all_candidate_full_genes.h5ad';inputs.append(p)
a=ad.read_h5ad(p,backed='r');ix=a.obs_names.get_indexer(m.cell_id);gi=a.var_names.get_indexer(features);assert (ix>=0).all() and (gi>=0).all()
v=np.column_stack([np.asarray(a.X[:,int(j)].toarray()).ravel()[ix] for j in gi]);a.file.close()
presence=sparse.load_npz(R/'phase41_full_transcriptome_setup_20261001/objects/sample_feature_presence.npz')
si=read('phase41_full_transcriptome_setup_20261001/tables/feature_presence_sample_index.tsv').set_index('sample').presence_row
assert np.asarray(presence[si.loc[m['sample'].unique()].to_numpy()][:,gi].toarray()).astype(bool).all()
m['class']=m.expanded_class_id
def attach(meta,x,lib):
 meta=meta.copy()
 assert x.shape==(len(meta),6) and (lib>0).all()
 for j,g in enumerate(genes):meta[g+'__raw']=x[:,j];meta[g+'__log']=np.log1p(x[:,j]/lib*1e4)
 return meta
parts=[attach(m,v,m.full_library_counts.to_numpy())]
n=read('phase102_NCI_author_CAF_20261020/tables/author_CAF_cell_registry.tsv.gz')
aliases=json.loads((R/'phase97_normal_background_20261020/resources/GENE_ALIAS_MAP.json').read_text())['aliases']
mapping=[]
for sm,z in n.groupby('sample'):
 p=R/'phase95_NCI_expanded_reanalysis_20261016/inputs'/sm/'available_counts.h5ad';inputs.append(p)
 a=ad.read_h5ad(p,backed='r');ix=a.obs_names.get_indexer(z.cell_id);assert (ix>=0).all()
 gi=[]
 for gene,feature in zip(genes,features):
  hit=np.flatnonzero((a.var_names==feature)|a.var.symbol.isin(aliases[gene])|a.var_names.isin(aliases[gene]))
  assert len(hit)==1,(sm,gene,len(hit))
  gi.append(int(hit[0]));mapping.append(dict(sample=sm,gene=gene,canonical_feature=feature,object_feature=a.var_names[hit[0]],object_symbol=a.var.symbol.iloc[hit[0]]))
 v=np.column_stack([np.asarray(a.X[:,int(j)].toarray()).ravel()[ix] for j in gi]);lib=a.obs.iloc[ix].n_counts.to_numpy();a.file.close()
 z=z.assign(source='NCI',donor=z.patient_id,site='Lung',cellset='expanded',**{'class':'Author CAF','full_library_counts':lib})
 parts.append(attach(z,v,lib))
cells=pd.concat(parts,ignore_index=True);assert len(cells)==6633 and cells.cell_id.is_unique
ref=read('phase103_pooled_NCI_sensitivities_20261020/tables/all_patient_expression.tsv');ref=ref[ref.cellset.eq('expanded')&ref.gene.isin(genes)]
err=[]
for (src,donor),q in cells.groupby(['source','donor']):
 for g in genes:
  expected=ref[(ref.source==src)&(ref.donor==donor)&(ref.gene==g)].expression
  assert len(expected)==1
  err.append(abs(np.log2(q[g+'__raw'].sum()/q.full_library_counts.sum()*1e4+.1)-expected.iloc[0]))
assert max(err)<1e-9
rows=[];donors=[]
for (src,cl),q in cells.groupby(['source','class'],sort=True):
 for g in genes:
  vals=[];det=[]
  for donor,z in q.groupby('donor'):
   u=float(z[g+'__log'].mean());v=float(z[g+'__raw'].gt(0).mean()*100)
   vals.append(u);det.append(v);donors.append(dict(source=src,group=cl,gene=g,donor=donor,cells=len(z),mean_log1pCP10k=u,RNA_detected_percent=v))
  rows.append(dict(source=src,group=cl,gene=g,donors=len(vals),cells=len(q),mean_log1pCP10k=float(np.mean(vals)),RNA_detected_percent=float(np.mean(det))))
df=pd.DataFrame(rows);df.to_csv(O/'tables/Fig1i_six_candidates_dotplot.tsv',sep='\t',index=False)
pd.DataFrame(mapping).to_csv(O/'resources/Fig1i_NCI_feature_mapping.tsv',sep='\t',index=False)
pd.DataFrame(donors).to_csv(O/'tables/Fig1i_patient_level_values.tsv',sep='\t',index=False)
cells.to_csv(O/'tables/Fig1i_cell_values.tsv.gz',sep='\t',index=False)
qa=dict(status='PASS',cells=len(cells),patients=cells.donor.nunique(),source_patient_counts=cells.groupby('source').donor.nunique().to_dict(),max_existing_patient_pseudobulk_error=max(err),source_class_groups=len(df)//6,genes=genes,biological_weighting='Each donor equal weight within each source/class; mean of donor means and detection fractions',inference='Descriptive; no cell-level P values',NCI_coverage='Features present in author merged count object; per-library premerge feature masks unavailable',scope='Phase103 expanded original candidate pool plus NCI author CAF; original mixed identity and doublet flags retained; not a purified fibroblast claim')
(O/'resources/Fig1i_QA.json').write_text(json.dumps(qa,indent=2),encoding='utf-8')
hashes=[]
for p in inputs:
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 hashes.append(dict(path=str(p),sha256=h.hexdigest()))
(O/'resources/Fig1i_INPUT_HASHES.json').write_text(json.dumps(hashes,indent=2),encoding='utf-8')
print(json.dumps(qa));print(df[['source','group','cells','donors']].drop_duplicates().to_string(index=False))
