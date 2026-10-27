from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
import json,hashlib,gc,traceback
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
R=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=R/'phase113_marker_function_heatmaps_20261021'
for s in ['code','tables','resources','logs','delivery']:(O/s).mkdir(parents=True,exist_ok=True)
def dump(p,v):(O/p).write_text(json.dumps(v,ensure_ascii=False,indent=2,default=str),encoding='utf-8')
def read(p,**kw):return pd.read_csv(R/p,sep='\t',**kw)
panels=json.loads((R/'phase27_identity_review_20260926/resources/lineage_panels.json').read_text())
cp={'Fibroblast':['DCN','LUM','PDGFRA'],'Mural':['RGS5','CSPG4','NOTCH3'],'Osteogenic':['RUNX2','SP7','ALPL'],'T_NK':['CD3D','TRAC','NKG7'],'B_Plasma':['CD79A','MS4A1','JCHAIN'],'Myeloid':['LST1','TYROBP','CD68'],'Endothelial':['PECAM1','VWF','CDH5'],'Epithelial':['EPCAM','KRT8','KRT18'],'Mast':['TPSAB1','TPSB2','CPA3'],'Myogenic':['MYOD1','MYOG','DES'],'Osteoclast':['ACP5','CTSK','DCSTAMP'],'Chondrocyte':['ACAN','COL2A1','COMP'],'Erythroid':['HBA1','HBB','GYPA']}
assert all(set(v)<=set(panels[k]) for k,v in cp.items())
cgenes=[g for gs in cp.values() for g in gs];assert len(set(cgenes))==39
cfg=json.loads((R/'phase66_doublet_repair_DSC1_20261005/resources/PROGRAMS.json').read_text())
desired={'collagen_components':['COL1A1','COL1A2','COL3A1','COL6A1'],'collagen_biosynthesis':['P4HA1','PLOD2','SERPINH1','PCOLCE'],'collagen_crosslinking':['LOX','LOXL1','LOXL2','LOXL4'],'ECM_degradation':['MMP2','MMP14','CTSK','FN1'],'contraction':['ACTA2','TAGLN','MYL9','MYH11'],'mineralization':['RUNX2','ALPL','IBSP','SPP1'],'inflammation':['IL6','CCL2','CXCL8','ICAM1']}
kp={k:[g for g in gs if g in cfg['sets'][k]] for k,gs in desired.items()};assert all(len(gs)>=2 for gs in kp.values()),kp
six=['PNKD','TOR1AIP1','GRID2','TAF10','BEST1','RGS2-AS1'];kp['Candidates']=six
dump('resources/PANELS.json',dict(C=cp,K=kp,K_selected_not_in_frozen_set={k:sorted(set(gs)-set(kp[k])) for k,gs in desired.items()},programs=cfg))
meta=read('phase112_allcell_3D_landscape_20261021/tables/all_cell_metadata.tsv.gz')
canon=read('phase58_F1F5_round1_20261004/resources/source_canonical_feature_map.tsv');maps={src:dict(zip(z.symbol,z.feature)) for src,z in canon.groupby('dataset')}
cm=read('phase81_allCAF_bone_lung_20261012/tables/clinical_sample_manifest.tsv').set_index('sample')
cr=[];maprows=[];sampleledger=[]
for sm,q in meta.groupby('sample',sort=True):
 cache=O/'tables'/('C_sample_'+sm+'.tsv.gz')
 if cache.exists():cr.append(pd.read_csv(cache,sep='\t'));continue
 src=q.dataset.iloc[0];p=Path(cm.loc[sm,'full_called_counts_h5ad']);a=ad.read_h5ad(p,backed='r')
 ix=a.obs_names.get_indexer(q.cell_id);assert (ix>=0).all() and len(q)==a.n_obs
 raw=a.X[:].tocsr();lib=np.asarray(raw.sum(1)).ravel();assert (lib>0).all();assert np.allclose(lib,a.obs.n_counts.to_numpy())
 rows=[]
 for gene in cgenes:
  feature=maps[src].get(gene);j=a.var_names.get_indexer([feature])[0] if feature else -1
  maprows.append(dict(source=src,sample=sm,gene=gene,feature=feature,present=j>=0))
  if j<0:
   vals=np.full(len(q),np.nan);det=vals.copy()
  else:
   counts=np.asarray(raw[:,j].toarray()).ravel()[ix];vals=np.log1p(counts/lib[ix]*1e4);det=(counts>0).astype(float)
  temp=q[['donor_record','cell_type']].copy();temp['log']=vals;temp['det']=det
  for (donor,typ),z in temp.groupby(['donor_record','cell_type']):
   rows.append(dict(source=src,sample=sm,donor=donor,cell_type=typ,gene=gene,cells=len(z),measured_cells=int(z['log'].notna().sum()),sum_log=float(z['log'].sum()),detected_cells=float(z.det.sum())))
 t=pd.DataFrame(rows);t.to_csv(cache,sep='\t',index=False);cr.append(t)
 sampleledger.append(dict(sample=sm,cells=len(q),library_UMIs=float(lib.sum()),path=str(p),bytes=p.stat().st_size))
 a.file.close();del raw,a;gc.collect();print('C_SAMPLE',sm,len(q),flush=True)
c=pd.concat(cr,ignore_index=True);pat=c.groupby(['source','donor','cell_type','gene'],as_index=False)[['cells','measured_cells','sum_log','detected_cells']].sum()
pat['mean_log1pCP10k']=np.where(pat.measured_cells.gt(0),pat.sum_log/pat.measured_cells,np.nan);pat['detection_percent']=np.where(pat.measured_cells.gt(0),pat.detected_cells/pat.measured_cells*100,np.nan)
# Incomplete measured coverage is flagged rather than silently interpreted as all cells measured.
pat['complete_coverage']=pat.cells.eq(pat.measured_cells)
pat.to_csv(O/'tables/C_patient_type_gene.tsv.gz',sep='\t',index=False)
csummary=pat.groupby(['source','cell_type','gene']).agg(mean_log1pCP10k=('mean_log1pCP10k','mean'),detection_percent=('detection_percent','mean'),patients=('donor','nunique'),measured_patients=('mean_log1pCP10k','count'),cells=('cells','sum'),measured_cells=('measured_cells','sum')).reset_index()
csummary.to_csv(O/'tables/C_dotplot_values.tsv',sep='\t',index=False)
if maprows:pd.DataFrame(maprows).to_csv(O/'resources/C_feature_mapping.tsv',sep='\t',index=False)
if sampleledger:pd.DataFrame(sampleledger).to_csv(O/'resources/C_input_ledger.tsv',sep='\t',index=False)
print('C_COMPLETE',len(csummary),flush=True)
# K uses precisely the old-source subset of the accepted 6633-cell candidate figures.
m=read('phase111_reference_figure_adaptation_20261021/tables/Fig1i_cell_values.tsv.gz');m=m[m.source.isin(['SCPCP000017','SCPCP000023'])].copy();assert len(m)==5103
kg=[g for gs in kp.values() for g in gs];assert len(set(kg))==len(kg)
registry=json.loads((R/'phase96_all14_literature_20261016/project_registry.json').read_text());sixmap={g['gene']:g['feature_ids'] for g in registry}
a=ad.read_h5ad(R/'phase41_full_transcriptome_setup_20261001/objects/all_candidate_full_genes.h5ad',backed='r');ix=a.obs_names.get_indexer(m.cell_id);assert (ix>=0).all()
kraw=a.X[:].tocsr()
presence=sparse.load_npz(R/'phase41_full_transcriptome_setup_20261001/objects/sample_feature_presence.npz').tocsr();si=read('phase41_full_transcriptome_setup_20261001/tables/feature_presence_sample_index.tsv').set_index('sample').presence_row
vals=np.full((len(m),len(kg)),np.nan);coverage=[]
for j,gene in enumerate(kg):
 for src in m.source.unique():
  mask=m.source.eq(src).to_numpy();fid=sixmap[gene] if gene in six else maps[src].get(gene);gi=a.var_names.get_indexer([fid])[0] if fid else -1
  assert gi>=0,(src,gene,fid)
  measured=presence[:,gi].toarray().ravel().astype(bool)[si.loc[m.loc[mask,'sample']].to_numpy()]
  count=kraw[:,gi].toarray().ravel()[ix[mask]];z=np.log1p(count/m.loc[mask,'full_library_counts'].to_numpy()*1e4);z[~measured]=np.nan;vals[mask,j]=z
  coverage.append(dict(source=src,gene=gene,feature=fid,cells=int(mask.sum()),measured=int(measured.sum())))
  if gene in six:assert np.allclose(z,m.loc[mask,gene+'__log'],equal_nan=True)
a.file.close();pd.DataFrame(coverage).to_csv(O/'resources/K_feature_coverage.tsv',sep='\t',index=False)
cell=m[['cell_id','source','donor','site','class']].reset_index(drop=True).copy()
for j,g in enumerate(kg):cell[g]=vals[:,j]
cell.to_csv(O/'tables/K_gene_cell_values.tsv.gz',sep='\t',index=False)
group=['source','class','donor'];profile=cell.groupby(group,sort=True)[kg].mean();info=cell.groupby(group,sort=True).agg(cells=('cell_id','size'),site=('site','first'));assert cell.groupby(group).site.nunique().eq(1).all()
profile=info.join(profile).reset_index();profile.to_csv(O/'tables/K_gene_patient_class_profiles.tsv',sep='\t',index=False)
scores=read('phase66_doublet_repair_DSC1_20261005/tables/cells_repaired_with_confidence.tsv.gz',usecols=['cell_id']+['raw__'+k for k in cfg['sets']]);scores=m[['cell_id','source','class','donor','site']].merge(scores,on='cell_id',validate='one_to_one');assert len(scores)==5103
sc=scores.groupby(group)[['raw__'+k for k in cfg['sets']]].mean();scoreprofile=info.join(sc).reset_index();scoreprofile.to_csv(O/'tables/K_program_patient_class_profiles.tsv',sep='\t',index=False)
assert profile[['source','class','donor']].equals(scoreprofile[['source','class','donor']])
qa=dict(status='data_complete',C_cells=len(meta),C_sources=meta.groupby('dataset').size().to_dict(),C_markers=len(cgenes),C_incomplete_patient_gene_rows=int((~pat.complete_coverage).sum()),K_cells=len(m),K_genes=len(kg),K_programs=len(cfg['sets']),K_profiles=len(profile),six_candidate_values_reproduced=True,original_data_modified=False,inference='descriptive expression only; no P values; existing functional scores not independent validation of constituent genes')
dump('DATA_QA.json',qa);print(json.dumps(qa),flush=True)
