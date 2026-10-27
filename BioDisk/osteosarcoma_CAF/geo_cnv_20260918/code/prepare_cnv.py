from pathlib import Path
import os,json,gzip,hashlib,datetime
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='4'
import numpy as np,pandas as pd,anndata as ad
from scipy.io import mmwrite
from threadpoolctl import threadpool_limits
threadpool_limits(4)
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/geo_cnv_20260918');Q=P.parent/'geo_preprocess_20260918'
M=pd.read_csv(Q/'manifest.tsv',sep='\t').fillna('');pack=pd.read_csv(Q/'review/cluster_review_packets.tsv',sep='\t',dtype={'cluster':str})
PANELS={
'GSE318840':{'OB_Fib':['COL1A1','CDH11','ITGA10','RUNX2','IBSP','FGFR3'],'Myeloid':['CD163','MS4A7','CD14'],'Endo':['EGFL7','CDH5','PECAM1'],'T_NK':['KLRK1','CD3E','NKG7'],'Mural':['PDGFRB','NOTCH3'],'OC':['ACP5','CTSK','MMP9'],'Myogenic':['DES','MYBPC2']},
'GSE162454':{'Osteo':['ALPL','RUNX2','IBSP'],'CAF_marker_program':['COL1A1','FAP','VIM'],'Myeloid':['LYZ','CD68'],'T_NK':['CD2','CD3D','CD3E','CD3G','GNLY','NKG7','KLRD1','KLRB1'],'Endo':['EGFL7','PLVAP'],'B':['MS4A1','CD79A'],'Plasma':['IGHG1','MZB1'],'OC':['ACP5','CTSK']},
'GSE250015':{'Osteo':['COL1A1','RUNX2'],'Fib':['COL3A1','LUM','DCN'],'Mural':['ACTA2','RGS5'],'Endo':['PECAM1','VWF'],'Myeloid':['CD14','LYZ','CD86'],'T_NK':['CD3D','CD3E','CD2'],'B_Plasma':['JCHAIN','MZB1'],'OC':['MMP9','CTSK','ACP5'],'author_Osteocyte_myogenic_markers':['MYF5','CHRNA1'],'Cycling':['MKI67','TOP2A']}}
(P/'review/author_marker_panels.json').write_text(json.dumps(PANELS,indent=2))
feature=next((P.parent/'phase7_review_states_20260917/download_expansion_20260917/files/GSE162454').glob('*features.tsv.gz'))
symids={}
with gzip.open(feature,'rt') as f:
 for line in f:
  z=line.rstrip().split('\t');symids.setdefault(z[1],set()).add(z[0])
s2e={s:next(iter(ids)) for s,ids in symids.items() if len(ids)==1}
go=pd.read_csv(P/'resources/infercnv_gene_order_GRCh38_98.tsv',sep='\t',header=None,names=['id','chr','start','end'],dtype={'chr':str}).set_index('id')
protocol=json.loads((P.parent/'phase6_nonpdx_cnv_20260916/resources/reference_protocol.json').read_text());protocol.update(minimum_total_holdout=20,holdout_reason='require independent heldout calibration; prohibit training-only calibration',doublet_rank='90th percentile of final chosen method score per sample',target_policy='all retained cells except sampled reference and holdout; no candidate-only gating',scPCA='not included')
(P/'resources/reference_protocol.json').write_text(json.dumps(protocol,indent=2))
summ=[];panel=[];transition=[];coverage=[]
for r in M.to_dict('records'):
 sm=r['sample'];a=ad.read_h5ad(Q/'accepted'/sm/'analysis_ready.h5ad');o=a.obs.copy();o['cell_id']=a.obs_names
 ev=pd.read_csv(Q/'review'/(sm+'_cell_evidence.tsv.gz'),sep='\t',index_col='cell_id',dtype={'cluster':str}).loc[a.obs_names]
 assert ev.index.equals(a.obs_names)
 for c in ev.columns:
  if c.startswith(('pctUMI_','det_')):o[c]=ev[c].to_numpy()
 pp=pack[pack['sample']==sm].set_index('cluster');o['CTA_context']=o.cluster.astype(str).map(pp.CTA_context).astype(bool)
 o['phase2_label']=o.lineage_prelim.astype(str);o['lineage_phase3']=o.lineage_reviewed_candidate.astype(str);o['short']=r['short']
 o.groupby(['phase2_label','lineage_phase3'],observed=True).size().rename('cells').reset_index().assign(sample=sm).pipe(lambda v:transition.extend(v.to_dict('records')))
 x=a.layers['counts'].tocsr();xc=x.tocsc();lx=a.X.tocsc();symbols=a.var.symbol.astype(str).to_numpy()
 def count(genes):return np.asarray(xc[:,np.flatnonzero(np.isin(symbols,genes))].sum(1)).ravel()
 for name,genes in PANELS[r['dataset']].items():
  ix=[np.flatnonzero(symbols==g) for g in genes];det=np.zeros(a.n_obs);expr=np.zeros(a.n_obs)
  for ids in ix:
   if len(ids):det+=np.asarray(xc[:,ids].sum(1)).ravel()>0;expr+=np.asarray(lx[:,ids].sum(1)).ravel()
  coverage.append(dict(sample=sm,panel=name,requested=len(genes),covered=sum(len(i)>0 for i in ix),missing=';'.join(g for g,i in zip(genes,ix) if not len(i))))
  for cl in o.cluster.astype(str).unique():
   mask=(o.cluster.astype(str)==cl).to_numpy();panel.append(dict(sample=sm,cluster=cl,current_label=str(o.loc[mask,'lineage_phase3'].iloc[0]),cells=int(mask.sum()),author_panel=name,pct_two_detected=float((det[mask]>=2).mean()*100),mean_log1p_program=float(expr[mask].mean()/len(genes)),comparison_level='published_marker_support_not_author_cell_labels'))
 o['pctUMI_OC_specific']=count(['OCSTAMP','DCSTAMP','ATP6V0D2','SIGLEC15','SUCNR1','CALCR'])/o.n_counts*100
 o['pctUMI_competing_mes']=o.pctUMI_Fib+o.pctUMI_Osteo+o.pctUMI_Mural
 ds=o.scDblFinder_score if o.scDblFinder_score.notna().all() else o.doublet_score
 base=(o.n_genes>=500)&(o.n_counts>=750)&ds.notna()&(ds<=ds.quantile(.9))&(o.pctUMI_competing_mes<=.2)&((o.pctUMI_Epi+o.pctUMI_Myogenic)<=.2)&(o.pctUMI_OC_specific<=.1)&~o.CTA_context
 o['ref_immune']=False
 for typ in ['T_NK','B_Plasma','Myeloid','DC']:
  o['ref_immune']|=(o.lineage_phase3==typ)&(o['det_'+typ+'_review']>=3)&(o['pctUMI_'+typ]>=.5)&base&(o.pctUMI_Endo<=.2)
 ipct=sum(o['pctUMI_'+k] for k in ['T_NK','B_Plasma','Myeloid','DC'])
 o['ref_endo']=(o.lineage_phase3=='Endo')&(o.det_Endo_review>=3)&(o.pctUMI_Endo>=.5)&base&(ipct<=.2)
 o['cnv_role']='target';rng=np.random.default_rng(270231)
 for typ,cap in [('immune',400),('endo',160)]:
  ids=np.flatnonzero(o['ref_'+typ].to_numpy());rng.shuffle(ids);ids=ids[:cap]
  if len(ids)<20:continue
  nh=max(10,int(len(ids)*.2));o.iloc[ids[:nh],o.columns.get_loc('cnv_role')]='holdout_'+typ;o.iloc[ids[nh:],o.columns.get_loc('cnv_role')]='ref_'+typ
 ref=o.cnv_role.str.startswith('ref_');hold=o.cnv_role.str.startswith('holdout_')
 # If just one small eligible lineage, increase holdout to 20 only when >=30 training remain.
 if hold.sum()<20 and ref.sum()+hold.sum()>=50:
  ix=np.flatnonzero(ref.to_numpy());rng.shuffle(ix)
  for j in ix[:20-int(hold.sum())]:o.iloc[j,o.columns.get_loc('cnv_role')]=str(o.iloc[j].cnv_role).replace('ref_','holdout_')
 ref=o.cnv_role.str.startswith('ref_');hold=o.cnv_role.str.startswith('holdout_')
 assert not (ref&hold).any() and not o.loc[ref|hold,'lineage_phase3'].isin(['FibLike','OsteoLike','Mural','MesAmb','Unresolved','OC_like','Myogenic_pending']).any()
 assert (o.loc[ref|hold,'pctUMI_competing_mes']<=.2).all()
 ids=[g if g.startswith('ENSG') else s2e.get(s,'UNMAPPED__'+g) for g,s in zip(a.var_names,symbols)];assert len(ids)==len(set(ids))
 expressed=np.asarray((x>0).sum(0)).ravel()>=3;autocov={str(ch):sum(bool(e) and g in go.index and go.at[g,'chr']==str(ch) for g,e in zip(ids,expressed)) for ch in range(1,23)}
 eligible=bool(ref.sum()>=30 and hold.sum()>=20 and min(autocov.values())>=30)
 row=dict(sample=sm,dataset=r['dataset'],site=r['site'],patient_id=r['patient_id'],run_cells=a.n_obs,targets=int((o.cnv_role=='target').sum()),ref_immune=int((o.cnv_role=='ref_immune').sum()),ref_endo=int((o.cnv_role=='ref_endo').sum()),holdout=int(hold.sum()),cnv_eligible=eligible,min_expressed_genes_per_autosome=min(autocov.values()),reason='PASS_pre_CNV_input_gate' if eligible else 'HOLD_insufficient_reference_holdout_or_gene_coverage',author_cell_labels_available=False)
 o['tool_id']=[sm+'_c'+str(i).zfill(6) for i in range(len(o))]
 o.to_csv(P/'review'/(sm+'_cells.tsv.gz'),sep='\t',index=False)
 if eligible:
  dest=P/'inputs'/sm;dest.mkdir(exist_ok=True)
  o.to_csv(dest/'cells.tsv',sep='\t',index=False)
  pd.DataFrame({'ensembl_id':ids,'symbol':symbols,'source_gene_id':a.var_names}).to_csv(dest/'genes.tsv',sep='\t',index=False)
  with gzip.open(dest/'counts.mtx.gz','wb',compresslevel=1) as f:mmwrite(f,x.T,field='integer')
  pd.DataFrame({'tool_id':o.tool_id,'annotation':np.where(ref,o.cnv_role,'C'+o.cluster.astype(str)+'_'+o.lineage_phase3)}).to_csv(dest/'annotations.tsv',sep='\t',index=False,header=False)
  (dest/'gene_coverage.json').write_text(json.dumps(autocov,indent=2));(dest/'INPUT_COMPLETE').write_text('CNV_INPUT_EXPORT_COMPLETE\n')
 summ.append(row);pd.DataFrame(summ).to_csv(P/'review/input_summary.tsv',sep='\t',index=False);print('INPUT_GATE',json.dumps(row),flush=True)
 del a,x,xc,lx,o
pd.DataFrame(panel).to_csv(P/'review/author_marker_support.tsv',sep='\t',index=False)
pd.DataFrame(transition).to_csv(P/'review/initial_vs_reviewed_labels.tsv',sep='\t',index=False)
pd.DataFrame(coverage).to_csv(P/'review/author_marker_coverage.tsv',sep='\t',index=False)
(P/'INPUTS_READY').write_text('GEO_CNV_INPUTS_READY\n');print('GEO_ANNOTATION_COMPARE_AND_CNV_INPUTS_READY',flush=True)
