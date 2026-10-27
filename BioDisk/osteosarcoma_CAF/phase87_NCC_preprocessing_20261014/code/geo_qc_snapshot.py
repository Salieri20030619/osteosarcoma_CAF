#!/usr/bin/env python3
"""SYNTAX_CHECK: python -m py_compile phase5.py
RUN_COMMAND: bash run_qc.sh
EXPECTED_SUCCESS_MARKER: GEO_QC_CLUSTERING_COMPLETE
Adapted from existing phase2.py; original outputs are preserved.
"""
import argparse,csv,gzip,hashlib,io,json,os,pathlib,platform,sys,tarfile,time,urllib.request
os.environ.setdefault('MPLBACKEND','Agg')
os.environ.setdefault('NUMBA_NUM_THREADS','4')
os.environ.setdefault('OMP_NUM_THREADS','4')
os.environ.setdefault('OPENBLAS_NUM_THREADS','4')
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmread
from threadpoolctl import threadpool_limits
threadpool_limits(limits=4)
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
from sklearn.metrics import adjusted_rand_score
from importlib.metadata import version
SEED=270231
SETS={
'T_NK':['CD3D','CD3E','TRAC','NKG7','KLRD1'],
'B_Plasma':['CD79A','MS4A1','CD79B','MZB1','JCHAIN'],
'Myeloid':['LYZ','LST1','TYROBP','AIF1','FCER1G'],
'OC':['ACP5','CTSK','ATP6V0D2','OCSTAMP','DCSTAMP'],
'Endo':['PECAM1','VWF','CLDN5','EMCN','KDR'],
'Epi':['EPCAM','KRT8','KRT18','KRT19','SFTPC','AGER'],
'Mural':['RGS5','CSPG4','NOTCH3','MCAM','MYH11'],
'FibLike':['DCN','LUM','PDGFRA','DPT','COL3A1','COL6A1'],
'OsteoLike':['RUNX2','SP7','ALPL','IBSP','BGLAP','DMP1']}
EXTRA=['FAP','PDPN','ACTA2','TAGLN','MYL9','TPM2','CNN1','SPP1','MGP','ENPP1','STC1','PTPRC','CD74','HLA-DRA','IL6','CXCL12','FOS','JUN']
DOT=['CD3D','NKG7','CD79A','MZB1','LYZ','LST1','ACP5','CTSK','PECAM1','VWF','EPCAM','SFTPC','RGS5','NOTCH3','DCN','LUM','PDGFRA','FAP','RUNX2','ALPL','IBSP','ACTA2','TAGLN']
def log(*s): print(time.strftime('%F %T'),*s,flush=True)
def dump(p,o): p.write_text(json.dumps(o,ensure_ascii=False,indent=2)+'\n')
def read_matrix(path):
 a=ad.read_h5ad(path/'counts.h5ad')
 assert sp.issparse(a.X) and a.obs_names.is_unique and a.var_names.is_unique
 return a,0

def mad(x):return float(np.median(np.abs(x-np.median(x))))
def stats(a):
 x=a.X;a.obs['n_counts']=np.asarray(x.sum(1)).ravel();a.obs['n_genes']=x.getnnz(1)
 mt=a.var.symbol.str.startswith('MT-').values
 assert mt.sum()>0,'Missing MT symbols: check feature annotation'
 a.obs['pct_mt']=np.divide(np.asarray(x[:,mt].sum(1)).ravel()*100,a.obs.n_counts,out=np.zeros(a.n_obs),where=a.obs.n_counts.values>0)
def get_gene(a,g,layer=None):
 ix=np.flatnonzero(a.var.symbol.values==g)
 x=a.X if layer is None else a.layers[layer]
 return np.asarray(x[:,ix].sum(1)).ravel() if len(ix) else np.zeros(a.n_obs)
def choose_label(row):
 gates={k:row['frac2_'+k]>=0.2 for k in SETS}
 # General myeloid genes also occur in osteoclasts; distinguish with osteoclast-specific markers.
 if gates['OC'] and row['frac2_OC']>=0.35:return 'OC_like'
 immune=[k for k in ['T_NK','B_Plasma','Myeloid'] if gates[k]]
 structural=[k for k in ['Endo','Epi','Mural','FibLike','OsteoLike'] if gates[k]]
 active=immune+structural
 if not active:return 'Unresolved'
 order=sorted(active,key=lambda k:row['score_'+k],reverse=True)
 best=order[0]
 if len(order)>1 and row['score_'+best]-row['score_'+order[1]]<0.15:return 'Mixed'
 if best in ['FibLike','OsteoLike','Mural']:
  if sum(gates[k] for k in ['FibLike','OsteoLike','Mural'])>1:return 'MesAmb'
  return best
 return best

def run_one(row,args):
 gsm=row['sample'];out=args.out/gsm;out.mkdir(parents=True,exist_ok=True)
 fingerprint=hashlib.sha256((pathlib.Path(__file__).read_text()+json.dumps(row,sort_keys=True)+str(args.doublet_rate)).encode()).hexdigest()
 done=out/'complete.json'
 if done.exists():
  old=json.loads(done.read_text())
  if old.get('fingerprint')==fingerprint and (out/'annotated.h5ad').exists():log('SAMPLE_CACHED',gsm);return old
  raise RuntimeError(f'{out}: existing result has different parameters; use a new --out directory')
 path=pathlib.Path(row['matrix_dir']);assert path.is_dir()
 log('SAMPLE_START',gsm)
 a,peaks=read_matrix(path);stats(a);n0=a.n_obs
 a.obs['barcode']=a.obs_names;a.obs_names=[gsm+'_'+s for s in a.obs_names]
 a.obs['sample']=gsm;a.obs['patient_id']=row['patient_id'];a.obs['site']=row['site'];a.obs['group']=row['group']
 for key in ['dataset','disease_stage_raw','treatment_raw','library_mapping_status','source_evidence']:a.obs[key]=row[key]
 a.obs['technology_reported']=row['technology'];a.obs['contains_peaks']=peaks>0
 for key in ['patient_independence_status','source_type','clinical_pairing','cross_cohort_patient_pooling']:a.obs[key]=row[key]
 a.obs['technology_review']='RNA_plus_Peaks_review' if peaks else 'RNA_only'
 a.obs['cnv_status']='not_run';a.obs['ambient_status']='not_corrected_filtered_matrix_only'
 lg=np.log10(a.obs.n_genes.values+1);lc=np.log10(a.obs.n_counts.values+1);mt=a.obs.pct_mt.values
 lowg=max(200,int(np.floor(10**(np.median(lg)-3*mad(lg))-1)))
 lowc=max(300,int(np.floor(10**(np.median(lc)-3*mad(lc))-1)))
 highmt=min(20.,max(5.,float(np.median(mt)+3*mad(mt))))
 a.obs['high_complexity_flag']=(lc>np.median(lc)+3*mad(lc))|(lg>np.median(lg)+3*mad(lg))
 a.obs['qc_low_genes']=a.obs.n_genes<lowg;a.obs['qc_low_counts']=a.obs.n_counts<lowc;a.obs['qc_high_mt']=a.obs.pct_mt>highmt
 a.obs['qc_pass_without_mt']=~(a.obs.qc_low_genes|a.obs.qc_low_counts)
 a.obs['qc_pass']=~(a.obs.qc_low_genes|a.obs.qc_low_counts|a.obs.qc_high_mt)
 allobs=a.obs.copy();q=a[a.obs.qc_pass].copy();assert q.n_obs>=30,('Too few cells after QC for clustering; preserve and review',gsm,q.n_obs)
 thresholds={'min_genes':lowg,'min_counts':lowc,'max_pct_mt':highmt,'MAD_scale':'unscaled median absolute deviation','lower_log_MAD':3,'upper_mt_MAD':3,'high_complexity':'flag_only','doublet_expected_rate':args.doublet_rate}
 dump(out/'qc_thresholds.json',thresholds)
 log('QC_PASS',gsm,q.n_obs,'thresholds',thresholds)
 threshold=float('nan');doublet_status='not_evaluated_low_cell_count'
 q.obs['doublet_score']=np.nan;q.obs['predicted_doublet']=pd.array([pd.NA]*q.n_obs,dtype='boolean')
 q.obs['doublet_exclude']=False
 if q.n_obs>=50:
  try:
   sc.pp.scrublet(q,expected_doublet_rate=args.doublet_rate,random_state=SEED,n_prin_comps=min(30,q.n_obs-2),use_approx_neighbors=False,verbose=False)
   threshold=float(q.uns['scrublet'].get('threshold',float('nan')))
   rate=float(q.obs.predicted_doublet.mean())
   valid=np.isfinite(threshold) and not q.obs.predicted_doublet.isna().any()
   doublet_status='automatic_threshold_pass' if valid and q.n_obs>=100 and 0.001<=rate<=0.25 else 'automatic_threshold_review_required'
   if doublet_status=='automatic_threshold_pass':q.obs['doublet_exclude']=q.obs.predicted_doublet.astype(bool)
  except (ValueError,RuntimeError,KeyError) as exc:
   doublet_status='scrublet_not_resolved_'+type(exc).__name__
   (out/'scrublet_issue.txt').write_text(repr(exc))
 q.obs['doublet_assessment_status']=doublet_status
 allobs['doublet_assessment_status']=doublet_status
 allobs['doublet_exclude']=q.obs.doublet_exclude.reindex(allobs.index).fillna(False).astype(bool)
 thresholds.update({'scrublet_threshold':threshold,'doublet_assessment_status':doublet_status,'threshold_review_policy':'review if <100 QC cells, missing threshold, predicted fraction <0.001 or >0.25; no fixed-fraction deletion','capture_mapping':row['library_mapping_status'],'doublet_prior_policy':'legacy 0.06 prior, exploratory sample-level assessment; not a measured loading rate','technology':row['technology'],'mt_policy':'sample-specific empirical thresholds; platform percentages not compared as biological effects'})
 dump(out/'qc_thresholds.json',thresholds)
 allobs['doublet_score']=q.obs.doublet_score.reindex(allobs.index)
 allobs['predicted_doublet']=q.obs.predicted_doublet.reindex(allobs.index).astype('boolean')
 allobs['analysis_keep']=False
 a=q[~q.obs.doublet_exclude.values].copy();assert a.n_obs>=30
 allobs.loc[a.obs_names,'analysis_keep']=True
 allobs.to_csv(out/'all_cells_qc.tsv.gz',sep='\t',index_label='cell_id')
 fig,axes=plt.subplots(1,3,figsize=(12,3.5))
 axes[0].scatter(allobs.n_counts,allobs.n_genes,c=allobs.analysis_keep.map({True:'#1764AB',False:'#B0B0B0'}),s=2,rasterized=True)
 axes[0].set(xscale='log',yscale='log',xlabel='UMI',ylabel='Genes',title=row['short'])
 axes[1].hist(allobs.pct_mt,bins=60,color='#1764AB');axes[1].axvline(highmt,color='#B22222');axes[1].set(xlabel='MT pct',ylabel='Cells')
 axes[2].hist(q.obs.doublet_score,bins=50,color='#1764AB');axes[2].axvline(threshold if np.isfinite(threshold) else 0,color='#B22222');axes[2].set(xlabel='DB score',ylabel='Cells')
 fig.tight_layout();fig.savefig(out/'QC.png',dpi=180);plt.close(fig)
 a.layers['counts']=a.X.copy()
 sc.pp.normalize_total(a,target_sum=1e4);sc.pp.log1p(a)
 # Score signatures on all genes, preserving original counts and Ensembl IDs.
 for name,genes in SETS.items():
  ids=a.var_names[a.var.symbol.isin(genes)].tolist();assert len(ids)>=2,(name,ids)
  sc.tl.score_genes(a,ids,score_name='score_'+name,random_state=SEED,ctrl_size=50,use_raw=False)
  a.obs['det_'+name]=sum((get_gene(a,g,'counts')>0).astype(int) for g in genes)
 sc.pp.highly_variable_genes(a,flavor='seurat',n_top_genes=3000)
 # Remove mitochondrial, ribosomal and sex-linked technical drivers from HVG embedding.
 exclude=a.var.symbol.str.match(r'^(MT-|RPL|RPS)').values | a.var.symbol.isin(['XIST','RPS4Y1','DDX3Y','EIF1AY']).values
 a.var.loc[exclude,'highly_variable']=False
 b=a[:,a.var.highly_variable].copy();sc.pp.scale(b,max_value=10)
 npc=min(40,b.n_obs-1,b.n_vars-1);assert npc>=10
 sc.tl.pca(b,n_comps=npc,svd_solver='arpack',random_state=SEED)
 a.obsm['X_pca']=b.obsm['X_pca'];del b
 sc.pp.neighbors(a,n_neighbors=15,n_pcs=npc,random_state=SEED)
 for res in [0.4,0.8,1.2]:sc.tl.leiden(a,resolution=res,key_added='leiden_'+str(res),random_state=SEED,flavor='igraph',n_iterations=2,directed=False)
 a.obs['cluster']=a.obs['leiden_0.8'].copy();sc.tl.umap(a,random_state=SEED)
 # Cluster marker ranking is descriptive, not a patient-level disease comparison.
 sizes=a.obs.cluster.value_counts();valid=sizes[sizes>=3].index.tolist()
 rows=[];labelmap={}
 for cluster in a.obs.cluster.cat.categories:
  o=a.obs[a.obs.cluster==cluster];rr={'sample':gsm,'cluster':cluster,'cells':len(o)}
  for name in SETS:
   rr['score_'+name]=float(o['score_'+name].mean());rr['frac2_'+name]=float((o['det_'+name]>=2).mean())
  rr['provisional_label']=choose_label(rr);labelmap[cluster]=rr['provisional_label'];rows.append(rr)
 a.obs['lineage_prelim']=a.obs.cluster.astype(str).map(labelmap).astype('category')
 mes=(a.obs.det_FibLike>=2)|(a.obs.det_OsteoLike>=2)|(a.obs.det_Mural>=2)
 a.obs['cnv_target_candidate']=mes
 immune=(a.obs.det_T_NK>=2)|(a.obs.det_B_Plasma>=2)|(a.obs.det_Myeloid>=2)
 a.obs['cnv_ref_candidate']=immune&~mes&(a.obs.det_OC<2)&(a.obs.det_Epi<2)&(a.obs.det_Endo<2)
 a.obs.to_csv(out/'cell_annotations.tsv.gz',sep='\t',index_label='cell_id')
 pd.DataFrame(rows).to_csv(out/'cluster_identity_evidence.tsv',sep='\t',index=False)
 # Export marker fractions and means for manual review, including all requested functional candidates.
 expr=[]
 for cl in a.obs.cluster.cat.categories:
  mask=(a.obs.cluster==cl).values
  for g in sorted(set(sum(SETS.values(),[])+EXTRA)):
   expr.append({'sample':gsm,'cluster':cl,'gene':g,'pct_detected':float((get_gene(a,g,'counts')[mask]>0).mean()*100),'mean_log1p':float(get_gene(a,g)[mask].mean())})
 pd.DataFrame(expr).to_csv(out/'marker_evidence.tsv',sep='\t',index=False)
 fig,axes=plt.subplots(1,2,figsize=(13,5))
 sc.pl.umap(a,color='cluster',ax=axes[0],show=False,title=row['short']+' clust',size=8,legend_loc='on data')
 sc.pl.umap(a,color='lineage_prelim',ax=axes[1],show=False,title=row['short']+' lineage',size=8,palette='tab20')
 fig.tight_layout();fig.savefig(out/'UMAP.png',dpi=180,bbox_inches='tight');plt.close(fig)
 e=pd.DataFrame(expr);clusters=list(a.obs.cluster.cat.categories)
 means=e.pivot(index='cluster',columns='gene',values='mean_log1p').reindex(index=clusters,columns=DOT)
 pct=e.pivot(index='cluster',columns='gene',values='pct_detected').reindex(index=clusters,columns=DOT)
 z=(means-means.mean())/means.std(ddof=0).replace(0,1)
 fig,ax=plt.subplots(figsize=(12,max(3,len(clusters)*.28+1.5)))
 xx,yy=np.meshgrid(np.arange(len(DOT)),np.arange(len(clusters)))
 pts=ax.scatter(xx.ravel(),yy.ravel(),s=pct.values.ravel()*1.2,c=z.values.ravel(),cmap='RdBu_r',vmin=-2,vmax=2,edgecolors='none')
 ax.set_xticks(range(len(DOT)),DOT,rotation=90);ax.set_yticks(range(len(clusters)),['C'+str(c)+' '+labelmap[c] for c in clusters]);ax.invert_yaxis();ax.set_title(row['short']+' markers');ax.set_xlabel('Gene');ax.set_ylabel('Cluster')
 fig.colorbar(pts,ax=ax,label='Mean z');fig.tight_layout();fig.savefig(out/'Markers.png',dpi=180,bbox_inches='tight');plt.close(fig)
 a.uns['phase2_thresholds']=thresholds;a.uns['annotation_status']='provisional_marker_supported_CNV_pending'
 a.write_h5ad(out/'annotated.h5ad',compression='gzip')
 summary={'sample':gsm,'short':row['short'],'dataset':row['dataset'],'patient_id':row['patient_id'],'site':row['site'],'clinical_group':row['group'],'input_cells':n0,'qc_pass':q.n_obs,'doublets':int(q.obs.doublet_exclude.sum()),'predicted_doublets':int(q.obs.predicted_doublet.sum()),'doublet_status':doublet_status,'retained':a.n_obs,'clusters':len(rows),'cnv_target_candidates':int(a.obs.cnv_target_candidate.sum()),'cnv_ref_candidates':int(a.obs.cnv_ref_candidate.sum()),'nonRNA_rows_removed':peaks,'scrublet_threshold':threshold,'ARI_04_08':adjusted_rand_score(a.obs['leiden_0.4'],a.obs['leiden_0.8']),'ARI_08_12':adjusted_rand_score(a.obs['leiden_0.8'],a.obs['leiden_1.2']),'fingerprint':fingerprint}
 dump(done,summary);log('SAMPLE_COMPLETE',json.dumps(summary));return summary

def main():
 p=argparse.ArgumentParser();p.add_argument('--manifest',type=pathlib.Path,required=True);p.add_argument('--raw',type=pathlib.Path,required=True);p.add_argument('--out',type=pathlib.Path,required=True);p.add_argument('--sample');p.add_argument('--doublet-rate',type=float,default=.06);args=p.parse_args()
 assert 0<args.doublet_rate<.5
 args.raw.mkdir(parents=True,exist_ok=True);args.out.mkdir(parents=True,exist_ok=True)
 sc.settings.n_jobs=4;sc.settings.verbosity=1;plt.rcParams.update({'font.family':'DejaVu Sans','pdf.fonttype':42})
 rows=list(csv.DictReader(args.manifest.open(),delimiter='\t'))
 if args.sample:rows=[r for r in rows if r['sample']==args.sample]
 assert rows,'No matching samples'
 dump(args.out/(('software_versions_'+args.sample+'.json') if args.sample else 'software_versions.json'),{'python':platform.python_version(),**{n:version(n) for n in ['scanpy','anndata','numpy','scipy','pandas','matplotlib','scikit-learn','igraph','leidenalg','numba']}})
 results=[run_one(row,args) for row in rows]
 if args.sample:return
 pd.DataFrame(results).drop(columns='fingerprint').to_csv(args.out/'sample_summary.tsv',sep='\t',index=False)
 tables=[pd.read_csv(args.out/r['sample']/'cluster_identity_evidence.tsv',sep='\t') for r in rows]
 combined=pd.concat(tables,ignore_index=True);combined.to_csv(args.out/'all_cluster_evidence.tsv',sep='\t',index=False)
 combined.pivot_table(index='sample',columns='provisional_label',values='cells',aggfunc='sum',fill_value=0).to_csv(args.out/'lineage_counts.tsv',sep='\t')
 (args.out/'SUCCESS.txt').write_text('GEO_QC_CLUSTERING_COMPLETE\n')
 log('GEO_QC_CLUSTERING_COMPLETE')
if __name__=='__main__':main()
