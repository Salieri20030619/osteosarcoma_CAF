from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
import json,datetime,gc
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
from scipy.stats import rankdata,binomtest
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase23_loose_first_20260921'
assert (O/'MODEL_COMPLETE').exists()
def read(n):return pd.read_csv(O/'tables'/n,sep='\t',low_memory=False)
def save(d,n):d.to_csv(O/'tables'/n,sep='\t',index=False)
def truth(s):return s.fillna(False).astype(str).str.lower().isin(['true','1','1.0'])
d=read('new_all_entry_labels.tsv.gz').set_index('cell_id');d['cluster']=d.cluster.astype(str);cls=sorted(d.cluster.unique(),key=int)
genes=read('consensus_HVG_rank.tsv').ensembl_id.astype(str).tolist();symbol={};ncl=len(cls);sums=np.zeros((ncl,len(genes)),np.float64);detected=np.zeros_like(sums);counts=np.zeros(ncl,np.int64)
for sm in d['sample'].unique():
 b=ad.read_h5ad(O/'inputs'/f'{sm}.h5ad',backed='r');ids=d.index[d['sample']==sm];ob=d.loc[ids];idx=b.obs_names.get_indexer(ids);assert (idx>=0).all()
 gm=pd.Series(np.arange(b.n_vars),index=b.var.ensembl_id.astype(str));ix=gm.loc[genes].to_numpy();x=b.X[idx].tocsr()[:,np.sort(ix)][:,np.argsort(np.argsort(ix))].astype(np.float64)
 nc=ob.n_counts.to_numpy();norm=x.multiply((1e4/nc)[:,None]).tocsr();norm.data=np.log1p(norm.data);calls=pd.Categorical(ob.cluster,categories=cls).codes
 G=sparse.csr_matrix((np.ones(len(calls)),(calls,np.arange(len(calls)))),shape=(ncl,len(calls)))
 sums+=(G@norm).toarray();binary=x.copy();binary.data[:]=1;detected+=(G@binary).toarray();counts+=np.bincount(calls,minlength=ncl)
 symbol.update(dict(zip(b.var.ensembl_id,b.var.symbol)));b.file.close();del b,x,norm,binary;gc.collect();print('MARKERS_SAMPLE',sm,flush=True)
rows=[];total=sums.sum(0);totaldet=detected.sum(0);N=counts.sum()
for k,cl in enumerate(cls):
 inside=sums[k]/counts[k];outside=(total-sums[k])/(N-counts[k]);delta=inside-outside;frac=detected[k]/counts[k];outfrac=(totaldet-detected[k])/(N-counts[k])
 for j in np.argsort(delta)[-50:][::-1]:rows.append({'cluster':cl,'gene_id':genes[j],'symbol':symbol.get(genes[j],''),'mean_logCP10k':inside[j],'delta_logCP10k':delta[j],'fraction_detected':frac[j],'outside_fraction':outfrac[j],'evidence':'descriptive; no cell-level P value'})
save(pd.DataFrame(rows),'all_gene_cluster_top50.tsv')
# An ordered attrition audit, not filters for the new entry model.
fib=(d.det_Fib>=2)|((d.det_Fib>=1)&(d.det_ECM>=3));comp=np.zeros(len(d),bool)
for name,ng in [('T_NK',3),('B_Plasma',3),('Myeloid',3),('DC',3),('Endo',3),('Epi',3),('Myogenic',3),('OC_specific',2)]:comp|=((d['det_'+name]>=ng)&(d['pct_'+name]>1)).to_numpy()
lo=d.historical_infer_call.eq('low_signal')|d.historical_cnv2_call.eq('diploid');hi=d.historical_infer_call.eq('elevated_signal')|d.historical_cnv2_call.eq('aneuploid')
rules=[('0_provider_called',np.ones(len(d),bool)),('1_prior_project_QC',truth(d.historical_QC_keep)),('2_prior_300_genes',d.n_genes>=300),('3_prior_Fib_ECM_marker_gate',fib),('4_prior_competing_lineage_gate',~comp),('5_prior_CNV_low_no_high',lo&~hi),('6_prior_reference_acceptance',truth(d.historical_reference_QC_pass)),('7_prior_strict_membership',truth(d.historical_strict))]
cum=np.ones(len(d),bool);audit=[];variants={}
for name,condition in rules:
 v=np.asarray(condition,bool);before=cum.copy();cum&=v;variants[name]=cum.copy()
 for sm,g in d.groupby('sample',sort=False):
  ix=d.index.get_indexer(g.index);audit.append({'sample':sm,'dataset':g.dataset.iloc[0],'step':name,'entry_cells':len(ix),'standalone_condition_pass':int(v[ix].sum()),'cumulative_retained':int(cum[ix].sum()),'incremental_removed':int((before[ix]&~cum[ix]).sum())})
save(pd.DataFrame(audit),'ordered_legacy_gate_attrition.tsv')
gatecols=pd.DataFrame({k:v for k,v in variants.items()},index=d.index);gatecols.to_csv(O/'tables/ordered_gate_cell_flags.tsv.gz',sep='\t',index_label='cell_id')
COHORTS=['GSE152048','GSE162454','GSE318840'];PROGS=['myCAF_customer','iCAF_customer','apCAF_customer','Osteogenic_customer','Mineral_customer']
def exact(v):
 z=np.round(np.asarray(v,float),12);z=z[z!=0]
 if not len(z):return 1.
 r=np.rint(rankdata(abs(z))*2).astype(int);dist=np.zeros(int(r.sum())+1,dtype=np.int64);dist[0]=1;last=0
 for rr in r:
  old=dist[:last+1].copy();dist[rr:rr+last+1]+=old;last+=rr
 assert int(dist.sum())==2**len(z)
 t=min(int(r[z>0].sum()),int(r[z<0].sum()));return min(1.,2*int(dist[:t+1].sum())/2**len(z))
assert exact([1,1,1,1])==.125 and exact([1,-1])==1
def stats(e,keys):
 rows=[]
 for key in keys:
  g=e
  for c,v in key.items():g=g[g[c]==v]
  v=g.delta.to_numpy(float);n=len(v);pos=int((v>1e-12).sum());neg=int((v< -1e-12).sum())
  rows.append({**key,'n_patients':n,'n_datasets':g.dataset.nunique(),'positive':pos,'negative':neg,'median_delta':float(np.median(v)) if n else np.nan,'P':exact(v) if n>=2 else np.nan,'P_sign':binomtest(pos,pos+neg).pvalue if n>=2 and pos+neg else np.nan,'min_inside_cells':g.inside_cells.min() if n else 0,'scope':'postselection_exploratory'})
 out=pd.DataFrame(rows);out['FDR']=np.nan
 for _,ix in out.groupby(['comparison_scope','minimum_cells']).groups.items():
  pv=out.loc[ix,'P'];out.loc[ix,'FDR']=np.where(pv.notna(),multipletests(pv.fillna(1),method='fdr_bh')[1],np.nan)
 return out
eligible=d[d.dataset.isin(COHORTS)&truth(d.known_patient)].copy();assert not eligible.donor_record.str.contains('UNKNOWN|UNCONFIRMED',case=False).any();effects=[];keys=[]
scopes={'all_entry_clusters':np.ones(len(d),bool)}
reviewfile=O/'tables/reviewed_cluster_lineage.tsv'
if reviewfile.exists():
 ann=read('reviewed_cluster_lineage.tsv');selected=ann.loc[ann.CAF_exploration.astype(str).str.lower().eq('true'),'cluster'].astype(str);scopes['within_mesenchymal_exploration']=d.cluster.isin(selected).to_numpy()
for scope,mask in scopes.items():
 dd=eligible[eligible.index.isin(d.index[mask])];cc=sorted(dd.cluster.unique(),key=int)
 for cl in cc:
  for donor,g in dd.groupby('donor_record',observed=True):
   inside=g[g.cluster==cl];outside=g[g.cluster!=cl]
   for pr in PROGS:
    # No 80% panel gate: observed coverage is recorded; scores use measured genes only.
    x=pd.to_numeric(inside['score_'+pr],errors='coerce').dropna();y=pd.to_numeric(outside['score_'+pr],errors='coerce').dropna()
    if len(x) and len(y):effects.append({'comparison_scope':scope,'cluster':cl,'program':pr,'donor_record':donor,'dataset':g.dataset.iloc[0],'inside_cells':len(x),'outside_cells':len(y),'delta':x.mean()-y.mean(),'panel_coverage':g['coverage_'+pr].min()})
  for minimum in [1,10,20]:
   for pr in PROGS:keys.append({'comparison_scope':scope,'minimum_cells':minimum,'cluster':cl,'program':pr})
E=pd.DataFrame(effects);save(E,'patient_program_differences_no_entry_gate.tsv');expanded=[]
for minimum in [1,10,20]:
 ee=E[(E.inside_cells>=minimum)&(E.outside_cells>=minimum)].copy();ee['minimum_cells']=minimum;expanded.append(ee)
tests=stats(pd.concat(expanded),keys);save(tests,'patient_program_tests_no_entry_gate.tsv')
save(eligible[['sample','dataset','patient_id','donor_record']].drop_duplicates(),'primary_patient_roster.tsv')
ag=pd.DataFrame(audit).groupby('step',sort=False).cumulative_retained.sum();fig,ax=plt.subplots(figsize=(10,5));ax.barh(ag.index,ag.to_numpy(),color='#337f9b');ax.invert_yaxis();ax.set(xlabel='Cells retained when historical gates are added in order',title='Entry model uses all cells; historical conditions shown only as sensitivity');fig.tight_layout();fig.savefig(O/'figures/ordered_gate_attrition.png',dpi=180);plt.close(fig)
qa={'all_entry_cells':len(d),'clusters':len(cls),'primary_independent_patients':eligible.donor_record.nunique(),'program_tests':len(tests),'no_entry_gene_UMI_MT_marker_CNV_reference_gates':True,'no_cell_level_inference':True,'old_C7_C9_not_used_to_fit_clusters':True,'minimum_10_and_20_cells_are_sensitivity_only':True,'source_files_preserved':True}
(O/'delivery/review_QA.json').write_text(json.dumps(qa,indent=2));(O/'REVIEW_COMPLETE').write_text('ALL_CLUSTER_MARKERS_AND_PATIENT_TESTS_COMPLETE\n');print('REVIEW_COMPLETE',qa,flush=True)
