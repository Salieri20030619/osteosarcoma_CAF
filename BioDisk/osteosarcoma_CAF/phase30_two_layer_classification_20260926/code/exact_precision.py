from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='4'
import json,datetime,itertools,shutil,hashlib
import numpy as np,pandas as pd,anndata as ad
from statsmodels.stats.multitest import multipletests
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase30_two_layer_classification_20260926'
back=O/'precision_audit_before_exact';back.mkdir(exist_ok=False)
file=O/'tables/expanded_SCPCP000017_all_gene_statistics.tsv.gz';gf=pd.read_csv(file,sep='\t')
pf=pd.read_csv(O/'tables/class_program_statistics.tsv',sep='\t')
for f in [file,O/'tables/class_program_statistics.tsv',O/'tables/class_surface_gene_statistics.tsv',O/'tables/class_top_marker_statistics.tsv',O/'tables/statistics_summary.tsv',O/'resources/CLASSIFICATION_AND_STATISTICS_PROTOCOL.json']:
 shutil.copy2(f,back/f.name)
a=ad.read_h5ad(O/'objects/expanded_integrated.h5ad')
labels=pd.read_csv(O/'tables/expanded_classified_cells.tsv.gz',sep='\t',low_memory=False).set_index('cell_id');a.obs['class_id']=labels.class_id.reindex(a.obs_names)
panels=[c for c in a.obs if c.startswith('program_')]
def exact(d):
 n,g=d.shape;signs=np.array(list(itertools.product([-1,1],repeat=n)),dtype=np.float64);p=np.empty(g);obs=abs(d.mean(0))
 for start in range(0,g,256):
  null=np.abs(signs@d[:,start:start+256]/n);p[start:start+256]=(null>=obs[start:start+256][None,:]-1e-12).mean(0)
 return p
audit=[]
for cls in sorted(gf.loc[gf.test.eq('paired_signflip_MC20000'),'class_id'].unique()):
 gene=[];prog=[]
 for dn,g in a.obs[a.obs.dataset.eq('SCPCP000017')].groupby('donor_record',observed=True):
  gd=[];pdiff=[]
  for site,z in g.groupby('site',observed=True):
   target=z.class_id.eq(cls);other=~target
   if not target.any() or not other.any():continue
   it=a.obs.index.get_indexer(z.index[target]);ir=a.obs.index.get_indexer(z.index[other])
   tc=np.asarray(a.layers['counts'][it].sum(0)).ravel();oc=np.asarray(a.layers['counts'][ir].sum(0)).ravel()
   td=float(a.obs.iloc[it].full_library_counts.sum());od=float(a.obs.iloc[ir].full_library_counts.sum())
   gd.append(np.log2(tc/td*1e4+.1)-np.log2(oc/od*1e4+.1));pdiff.append(z.loc[target,panels].mean().to_numpy()-z.loc[other,panels].mean().to_numpy())
  if gd:gene.append(np.mean(gd,axis=0));prog.append(np.mean(pdiff,axis=0))
 gd=np.vstack(gene).astype(float);pdiff=np.vstack(prog).astype(float);assert len(gd)<=15
 idx=gf.class_id.eq(cls);assert gf.loc[idx,'feature'].tolist()==list(a.var_names)
 assert np.allclose(gf.loc[idx,'median_delta'],np.median(gd,axis=0))
 gp=exact(gd);pp=exact(pdiff)
 oldp=gf.loc[idx,'P'].to_numpy();gf.loc[idx,'P']=gp;gf.loc[idx,'test']='exact_paired_signflip'
 for j,name in enumerate(panels):
  ix=pf.layer.eq('expanded')&pf.pool.eq('SCPCP000017')&pf.class_id.eq(cls)&pf.feature.eq(name);assert ix.sum()==1
  pf.loc[ix,'P']=pp[j];pf.loc[ix,'test']='exact_paired_signflip'
 audit.append({'class_id':cls,'n_donors':len(gd),'sign_assignments':2**len(gd),'max_abs_gene_P_change':float(abs(gp-oldp).max())})
 print('EXACT_COMPLETED',audit[-1],flush=True)
def fdr(d):
 for _,g in d.groupby('family'):
  q=multipletests(g.P.fillna(1),method='fdr_bh')[1];d.loc[g.index,'FDR']=np.where(g.P.notna(),q,np.nan)
 d['nominal_exploratory']=d.P.notna()&(d.P<.05)&(d.median_delta>0)
 return d
gf=fdr(gf);pf=fdr(pf);gf.to_csv(file,sep='\t',index=False);pf.to_csv(O/'tables/class_program_statistics.tsv',sep='\t',index=False)
surfaces=[];top=[];summary=[]
for f in sorted((O/'tables').glob('*_all_gene_statistics.tsv.gz')):
 d=pd.read_csv(f,sep='\t');surfaces.append(d[d.predefined_surface_candidate])
 for cls,g in d.groupby('class_id'):
  q=g[g.median_delta>0].sort_values(['nominal_exploratory','P','median_delta'],ascending=[False,True,False])
  top.append(q.head(30))
 layer=d.layer.iloc[0];pool=d.pool.iloc[0];pr=pf[pf.layer.eq(layer)&pf.pool.eq(pool)]
 summary.append({'layer':layer,'pool':pool,'gene_hypotheses':len(d),'gene_tests_with_replication':int(d.P.notna().sum()),'nominal_positive_gene_markers':int(d.nominal_exploratory.sum()),'FDR_positive_gene_markers':int(((d.FDR<.05)&(d.median_delta>0)).sum()),'program_tests_with_replication':int(pr.P.notna().sum()),'nominal_positive_programs':int(pr.nominal_exploratory.sum()),'min_gene_P':float(d.P.min()) if d.P.notna().any() else None,'max_independent_donors':int(d.n_independent_donors.max())})
pd.concat(surfaces).to_csv(O/'tables/class_surface_gene_statistics.tsv',sep='\t',index=False);pd.concat(top).to_csv(O/'tables/class_top_marker_statistics.tsv',sep='\t',index=False);pd.DataFrame(summary).to_csv(O/'tables/statistics_summary.tsv',sep='\t',index=False)
protocol=json.loads((O/'resources/CLASSIFICATION_AND_STATISTICS_PROTOCOL.json').read_text())
protocol['numerical_precision_update']='All tests with 13-15 donors upgraded from Monte Carlo to full sign enumeration, for all genes/programs in the affected families. Same statistic, donors, classes and hypotheses. No model/selection changes; original estimates preserved.'
protocol['test_statistic']='Two-sided paired sign-flip mean, exact for <=15 donors, 19999 random flips+observed beyond15'
(O/'resources/CLASSIFICATION_AND_STATISTICS_PROTOCOL.json').write_text(json.dumps(protocol,ensure_ascii=False,indent=2))
qa={'completed_at':datetime.datetime.now().astimezone().isoformat(),'all_MC_families_upgraded':True,'class_selection_unchanged':True,'statistic_unchanged':True,'affected':audit}
(O/'delivery/EXACT_PRECISION_QA.json').write_text(json.dumps(qa,indent=2))
print('NOMINAL_PROGRAMS',pf[pf.nominal_exploratory][['layer','pool','class_id','feature','n_independent_donors','median_delta','P','FDR']].to_csv(sep='\t',index=False))
