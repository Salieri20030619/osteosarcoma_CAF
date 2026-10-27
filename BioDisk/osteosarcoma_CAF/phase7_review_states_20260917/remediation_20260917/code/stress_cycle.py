from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd,scanpy as sc,harmonypy as hm
from sklearn.metrics import adjusted_rand_score
from threadpoolctl import threadpool_limits
threadpool_limits(4)
p=Path.cwd();root=p/'remediation_20260917';out=root/'stress_cycle';out.mkdir(exist_ok=True)
a=sc.read_h5ad(p/'integration_review/strict_CAF_harmony_review.h5ad');main='harmony_theta2_r0.5'
st=pd.read_csv(root/'resources/stress_core_author.csv');cc=pd.read_csv(root/'resources/seurat_cc_updated2019.tsv',sep='\t');assert len(st)==512 and (st.logFC>0).sum()==507 and (st.logFC<0).sum()==5
up=[g for g in st.loc[st.logFC>0,'ensembl_gene_id'] if g in a.var_names];down=[g for g in st.loc[st.logFC<0,'ensembl_gene_id'] if g in a.var_names]
sc.tl.score_genes(a,up,score_name='stress507_up',random_state=20260917,use_raw=False);sc.tl.score_genes(a,down,score_name='stress5_down',random_state=20260917,use_raw=False);a.obs['stress512_signed']=a.obs.stress507_up-a.obs.stress5_down
s=a.var_names[a.var.symbol.isin(cc.loc[cc.phase=='S','gene'])].tolist();g=a.var_names[a.var.symbol.isin(cc.loc[cc.phase=='G2M','gene'])].tolist();sc.tl.score_genes_cell_cycle(a,s_genes=s,g2m_genes=g,random_state=20260917,use_raw=False)
a.obs.to_csv(out/'cell_scores.tsv.gz',sep='\t',index_label='cell_id')
a.obs.groupby(['patient_id',main],observed=True)[['stress507_up','stress5_down','stress512_signed','S_score','G2M_score']].median().to_csv(out/'patient_state_scores.tsv',sep='\t');pd.crosstab(a.obs[main],a.obs.phase).to_csv(out/'state_cell_cycle_counts.tsv',sep='\t')
params={'source_commit':'893f4d54c6577f2aeb0e85cb79e60bc17cf12dc1','stress_original':512,'up_original':507,'down_original':5,'up_present':len(up),'down_present':len(down),'S_original':43,'G2M_original':54,'S_present':len(s),'G2M_present':len(g),'cycle_source':'Seurat 5.4.0 cc.genes.updated.2019','score':'expression matched control; positive core minus negative core','seed':20260917}
params['source_SHA256']={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in [root/'resources/stress_core_author.csv',root/'resources/seurat_cc_updated2019.tsv']}
results=[]
def run(mode,keep,regress=False):
 b=a[keep].copy();excluded=set(up+down+s+g);mask=b.var.highly_variable.to_numpy()&~b.var_names.isin(excluded);b=b[:,mask].copy()
 if regress:sc.pp.regress_out(b,['stress512_signed','S_score','G2M_score'],n_jobs=1)
 sc.pp.scale(b,max_value=10);sc.tl.pca(b,n_comps=30,random_state=20260917,svd_solver='arpack');h=hm.run_harmony(b.obsm['X_pca'],b.obs,'patient_id',theta=2,random_state=20260917,ncores=4,max_iter_harmony=20)
 x=np.asarray(h.Z_corr);b.obsm['X_harmony']=x if x.shape[0]==b.n_obs else x.T;sc.pp.neighbors(b,use_rep='X_harmony',n_neighbors=15,random_state=20260917);sc.tl.leiden(b,key_added='sensitivity_state',resolution=.5,flavor='igraph',n_iterations=2,directed=False,random_state=20260917)
 b.obs[['patient_id',main,'sensitivity_state','phase','stress512_signed']].to_csv(out/(mode+'_cells.tsv.gz'),sep='\t',index_label='cell_id')
 original=b.obs[main].astype(str).to_numpy();new=b.obs.sensitivity_state.astype(str).to_numpy()
 for state in ['0','2']:
  t=original==state;matches=[]
  for c in set(new):
   q=new==c;inter=int((t&q).sum());matches.append((inter/int((t|q).sum()),c,inter/int(t.sum()),inter/int(q.sum())))
  j,c,rec,precision=max(matches);q=new==c;counts=pd.crosstab(b.obs.patient_id,new);eligible=((counts[c]>=20)&((counts.sum(1)-counts[c])>=20)).sum()
  results.append({'mode':mode,'cells':b.n_obs,'HVG':b.n_vars,'ARI':adjusted_rand_score(original,new),'target_state':state,'original_state_cells_retained':int(t.sum()),'best_match':c,'Jaccard':j,'recall':rec,'precision':precision,'match_paired_patients':int(eligible)})
 print(mode,'DONE',flush=True)
run('exclude_stress_cycle_genes',np.ones(a.n_obs,dtype=bool))
run('exclude_and_regress_scores',np.ones(a.n_obs,dtype=bool),True)
keep=(a.obs.phase=='G1')&(a.obs.stress512_signed<=a.obs.groupby('patient_id',observed=True).stress512_signed.transform(lambda x:x.quantile(.75)))
run('G1_lower_stress_sensitivity',keep.to_numpy())
pd.DataFrame(results).to_csv(out/'reclustering_sensitivity.tsv',sep='\t',index=False);(out/'parameters.json').write_text(json.dumps(params,indent=2));(out/'COMPLETE').write_text('STRESS_CYCLE_SENSITIVITY_COMPLETE\n');print(json.dumps(params));print(pd.DataFrame(results).to_string(index=False))
