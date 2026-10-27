from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='2'
import numpy as np,pandas as pd,json,hashlib
from scipy.stats import false_discovery_control
R=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=R/'phase115_incremental_function_20261023'
def dump(n,x):(O/n).write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str)+'\n',encoding='utf-8')
def save(x,n):x.to_csv(O/'tables'/n,sep='\t',index=False)
m=pd.read_csv(O/'tables/cell_mapping.tsv.gz',sep='\t',low_memory=False);s=pd.read_csv(O/'tables/cell_function_scores.tsv.gz',sep='\t');assert m.cell_id.equals(s.cell_id)
cluster_file=O/'tables/query_cluster_membership.tsv.gz'
if cluster_file.exists():
 c=pd.read_csv(cluster_file,sep='\t').set_index('cell_id');m['cluster_compatible']=m.cell_id.map(c.cluster_compatible).fillna('not_applicable')
programs=list(json.loads((O/'resources/PROGRAMS.json').read_text()));cols=programs+[c for c in s if c.startswith('gene__')]
coverage=pd.read_csv(O/'tables/program_coverage.tsv',sep='\t').set_index('program');validprog=set(coverage[coverage.coverage.ge(.8)].index)
values=s[cols].to_numpy(float);rng=np.random.default_rng(20261023);signs={}
def pvals(e):
 n=len(e);out=np.full(e.shape[1],np.nan)
 if n<2:return out
 if n not in signs:
  if n<=16:signs[n]=1-2*((np.arange(2**n,dtype=np.uint32)[:,None]>>np.arange(n))&1).astype(np.int8)
  else:signs[n]=rng.choice([-1,1],size=(65536,n)).astype(np.int8)
 sig=signs[n];test=np.abs(sig@e/n);obs=np.abs(e.mean(0));ext=(test>=obs[None,:]-1e-12).sum(0)
 out=ext/len(sig) if n<=16 else (ext+1)/(len(sig)+1)
 out[np.all(np.abs(e)<1e-12,axis=0)]=np.nan;return out
rows=[];effects=[];counts=[]
for src,q in m.groupby('source',sort=True):
 print('SOURCE',src,len(q),flush=True)
 for arm in ['all','singlet','singlet_min5']:
  base=q if arm=='all' else q[~q.doublet];minimum=5 if arm.endswith('min5') else 1
  levels=[('class','class_compatible'),('state','state_compatible')]+([('cluster','cluster_compatible')] if src in ['NCC','NCI'] and cluster_file.exists() else [])
  for level,col in levels:
   labels=['F1','F2','F3','F4','F5','Unresolved'] if level=='class' else sorted(m[m.cohort.eq('baseline')].state_label.unique())+['Unresolved']
   if level=='cluster':labels=sorted(q.cluster_compatible.unique())
   for label in labels:
    islabel=base[col].eq(label)
    cells=int(islabel.sum());donors=base.loc[islabel,'donor'].nunique();counts.append(dict(source=src,arm=arm,level=level,label=label,cells=cells,source_patients=donors))
    per=[]
    for (d,site),g in base.groupby(['donor','site'],dropna=False):
     a=g.index[g[col].eq(label)].to_numpy();b=g.index[g[col].ne(label)].to_numpy()
     if len(a)<minimum or len(b)<minimum:continue
     delta=values[a].mean(0)-values[b].mean(0)
     per.append((d,site,delta,len(a),len(b)))
    donorrows=[]
    for d in sorted(set(x[0] for x in per)):
     z=[x for x in per if x[0]==d];dd=np.mean([x[2] for x in z],axis=0);donorrows.append(dd)
     for j,c in enumerate(cols):effects.append(dict(source=src,arm=arm,level=level,label=label,donor=d,endpoint=c,effect=dd[j],sites=len(z),target_cells=sum(x[3] for x in z),rest_cells=sum(x[4] for x in z)))
    e=np.array(donorrows) if donorrows else np.empty((0,len(cols)));n=len(e);ps=pvals(e) if n else np.full(len(cols),np.nan)
    for j,c in enumerate(cols):
     valid=c.startswith('gene__') or c in validprog
     allzero=n>0 and np.all(np.abs(e[:,j])<1e-12)
     rows.append(dict(source=src,cohort='added' if src in ['NCC','NCI'] else 'baseline',arm=arm,level=level,label=label,endpoint=c,family_kind='gene' if c.startswith('gene__') else 'program',cells=cells,source_patients=donors,N=n,mean_effect=float(e[:,j].mean()) if n else np.nan,median_effect=float(np.median(e[:,j])) if n else np.nan,positive_patients=int((e[:,j]>1e-12).sum()) if n else 0,P=float(ps[j]) if valid else np.nan,test='not_tested_N_lt2' if n<2 else 'not_tested_low_panel_coverage' if not valid else 'all_zero_effect_not_tested' if allzero else 'exact_patient_signflip' if n<=16 else 'MC65536_plus_observed',minimum_cells_per_side=minimum))
df=pd.DataFrame(rows);df['FDR']=np.nan
for key,idx in df.groupby(['source','arm','level','family_kind']).groups.items():
 ii=df.loc[idx].index[df.loc[idx].P.notna()];df.loc[ii,'FDR']=false_discovery_control(df.loc[ii,'P'].to_numpy(),method='bh')
df['nominal']=df.P.lt(.05);df['FDR_significant']=df.FDR.lt(.05);save(df,'functional_statistics.tsv');save(pd.DataFrame(effects),'patient_function_effects.tsv.gz');save(pd.DataFrame(counts),'label_counts_by_arm.tsv')
# Patient-level class profiles describe the full set of predeclared programs, including negative results.
joined=pd.concat([m[['cell_id','source','donor','site','cohort','class_compatible','state_compatible','doublet']],s[cols+[c for c in s if c.startswith('lineage__')]]],axis=1)
profile=joined[~joined.doublet].groupby(['source','donor','class_compatible'])[cols].mean().reset_index();save(profile,'patient_class_function_profiles.tsv')
# Descriptive functional resemblance, separate from transcriptome mapping and never used to claim significance.
axes=['myCAF_customer','iCAF_customer','apCAF_customer','fibrillar_I_III_V','collagen_IV','collagen_crosslinking','ECM_degradation','Osteogenic_customer','cartilage_matrix','stress_descriptive','cycle_descriptive']
axes=[a for a in axes if a in validprog]
z=joined[(joined.cohort=='baseline')&~joined.doublet];p=z.groupby(['source','donor','state_compatible'])[axes].mean().reset_index();templates=p.groupby(['source','state_compatible'])[axes].mean().groupby('state_compatible').mean()
scale=np.maximum(p[axes].std().to_numpy(),.05);center=p[axes].mean().to_numpy();T=(templates.to_numpy()-center)/scale
T=T-T.mean(1,keepdims=True);T=T/np.maximum(np.linalg.norm(T,axis=1,keepdims=True),1e-12)
new=joined[joined.cohort.eq('added')].copy();v=(new[axes].to_numpy()-center)/scale;v=v-v.mean(1,keepdims=True);v=v/np.maximum(np.linalg.norm(v,axis=1,keepdims=True),1e-12);sim=v@T.T;ii=np.argmax(sim,axis=1)
new['function_nearest_state']=templates.index.to_numpy()[ii];new['function_similarity']=sim[np.arange(len(new)),ii];new['F2C3_function_similarity']=sim[:,templates.index.get_loc('F2_C3')]
save(new[['cell_id','source','donor','class_compatible','state_compatible','doublet','function_nearest_state','function_similarity','F2C3_function_similarity']+axes],'descriptive_function_resemblance.tsv.gz')
save(new.groupby(['source','function_nearest_state']).agg(cells=('cell_id','size'),source_patients=('donor','nunique'),median_similarity=('function_similarity','median')).reset_index(),'descriptive_function_resemblance_counts.tsv')
save(templates.reset_index(),'reference_state_function_templates.tsv')
# Whether unmatched cells preferentially carry a function is tested under the same patient-level framework above, not silently discarded.
key=df[df.cohort.eq('added')&df.arm.eq('singlet')&df.nominal&df.mean_effect.gt(0)];save(key,'added_nominal_positive_results.tsv')
qa=dict(status='PASS',statistics=len(df),valid_tests=int(df.P.notna().sum()),added_valid_tests=int(df[df.cohort.eq('added')].P.notna().sum()),added_nominal_positive_singlet=len(key),effect_unit='Mean cell log1pCP10k difference, equal site then donor weights',no_cell_level_inference=True,unassessed_P_missing=True,FDR_families='source/arm/level/program_or_gene',low_panel_coverage_excluded_from_tests=sorted(set(programs)-validprog),baseline_statistics_new_comparable_raw_common_gene_layer='Original corrected Phase69 statistics are preserved and not overwritten')
dump('STATISTICS_COMPLETE.json',qa);print(json.dumps(qa));print(key[['source','level','label','endpoint','N','mean_effect','positive_patients','P','FDR']].head(55).to_string(index=False))
