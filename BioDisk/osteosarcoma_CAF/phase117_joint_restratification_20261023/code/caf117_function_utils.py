from caf117_common import *
from scipy.stats import false_discovery_control
def endpoints(m,X,genes):
 s=pd.read_csv(B/'tables/cell_function_scores.tsv.gz',sep='\t');assert s.cell_id.equals(m.cell_id)
 cov=pd.read_csv(B/'tables/program_coverage.tsv',sep='\t');programs=cov.loc[cov.coverage.ge(.8),'program'].tolist();names=programs+[c for c in s if c.startswith('gene__')];family={n:('target_gene' if n.startswith('gene__') else 'functional_program') for n in names};v=s[names].to_numpy()
 p=pd.read_csv(R/'phase116_original33_analogs_20261023/tables/inherited_pathway_cell_scores.tsv.gz',sep='\t');assert p.cell_id.equals(m.cell_id);pn=[c for c in p if '::' in c];v=np.column_stack([v,p[pn].to_numpy()]);names+=pn;family.update({n:'inherited_pathway' for n in pn})
 dictionary=read('original33_dictionary.tsv').fillna('');gi={g:i for i,g in enumerate(genes)};rows=[]
 for _,r in dictionary.iterrows():
  gs=r.stable_markers.split(';') if r.stable_markers else r.display_markers.split(';');gs=[g for g in gs if g];meas=[g for g in gs if g in gi];n='signature__'+r.subgroup_id;valid=len(meas)>=3 and len(meas)/max(len(gs),1)>=.8
  vec=np.asarray(X[:,[gi[g] for g in meas]].mean(1)).ravel() if meas else np.full(len(m),np.nan)
  v=np.column_stack([v,vec if valid else np.full(len(m),np.nan)]);names.append(n);family[n]='original_signature';rows.append(dict(endpoint=n,defined=len(gs),measured=len(meas),valid_for_test=valid,genes=';'.join(meas)))
 save(pd.DataFrame(rows),'original_signature_coverage.tsv');return v,names,family
signcache={}
def pvalues(e,key):
 n=len(e);out=np.full(e.shape[1],np.nan)
 if n<2:return out
 valid=np.isfinite(e).all(0)&np.any(np.abs(e)>1e-12,axis=0)
 if not valid.any():return out
 if n<=16:
  if n not in signcache:signcache[n]=1-2*((np.arange(2**n,dtype=np.uint32)[:,None]>>np.arange(n))&1).astype(np.int8)
  signs=signcache[n];sim=abs(signs@e[:,valid]/n);out[valid]=(sim>=abs(e[:,valid].mean(0))[None,:]-1e-12).mean(0)
 else:
  seed=int(hashlib.sha256(str(key).encode()).hexdigest()[:8],16);signs=np.random.default_rng(seed).choice([-1,1],size=(65536,n)).astype(np.int8);sim=abs(signs@e[:,valid]/n);out[valid]=((sim>=abs(e[:,valid].mean(0))[None,:]-1e-12).sum(0)+1)/65537
 return out
def evaluate(m,v,names,family,label_column,arms,prefix):
 rows=[];effects=[]
 for src,q in m.groupby('source'):
  for arm in arms:
   base=q if arm=='all' else q[~q.doublet];minimum=5 if arm=='singlet_min5' else 1
   for cl in sorted(q[label_column].unique()):
    per={}
    for (donor,site),g in base.groupby(['donor','site'],dropna=False):
     a=g.index[g[label_column].eq(cl)].to_numpy();b=g.index[g[label_column].ne(cl)].to_numpy()
     if min(len(a),len(b))<minimum:continue
     delta=v[a].mean(0)-v[b].mean(0);per.setdefault(donor,[]).append((delta,len(a),len(b)))
    ds=sorted(per);e=np.array([np.mean([a[0] for a in per[d]],axis=0) for d in ds]) if ds else np.empty((0,len(names)));n=len(ds);ps=pvalues(e,(src,arm,cl,prefix))
    for i,d in enumerate(ds):
     for j,name in enumerate(names):effects.append(dict(source=src,arm=arm,joint_cluster=cl,donor=d,endpoint=name,effect=e[i,j],target_cells=sum(a[1] for a in per[d]),rest_cells=sum(a[2] for a in per[d]),sites=len(per[d])))
    for j,name in enumerate(names):
     ok=n>0 and np.isfinite(e[:,j]).all();rows.append(dict(source=src,arm=arm,joint_cluster=cl,endpoint=name,family=family[name],cells=int(base[label_column].eq(cl).sum()),source_patient_records=int(base.loc[base[label_column].eq(cl),'donor'].nunique()),N=n,mean_effect=float(e[:,j].mean()) if ok else np.nan,positive_patients=int((e[:,j]>1e-12).sum()) if ok else 0,P=ps[j],test='insufficient_N' if n<2 else 'coverage_missing' if not ok else 'all_zero_not_tested' if np.all(abs(e[:,j])<1e-12) else 'exact_signflip' if n<=16 else 'MC65536_plus_observed'))
  event(prefix,'source_statistics_complete',source=src)
 df=pd.DataFrame(rows);df['FDR']=np.nan
 for _,idx in df.groupby(['source','arm','family']).groups.items():
  ii=df.loc[idx].index[df.loc[idx].P.notna()];df.loc[ii,'FDR']=false_discovery_control(df.loc[ii,'P'].to_numpy())
 df['nominal']=df.P.lt(.05);df['FDR_significant']=df.FDR.lt(.05);save(df,prefix+'_statistics.tsv');save(pd.DataFrame(effects),prefix+'_patient_effects.tsv.gz');return df
