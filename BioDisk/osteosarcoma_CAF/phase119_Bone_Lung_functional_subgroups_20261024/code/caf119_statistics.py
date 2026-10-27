from caf119_common import *
from scipy.stats import false_discovery_control
SIGN_CACHE={}
def summarize_effects(e,key):
 p=np.full(e.shape[1],np.nan);n=np.isfinite(e).sum(0);mu=np.full(e.shape[1],np.nan);pos=np.zeros(e.shape[1],int);groups={}
 for j in range(e.shape[1]):
  mask=np.isfinite(e[:,j]);v=e[mask,j]
  if len(v):mu[j]=v.mean();pos[j]=int((v>1e-12).sum())
  if len(v)<2 or np.all(abs(v)<1e-12):continue
  groups.setdefault(mask.tobytes(),(mask,[]))[1].append(j)
 for mask,cols in groups.values():
  z=e[mask][:,cols];nn=len(z)
  if nn<=16:
   if nn not in SIGN_CACHE:SIGN_CACHE[nn]=1-2*((np.arange(2**nn,dtype=np.uint32)[:,None]>>np.arange(nn))&1).astype(np.int8)
   a=SIGN_CACHE[nn];p[cols]=(abs(a@z)>=abs(z.sum(0))[None,:]-1e-10).mean(0)
  else:
   seed=int(hashlib.sha256((str(key)+str(cols)).encode()).hexdigest()[:8],16);a=np.random.default_rng(seed).choice([-1,1],size=(65536,nn)).astype(np.int8);p[cols]=((abs(a@z)>=abs(z.sum(0))[None,:]-1e-10).sum(0)+1)/65537
 return n,mu,pos,p
def evaluate(m,labels,nodes,score,prefix,holdout=False):
 assert score.cell_id.equals(m.cell_id);terms=[c for c in score if '__' in c];values=score[terms].to_numpy();rows=[];effects=[]
 for node in nodes.itertuples():
  background_options=['within_site','fibro_context'] if node.level=='parent' else ['within_parent','within_site','fibro_context']
  for src in m.loc[m.site.eq(node.site),'source'].unique():
   for scope in background_options:
    arms=['singlet'] if holdout else ['all','singlet','singlet_min5']
    for arm in arms:
     if scope=='fibro_context' and arm=='all':continue
     base=m[m.site.eq(node.site)&m.source.eq(src)]
     if arm!='all':base=base[~base.doublet]
     if holdout:base=base[base.holdout]
     if scope=='within_parent':base=base[labels.loc[base.index,'parent'].eq(node.parent)]
     if scope=='fibro_context':base=base[base.fibro_supported]
     col='parent' if node.level=='parent' else 'leaf';minimum=5 if arm=='singlet_min5' else 1;per=[];ds=[];tn=[];rn=[]
     for donor,g in base.groupby('donor'):
      target=labels.loc[g.index,col].eq(node.node_id);a=g.index[target].to_numpy();b=g.index[~target].to_numpy()
      if min(len(a),len(b))<minimum:continue
      per.append(values[a].mean(0)-values[b].mean(0));ds.append(donor);tn.append(len(a));rn.append(len(b))
     e=np.asarray(per) if per else np.empty((0,len(terms)));n,mu,pos,p=summarize_effects(e,(prefix,node.node_id,src,scope,arm))
     for j,term in enumerate(terms):
      family=term.split('__')[0];rows.append(dict(site=node.site,source=src,node_id=node.node_id,parent=node.parent,level=node.level,background=scope,arm=arm,endpoint=term,family=family,cells=int(labels.loc[base.index,col].eq(node.node_id).sum()),eligible_patient_pairs=len(ds),N=int(n[j]),mean_effect=mu[j],positive_patients=int(pos[j]),P=p[j],test='not_tested_N_lt2' if n[j]<2 else 'not_tested_allzero' if np.all(abs(e[np.isfinite(e[:,j]),j])<1e-12) else 'exact_signflip' if n[j]<=16 else 'MC65536_plus_observed'))
     for i,donor in enumerate(ds):
      for j,term in enumerate(terms):
       if np.isfinite(e[i,j]):effects.append(dict(site=node.site,source=src,node_id=node.node_id,parent=node.parent,level=node.level,background=scope,arm=arm,donor=donor,endpoint=term,effect=e[i,j],target_cells=tn[i],rest_cells=rn[i]))
  event(prefix,'node_complete',node=node.node_id)
 df=pd.DataFrame(rows);df['FDR']=np.nan
 for _,ix in df.groupby(['site','source','level','background','arm','family']).groups.items():
  good=df.loc[ix].index[df.loc[ix].P.notna()];df.loc[good,'FDR']=false_discovery_control(df.loc[good,'P'].to_numpy())
 df['nominal']=df.P.lt(.05);df['FDR_significant']=df.FDR.lt(.05);save(df,prefix+'_statistics.tsv');save(pd.DataFrame(effects),prefix+'_patient_effects.tsv.gz');return df
