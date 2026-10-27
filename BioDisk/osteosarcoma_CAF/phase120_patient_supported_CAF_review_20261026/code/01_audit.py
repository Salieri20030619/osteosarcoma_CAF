from pathlib import Path
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']: os.environ[k]='2'
import pandas as pd,numpy as np,json,hashlib,datetime
from scipy import sparse
from scipy.stats import false_discovery_control
R=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF'); I=R/'phase119_Bone_Lung_functional_subgroups_20261024'; O=R/'phase120_patient_supported_CAF_review_20261026'
for sub in ['code','resources','tables','logs','delivery']: (O/sub).mkdir(parents=True,exist_ok=True)
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
def dump(n,x): (O/n).write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str)+'\n',encoding='utf-8')
def read(n):return pd.read_csv(I/'tables'/n,sep='\t',low_memory=False)
def save(x,n):x.to_csv(O/'tables'/n,sep='\t',index=False,encoding='utf-8')
def event(s,**kw):
 x=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),stage=s,**kw)
 with (O/'events.jsonl').open('a',encoding='utf-8') as f:f.write(json.dumps(x,default=str)+'\n')
 print(json.dumps(x,default=str),flush=True)
rule=dict(version='v1_frozen_before_new_evaluation',scope='Fixed6633; no reclustering or scientific figures',minimum_evaluable_patients=6,preferred_patients=10,minimum_target_and_control_cells=5,positive_patient_fraction=.75,maximum_patient_cell_fraction=.5,seed_jaccard_min=.8,primary_support='Within one source: >=6 independent donors with >=5 target AND control cells; no pooling source-local records into global N',replication='Additional source >=3 evaluable donors and >=75% same direction; this is replication support, not external independent validation',functional='At least one positive nominal P<.05 collagen production/matrix endpoint, >=75% positive, leave-one-donor-out mean stays positive; a non-collagen processing panel must also have positive mean and >=75% positive direction. FDR reported, not sole gate.',identity='Retain historical fibro-supported cells; NCI additionally restrict to author diploid; remove only explicit immune/endothelial coexpression in sensitivity, not source data. This is a sensitivity context, not validated purity.',identity_sensitivity='Same qualifying primary endpoint must remain positive with >=75% positive in >=6 min5 donor pairs in identity sensitivity; nominal P not required twice.',nested='Children primary contrast with siblings; parents with site remainder; separately retain fibro-context comparison',effect_threshold='No universal numeric score cut-off invented; report effect, median, consistency, LOO and min10 sensitivity',selection_bias='Post-discovery review, not independent confirmation; existing heldout findings retained; no threshold tuning after outputs')
if (O/'PROTOCOL_V1.json').exists():assert json.loads((O/'PROTOCOL_V1.json').read_text())==rule
else:dump('PROTOCOL_V1.json',rule)
files=[I/'tables'/n for n in ['cohort_registry.tsv.gz','hierarchical_cell_labels.tsv.gz','hierarchical_nodes.tsv','main_patient_effects.tsv.gz','main_statistics.tsv','cell_function_scores.tsv.gz','genes.tsv','target_collagen_candidate_review.tsv','node_inherited_CNV_records.tsv','holdout_statistics.tsv']]+[I/'objects/log1p_cp10k.npz',I/'resources/PROGRAMS.json']
before={str(f):sha(f) for f in files};dump('INPUT_HASHES.json',before);event('inputs_hashed')
m=read('cohort_registry.tsv.gz');l=read('hierarchical_cell_labels.tsv.gz');nodes=read('hierarchical_nodes.tsv');s=read('cell_function_scores.tsv.gz');rv=read('target_collagen_candidate_review.tsv')
assert len(m)==6633 and m.cell_id.equals(s.cell_id) and m.cell_id.equals(l.cell_id)
assert (~m.doublet).sum()==6423
print('REGISTRY',m.columns.tolist(),'LABELS',l.columns.tolist(),flush=True)
ncpath=R/'phase102_NCI_author_CAF_20261020/tables/author_CAF_cell_registry.tsv.gz'
nc=pd.read_csv(ncpath,sep='\t').set_index('cell_id');m['nci_cnv']=m.cell_id.map(nc.native_copykat_call)
genes=read('genes.tsv').symbol;X=sparse.load_npz(I/'objects/log1p_cp10k.npz');gi={g:i for i,g in enumerate(genes)}
panels={
 'FNDC1_paper':['FNDC1','COMP','MMP2','COL11A1','COL10A1','P4HA3'],
 'MYH11_paper':['MYH11','RGS5','MCAM','MYH9','TGFBI'],
 'CD36_paper':['CD36','FABP4','APOE','STEAP4','APOLD1'],
 'C7_paper':['C7','PTGDS','TNXB','FGL2','IGF2','CXCL12'],
 'fibro_identity':['DCN','LUM','PDGFRA','COL1A1','COL1A2','COL3A1'],
 'contractile':['ACTA2','TAGLN','MYL9','MYH11'],
 'collagen_processing_nonCOL':['SERPINH1','P4HA1','P4HA2','P4HA3','PLOD1','PLOD2','PLOD3','PCOLCE','PCOLCE2','BMP1','LOX','LOXL1','LOXL2'],
 'immune_identity':['PTPRC','LST1','TYROBP','CD3D','CD3E'],
 'endothelial_identity':['PECAM1','VWF','CDH5','KDR'],
 'osteogenic_identity':['RUNX2','SP7','ALPL','IBSP','BGLAP','PHOSPHO1'],
 'CXCL14_ITGA11_context':['ITGA11','ITGB1','TGFB1'],
}
dump('resources/NEW_PANELS.json',dict(panels=panels,scope='Paper panels exactly from primary text; remaining panels operational sensitivity panels, not literature-derived universal classifiers. CXCL14 receptor axis is contextual only, not evidence of cell-cell mechanism.',coverage={k:dict(present=[g for g in v if g in gi],missing=[g for g in v if g not in gi]) for k,v in panels.items()}))
for name,gs in panels.items():s['audit__'+name]=np.asarray(X[:,[gi[g] for g in gs if g in gi]].mean(1)).ravel()
def positive(g):return X[:,gi[g]].toarray().ravel()>0 if g in gi else np.zeros(len(m),bool)
immune=positive('PTPRC')&(positive('LST1')|positive('TYROBP')|positive('CD3D')|positive('CD3E'))
endo=positive('PECAM1')&(positive('VWF')|positive('CDH5'))
ident=m.fibro_supported.to_numpy()&~immune&~endo&(~m.source.eq('NCI')|m.nci_cnv.eq('diploid')).to_numpy()
m['identity_sensitivity']=ident;m['immune_coexpression']=immune;m['endothelial_coexpression']=endo
save(m,'cell_identity_audit.tsv.gz')
terms=[c for c in s if c.startswith(('program__','rank__','audit__'))]; vals=s[terms].to_numpy(); out=[]; effects=[];support=[]; composition=[]; panelgenes=sorted(set(sum(panels.values(),[])))
signcache={}
def stat(v):
 v=np.asarray(v);v=v[np.isfinite(v)];n=len(v)
 if not n:return (0,np.nan,np.nan,0,np.nan,np.nan)
 if n<2 or np.all(abs(v)<1e-12):p=np.nan
 elif n<=16:
  if n not in signcache:signcache[n]=1-2*((np.arange(2**n)[:,None]>>np.arange(n))&1)
  p=np.mean(abs(signcache[n]@v)>=abs(v.sum())-1e-10)
 else:
  a=np.random.default_rng(12026+n).choice([-1,1],size=(65536,n));p=(np.sum(abs(a@v)>=abs(v.sum())-1e-10)+1)/65537
 return n,v.mean(),np.median(v),int((v>1e-12).sum()),p,float(((v.sum()-v)/(n-1)).min()) if n>1 else np.nan
for node in nodes.itertuples():
 col='parent' if node.level=='parent' else 'leaf';target=l[col].eq(node.node_id).to_numpy();base=m.site.eq(node.site).to_numpy()&~m.doublet.to_numpy()
 for (src,donor),g in m[base&target].groupby(['source','donor']):
  ix=g.index;support.append(dict(node_id=node.node_id,site=node.site,level=node.level,source=src,donor=donor,cells=len(g),identity_cells=int(ident[ix].sum()),immune_cells=int(immune[ix].sum()),endothelial_cells=int(endo[ix].sum()),nci_aneuploid=int(g.nci_cnv.eq('aneuploid').sum()),old_osteo=int(g.old_lineage.eq('Osteo').sum()),old_fibro=int(g.old_lineage.eq('Fibro').sum())))
 for src,g in m[base&target].groupby('source'):
  for gene in panelgenes:
   if gene not in gi:continue
   a=X[g.index,gi[gene]].toarray().ravel();composition.append(dict(node_id=node.node_id,source=src,gene=gene,cells=len(g),mean=float(a.mean()),detection=float((a>0).mean())))
 bg=base.copy()
 if node.level=='child':bg&=l.parent.eq(node.parent).to_numpy()
 for mode,mask,minimum in [('main_min1',bg,1),('main_min5',bg,5),('main_min10',bg,10),('identity_min5',bg&ident,5),('fibro_context_min5',base&ident,5)]:
  for src in m.loc[base,'source'].unique():
   records=[]
   for donor,g in m[mask&m.source.eq(src).to_numpy()].groupby('donor'):
    a=g.index[target[g.index]];b=g.index[~target[g.index]]
    if min(len(a),len(b))<minimum:continue
    delta=vals[a].mean(0)-vals[b].mean(0);records.append(delta)
    for j,t in enumerate(terms):effects.append(dict(node_id=node.node_id,site=node.site,level=node.level,source=src,mode=mode,donor=donor,endpoint=t,effect=delta[j],target_cells=len(a),control_cells=len(b)))
   for j,t in enumerate(terms):
    v=[a[j] for a in records];n,mu,median,pos,p,loo=stat(v)
    out.append(dict(node_id=node.node_id,site=node.site,level=node.level,source=src,mode=mode,endpoint=t,family=t.split('__')[0],N=n,mean_effect=mu,median_effect=median,positive_patients=pos,P=p,LOO_min_mean=loo))
 event('node_complete',node=node.node_id)
df=pd.DataFrame(out);df['FDR']=np.nan
for key,ix in df.groupby(['site','source','level','mode','family']).groups.items():
 good=df.loc[ix].index[df.loc[ix].P.notna()];df.loc[good,'FDR']=false_discovery_control(df.loc[good].P.to_numpy()) if len(good) else np.nan
save(df,'functional_audit.tsv');save(pd.DataFrame(effects),'patient_effects.tsv.gz');sp=pd.DataFrame(support);save(sp,'node_donor_support.tsv');save(pd.DataFrame(composition),'marker_detection.tsv')
production=['program__collagen_biosynthesis','program__collagen_assembly','program__collagen_crosslinking','program__fibrillar_I_III_V','program__mCAF_matrix','program__F2C3_collagen_support']
dec=[]
for x in rv.itertuples():
 z=df[df.node_id.eq(x.node_id)];counts=sp[sp.node_id.eq(x.node_id)];best=[]
 for src,g in z[z['mode'].eq('main_min5')].groupby('source'):
  q=g[g.endpoint.isin(production)];passing=q[q.N.ge(6)&q.mean_effect.gt(0)&q.P.lt(.05)&(q.positive_patients/q.N>=.75)&q.LOO_min_mean.gt(0)]
  p=g[g.endpoint.eq('audit__collagen_processing_nonCOL')].iloc[0];processing=p.N>=6 and p.mean_effect>0 and p.positive_patients/p.N>=.75 and p.LOO_min_mean>0
  identity=[]
  for end in passing.endpoint:
   a=z[z.source.eq(src)&z['mode'].eq('identity_min5')&z.endpoint.eq(end)].iloc[0]
   if a.N>=6 and a.mean_effect>0 and a.positive_patients/a.N>=.75 and a.LOO_min_mean>0:identity.append(end)
  c=counts[counts.source.eq(src)].cells;dom=c.max()/c.sum();stable=x.effective_seed_jaccard>=.8
  passed=bool(len(identity)>0 and processing and dom<=.5 and stable)
  best.append(dict(source=src,N_min5=int(q.N.max()),N_min10=int(z[z.source.eq(src)&z['mode'].eq('main_min10')].N.max()),N_identity5=int(z[z.source.eq(src)&z['mode'].eq('identity_min5')].N.max()),nominal_endpoints=';'.join(passing.endpoint),identity_supported_endpoints=';'.join(identity),processing_support=bool(processing),source_patient_dominance=float(dom),stable=bool(stable),qualified=passed))
 for a in best:dec.append(dict(node_id=x.node_id,site=x.site,level=x.level,fit_cells=x.fit_cells,total_source_patient_records=x.source_patients,**a))
decision=pd.DataFrame(dec);save(decision,'source_candidate_gate.tsv')
old=read('main_statistics.tsv');check=df[df['mode'].isin(['main_min1','main_min5'])&df.family.isin(['program','rank'])].copy();check['arm']=check['mode'].map({'main_min1':'singlet','main_min5':'singlet_min5'});check['background']=np.where(check.level.eq('parent'),'within_site','within_parent');joined=check.merge(old,on=['node_id','source','endpoint','arm','background'],suffixes=('_new','_old'),validate='one_to_one');assert (joined.N_new==joined.N_old).all();err=float((joined.mean_effect_new-joined.mean_effect_old).abs().max());assert err<1e-10
exact=joined[joined.N_new.le(16)];perr=float((exact.P_new-exact.P_old).abs().max());assert perr<1e-9
assert all(sha(p)==h for p,h in before.items())
dump('AUDIT_COMPLETE.json',dict(status='PASS',nodes=len(nodes),parents=int(nodes.level.eq('parent').sum()),cells=len(m),main_cells=int((~m.doublet).sum()),recomputed_comparisons=len(joined),max_effect_error=err,max_exact_P_error=perr,input_hashes_unchanged=True,qualified_source_nodes=decision[decision.qualified].to_dict('records'),remaining=['literature synthesis','candidate review','report and delivery QA']))
print(decision[decision.node_id.isin(['Bone_07','Bone_07_s1','Lung_04','Lung_04_s2'])].to_string(index=False),flush=True)
event('audit_complete',qualified=int(decision.qualified.sum()))
