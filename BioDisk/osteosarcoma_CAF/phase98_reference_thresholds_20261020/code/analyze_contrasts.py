import pathlib,json,hashlib
import numpy as np,pandas as pd
from scipy.stats import t as tdist
from statsmodels.stats.multitest import multipletests
p=pathlib.Path(__file__).resolve().parents[1];root=p.parent
genes=json.loads((p/'resources/PROTOCOL_FROZEN.json').read_text())['genes']
registry=json.loads((root/'phase96_all14_literature_20261016/project_registry.json').read_text())
ids={r['gene']:r['feature_ids'] for r in registry}
lung=pd.read_csv(p/'tables/lung_reference_pseudobulk.tsv',sep='\t',encoding='utf-8')
bone=pd.read_csv(p/'tables/bone_reference_with_coverage_sensitivity.tsv',sep='\t',encoding='utf-8')
specs=[('min'+str(n),'min'+str(n),'min'+str(n)) for n in [1,20,50,100]]
rows=[];tumor_rows=[];meanchecks=[]
def contrast(groups,signs):
 ns=np.array([len(x) for x in groups])
 if (ns<2).any() or any(not np.isfinite(x).all() for x in groups):return dict(status='not_evaluable',effect=np.nan,P=np.nan)
 means=np.array([np.mean(x) for x in groups]);v=np.array([np.var(x,ddof=1)/len(x) for x in groups]);effect=float(means@signs);se=float(np.sqrt(v.sum()))
 if se==0:return dict(status='zero_variance_not_testable',effect=effect,P=np.nan)
 df=float(v.sum()**2/np.sum(v*v/(ns-1)));z=effect/se;pv=float(2*tdist.sf(abs(z),df));hw=float(tdist.ppf(.975,df)*se)
 return dict(status='exploratory_external_reference',effect=effect,SE=se,df=df,P=pv,CI_low=effect-hw,CI_high=effect+hw)
for source in ['SCPCP000017','SCPCP000023']:
 f=root/'phase81_allCAF_bone_lung_20261012/objects'/f'all_cells__expanded__{source}__All_CAF__independent.npz'
 x=np.load(f,allow_pickle=False)
 for gene in genes:
  ix=np.flatnonzero(x['features']==ids[gene]);assert len(ix)==1,(gene,ids[gene])
  tb=x['bone_expression'][:,ix[0]];tl=x['lung_expression'][:,ix[0]]
  for tissue,vals,donors in [('Bone',tb,x['bone_patients']),('Lung',tl,x['lung_patients'])]:
   for donor,val in zip(donors,vals):tumor_rows.append(dict(source=source,donor=str(donor),site=tissue,gene=gene,expression=val,feature=ids[gene]))
  original=next(r for r in registry if r['gene']==gene)
  uncorrected=contrast([tl,tb],np.array([1,-1]))
  meanchecks.append(abs((tl.mean()-tb.mean())-original[source+'__mean_lung_minus_bone']))
  for spec,bs,ls in specs:
   b=bone[(bone.gene==gene)&(bone.specification==bs)]
   l=lung[(lung.gene==gene)&(lung.specification==ls)]
   nb=b.expression.dropna().to_numpy();nl=l.expression.dropna().to_numpy()
   c=contrast([tl,tb,nl,nb],np.array([1,-1,-1,1]))
   refdelta=nl.mean()-nb.mean() if len(nl) and len(nb) else np.nan
   if len(nl)==0 or len(nb)==0:c['status']='reference_gene_unmeasured'
   rows.append(dict(gene=gene,feature=ids[gene],source=source,specification=spec,bone_reference=bs,lung_reference=ls,tumor_bone_N=len(tb),tumor_lung_N=len(tl),reference_bone_N=len(nb),reference_lung_N=len(nl),original_effect=tl.mean()-tb.mean(),original_P=original[source+'__P'],uncorrected_Welch_P=uncorrected['P'],original_FDR_genomewide=original[source+'__FDR'],normal_background_effect=refdelta,**c))
df=pd.DataFrame(rows)
df['FDR_14']=np.nan;df['family']=''
for (src,spec),g in df.groupby(['source','specification']):
 ix=g.index[g.P.notna()]
 if len(ix):df.loc[ix,'FDR_14']=multipletests(df.loc[ix,'P'],method='fdr_bh')[1]
 df.loc[g.index,'family']=src+'|'+spec+'|prespecified14_external_sensitivity'
 df.loc[g.index,'family_valid_tests']=len(ix)
df['nominal']=df.P.lt(.05);df['FDR_significant']=df.FDR_14.lt(.05)
df['same_direction']=np.where(df.effect.notna(),df.effect*df.original_effect>0,np.nan)
df.to_csv(p/'tables/normal_background_all_results.tsv',sep='\t',encoding='utf-8',index=False,na_rep='NA')
df[df.specification.eq('min20')].to_csv(p/'tables/normal_background_primary14.tsv',sep='\t',encoding='utf-8',index=False,na_rep='NA')
pd.DataFrame(tumor_rows).to_csv(p/'tables/original_tumor_patient_expression.tsv',sep='\t',encoding='utf-8',index=False)
qa=dict(original_effect_max_error=float(max(meanchecks)),original_patient_N={s:{'bone':int(g.tumor_bone_N.iloc[0]),'lung':int(g.tumor_lung_N.iloc[0])} for s,g in df.groupby('source')},source_separate=True,reference_variance_included=True,missing_reference_genes_not_filled_with_zero=True,inference='unadjusted_external_fourgroup_contrast; not study-adjusted causal effect',primary_rows=int(df.specification.eq('min20').sum()))
assert max(meanchecks)<1e-10
(p/'resources/CONTRAST_QA.json').write_text(json.dumps(qa,indent=2),encoding='utf-8')
print(df[df.specification.eq('min20')][['gene','source','reference_bone_N','reference_lung_N','original_effect','normal_background_effect','effect','P','FDR_14','status']].to_string(index=False))
