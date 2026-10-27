import pathlib,json,collections,hashlib
import numpy as np,pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
R=pathlib.Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');P=R/'phase104_full_gene_four_models_20261020'
base=R/'phase81_allCAF_bone_lung_20261012/objects'
olds={s:np.load(base/('all_cells__expanded__'+s+'__All_CAF__independent.npz')) for s in ['SCPCP000017','SCPCP000023']}
nci=np.load(P/'objects/NCI_allgene_author_CAF.npz');lung=np.load(P/'objects/normal_lung_allgene.npz')
bone=pd.read_csv(P/'tables/normal_bone_allgene_raw.tsv',sep='\t');bm=pd.read_csv(P/'tables/normal_bone_donor_totals.tsv',sep='\t')
braw=bone[bm.donor.tolist()].to_numpy().T;bexp=np.log2(braw/bm.library_UMIs.to_numpy()[:,None]*1e4+.1)
reg=json.loads((R/'phase96_all14_literature_20261016/project_registry.json').read_text());target={g['feature_ids']:g['gene'] for g in reg}
aliases=json.loads((R/'phase97_normal_background_20261020/resources/GENE_ALIAS_MAP.json').read_text())
f=olds['SCPCP000017']['features'];symbol=olds['SCPCP000017']['symbols'];idset=set(f);reverse=collections.defaultdict(set)
for a in olds.values():
 for fid,sy in zip(a['features'],a['symbols']):
  if sy and sy!='nan':reverse[str(sy)].add(str(fid))
can=pd.read_csv(R/'phase58_F1F5_round1_20261004/resources/source_canonical_feature_map.tsv',sep='\t')
for r in can[can.dataset.isin(olds)].itertuples():
 if isinstance(r.symbol,str) and r.feature in idset:reverse[r.symbol].add(r.feature)
for g,aa in aliases['aliases'].items():
 for name in aa:reverse[name].add(aliases['stable_ids'][g])
preferred=dict(zip(f,symbol))
for fid,g in target.items():preferred[fid]=g
audits=[]
def map_rows(labels,dataset):
 out={}
 for i,label in enumerate(labels):
  tokens=str(label).split('|')
  stable={x for x in tokens if x in idset}
  hits=stable if stable else set().union(*(reverse.get(x,set()) for x in tokens))
  fid=next(iter(hits)) if len(hits)==1 else None
  audits.append(dict(dataset=dataset,row=i,original_label=str(label),feature=fid or '',status='unique_match' if fid else 'ambiguous' if hits else 'unmapped'))
  if fid:out.setdefault(fid,[]).append(i)
 dup={k for k,v in out.items() if len(v)>1}
 for row in audits:
  if row['dataset']==dataset and row['feature'] in dup:row['status']='multiple_rows_for_feature_excluded'
 return {k:v[0] for k,v in out.items() if len(v)==1}
nm=map_rows(nci['features'],'NCI')
lm=map_rows(lung['features'],'normal_lung')
bmap=map_rows(bone.feature,'normal_bone')
oldmaps={s:{v:i for i,v in enumerate(a['features'])} for s,a in olds.items()}
common=[str(k) for k in f if k in nm and all(k in m for m in oldmaps.values())]
ids=[];sources=[];sites=[];mat=[]
for src,a in olds.items():
 ix=[oldmaps[src][k] for k in common]
 for site in ['Bone','Lung']:
  mat.append(a[site.lower()+'_expression'][:,ix]);n=len(a[site.lower()+'_patients'])
  ids.extend(a[site.lower()+'_patients'].tolist());sources.extend([src]*n);sites.extend([site]*n)
mat.append(nci['expression'][:,[nm[k] for k in common]])
ids.extend(nci['donors'].tolist());sources.extend(['NCI']*17);sites.extend(['Lung']*17)
Y=np.vstack(mat);ids=np.array(ids);sources=np.array(sources);sites=np.array(sites)
finite=np.isfinite(Y).all(0);expressed=(Y>np.log2(.1)+1e-12).any(0)
pd.DataFrame({'feature':common,'complete_patient_coverage':finite,'any_tumor_expression':expressed}).to_csv(P/'tables/tumor_gene_eligibility.tsv',sep='\t',index=False,encoding='utf-8')
keep=finite&expressed;common=np.array(common)[keep];Y=Y[:,keep]
assert set(target)<=set(common) and len(ids)==39 and sum(sites=='Bone')==16
np.savez_compressed(P/'objects/tumor_common_patient_expression.npz',features=common,expression=Y,donors=ids,sources=sources,sites=sites)
pd.DataFrame(audits).to_csv(P/'tables/gene_mapping_audit.tsv',sep='\t',encoding='utf-8',index=False)
protocol=dict(scope='All unambiguously mapped, common, completely observed genes with at least one nonzero tumor pseudobulk; no candidate-only filtering. Four separate test families.',
mapping='Exact stable IDs first, otherwise unique frozen original/canonical symbols or documented target aliases; ambiguous and multiple-row matches excluded; absent != zero.',
models=['pooled_tumor','source_adjusted_tumor','pooled_normal_min20','source_adjusted_normal_min20'],
normal_reference='5 bone Fibro-MSC donors (1439 cells), 3 lung fibroblast donors (2136 cells), >=20 cells each',
patient_scope='Original two main cohorts independent roster +17 NCI patients: 16 Bone,23 Lung,6633 cells; source-specific paired exclusions retained.',
normalization='Raw full-library pseudobulk log2(CP10k+0.1); library denominators not restricted to intersected gene set.',
pooled='Welch two-group unpaired inference; no FC gate',
adjusted='Intercept + Lung + source dummies; HC3 variance; informative-source df=19; NCI source only Lung contributes no within-source site contrast',
normal='Subtract normal Lung minus normal Bone; independent variance sum and Satterthwaite df; require a valid underlying tumor test',
FDR='BH across all valid genes separately in each of four model families; no correction shared across models; original14 highlighted only',
candidate_origin='Previously selected by source-specific nominal P<0.05 with concordant direction, original genomewide FDR not significant',
not_new_strict_cell_analysis=True)
(P/'resources/PROTOCOL_FROZEN.json').write_text(json.dumps(protocol,indent=2),encoding='utf-8')
b=Y[sites=='Bone'];l=Y[sites=='Lung']
bv=b.var(0,ddof=1)/len(b);lv=l.var(0,ddof=1)/len(l)
e_pool=l.mean(0)-b.mean(0);v_pool=bv+lv;d_pool=np.divide(v_pool**2,bv**2/(len(b)-1)+lv**2/(len(l)-1),out=np.full_like(v_pool,np.nan),where=v_pool>1e-24)
X=np.column_stack([np.ones(len(ids)),(sites=='Lung').astype(float),pd.get_dummies(sources,drop_first=True,dtype=float).to_numpy()])
inv=np.linalg.inv(X.T@X);beta=inv@X.T@Y;res=Y-X@beta;h=np.einsum('ij,jk,ik->i',X,inv,X)
w=(inv@X.T)[1];v_adj=((w/(1-h))[:,None]**2*res**2).sum(0);e_adj=beta[1]
assert np.max(np.abs(w[sources=='NCI']))<1e-12
rows=[];results={}
def evaluate(e,v,df,valid):
 se=np.sqrt(v);pv=np.full(len(e),np.nan);lo=np.full(len(e),np.nan);hi=np.full(len(e),np.nan)
 valid=valid&np.isfinite(e)&np.isfinite(v)&np.isfinite(df)&(v>1e-24)&(df>0)
 pv[valid]=2*stats.t.sf(np.abs(e[valid]/se[valid]),df[valid]);ci=stats.t.ppf(.975,df[valid])*se[valid]
 lo[valid]=e[valid]-ci;hi[valid]=e[valid]+ci
 q=np.full(len(e),np.nan);q[valid]=multipletests(pv[valid],method='fdr_bh')[1]
 return se,pv,q,lo,hi,valid
normal_mask=np.array([k in bmap and k in lm for k in common])
NE=np.full(len(common),np.nan);NVL=NE.copy();NVB=NE.copy()
for j,k in enumerate(common):
 if normal_mask[j]:
  bb=bexp[:,bmap[k]];ll=lung['expression'][:,lm[k]]
  NE[j]=ll.mean()-bb.mean();NVL[j]=ll.var(ddof=1)/3;NVB[j]=bb.var(ddof=1)/5
for model,e,v,df in [('pooled',e_pool,v_pool,d_pool),('source_adjusted',e_adj,v_adj,np.full(len(common),19.))]:
 valid0=np.isfinite(v)&(v>1e-24)
 for normal in [False,True]:
  name=model+('_normal_min20' if normal else '_tumor')
  effect=e-NE if normal else e
  var=v+NVL+NVB if normal else v
  if normal:
   den=v*v/df+NVL*NVL/2+NVB*NVB/4
   deg=np.divide(var*var,den,out=np.full_like(var,np.nan),where=den>1e-48)
  else:deg=df
  se,pv,q,lo,hi,valid=evaluate(effect,var,deg,valid0&(normal_mask if normal else True))
  data=pd.DataFrame(dict(model=name,feature=common,gene=[preferred.get(k,k) for k in common],effect=effect,SE=se,df=deg,P=pv,FDR=q,CI_low=lo,CI_high=hi,valid_test=valid,normal_feature_available=normal_mask,is_original14=[k in target for k in common]))
  data['significance']=np.select([data.FDR.lt(.05),data.P.lt(.05)],['FDR<0.05','P<0.05_only'],default='not_significant')
  data.loc[~valid,'significance']='not_evaluable';data['family_valid_tests']=int(valid.sum())
  data.to_csv(P/'tables'/(name+'_allgenes.tsv.gz'),sep='\t',index=False,encoding='utf-8',na_rep='NA')
  results[name]=data
  rows.append(dict(model=name,common_tumor_eligible=len(common),normal_covered=int(normal_mask.sum()),valid_tests=int(valid.sum()),nominal=int((pv<.05).sum()),FDR_significant=int((q<.05).sum()),P_only=int(((pv<.05)&~(q<.05)).sum()),FDR_lung_higher=int(((q<.05)&(effect>0)).sum()),FDR_bone_higher=int(((q<.05)&(effect<0)).sum())))
out=pd.concat(results.values(),ignore_index=True);c=out[out.is_original14].copy();c['same_original_direction']=np.where(c.gene.eq('EBF1'),c.effect<0,c.effect>0)
c.to_csv(P/'tables/original14_full_family_results.tsv',sep='\t',index=False,encoding='utf-8',na_rep='NA')
pd.DataFrame(rows).to_csv(P/'tables/four_analysis_summary.tsv',sep='\t',encoding='utf-8',index=False)
# Previous14 estimates and raw P must reproduce; only FDR family changes.
old=pd.read_csv(R/'phase103_pooled_NCI_sensitivities_20261020/tables/all_14gene_results.tsv',sep='\t')
old=old[old.cellset.eq('expanded')&old.cohort.eq('with_NCI')&~old.model.eq('source_adjusted_OLS_HC3_global_df')]
errs=[]
for r in c.itertuples():
 mode='pooled_Welch' if r.model.startswith('pooled') else 'source_adjusted_OLS_HC3'
 ct='normal_min20' if 'normal' in r.model else 'tumor_only'
 prev=old[old.gene.eq(r.gene)&old.model.eq(mode)&old.contrast.eq(ct)].iloc[0]
 for key in ['effect','P','SE']:
  v=getattr(r,key);v0=prev[key]
  if np.isnan(v0):assert np.isnan(v)
  else:errs.append(abs(v-v0));assert abs(v-v0)<1e-10,(r.gene,key,v,v0)
qa=dict(status='PASS',common_tumor_eligible=len(common),normal_common=int(normal_mask.sum()),all39_complete_for_selected_genes=True,previous14_reproduced_checks=len(errs),previous14_max_error=float(max(errs)),NCI_total_cells=int(nci['cells'].sum()),normal_bone_cells=int(bm.cells.sum()),normal_lung_cells=int(lung['cells'].sum()),library_denominators_preserved=True,strict_subset_not_used=True)
(P/'resources/ANALYSIS_QA.json').write_text(json.dumps(qa,indent=2),encoding='utf-8')
print(pd.DataFrame(rows).to_string(index=False));print(c[['model','gene','effect','P','FDR','same_original_direction']].to_string(index=False));print(json.dumps(qa))
