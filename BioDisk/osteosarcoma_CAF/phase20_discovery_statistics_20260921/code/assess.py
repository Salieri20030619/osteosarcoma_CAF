from pathlib import Path
import json,hashlib,datetime
import numpy as np,pandas as pd
from scipy.stats import rankdata,binomtest,wilcoxon
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
B=Path('C:/codex_test');I=B/'CAF_evidence_20260921';O=B/'CAF_discovery_stats_20260921'
COHORTS=['GSE152048','GSE162454','GSE318840'];PROGS=['myCAF_customer','iCAF_customer','apCAF_customer','Osteogenic_customer','Mineral_customer']
protocol={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'primary_question':'Are frozen discovery-cluster programs elevated within patients, especially supported C7/C9? External outcome significance is not an acceptance condition.',
 'primary_cohorts':COHORTS,'cohort_rationale':'separate reported recruitment: Shanghai Sixth Peoples Hospital11 patients; Guangxi Medical University6 patients; Vinmec Hanoi8 patients. Source-specific IDs used once per comparison. Cohort independence inferred from separate recruitment descriptions, not identifiable-person matching.',
 'other_cohorts':'retained descriptively and per-study tests; excluded from pooled primary N because cross-source reuse/alias mapping not resolved',
 'primary_layer':'supported','primary_test':'two-sided exact paired Wilcoxon signed-rank, enumerated integer-rank distribution including tied ranks; drop differences rounded to zero at12 decimals; n>=2 independent patients',
 'unit':'one within-patient cluster-minus-rest difference in full-source logCP10k panel score; at least10 cells on both sides; panel coverage>=80%; no cells-as-replicates',
 'FDR':'BH across all18 clusters x5 programs=90 in primary supported layer; separate sensitivity families for strict50 and expanded120; unassessable tests use1 in BH pool then display NA. Genes corrected in separate per-layer family.',
 'sensitivity':['exact two-sided sign test','per-cohort tests','leave-one-primary-cohort-out','restrict target and comparator to raw-panel-positive myCAF cells, same10-cell rule and per-layer BH'],
 'effects':'median patient difference, positive/negative counts, 95% percentile CI via5000 within-cohort patient bootstraps',
 'interpretation':'discovery/post-selection exploratory tests: clusters and marker candidates were selected from same expression data; not a selective-inference-corrected cluster-existence test, not proof of a novel stable lineage or metastatic association',
 'sources':{'GSE152048':'https://pmc.ncbi.nlm.nih.gov/articles/PMC7730477/','GSE162454':'https://www.frontiersin.org/journals/oncology/articles/10.3389/fonc.2021.709210/full','GSE318840':'https://pmc.ncbi.nlm.nih.gov/articles/PMC13330631/','signed_rank':'https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.wilcoxon.html'}}
(O/'PROTOCOL.json').write_text(json.dumps(protocol,ensure_ascii=False,indent=2),encoding='utf-8')
def read(n):return pd.read_csv(I/'tables'/n,sep='\t',low_memory=False)
def save(d,n):d.to_csv(O/'tables'/n,sep='\t',index=False)
def exact(v):
 d=np.round(np.asarray(v,dtype=float),12);d=d[d!=0]
 if not len(d):return 1.,0.,0
 ranks=np.rint(rankdata(abs(d))*2).astype(int);total=int(ranks.sum());dist=np.zeros(total+1,dtype=np.int64);dist[0]=1;last=0
 for r in ranks:
  old=dist[:last+1].copy();dist[r:r+last+1]+=old;last+=r
 assert int(dist.sum())==2**len(d)
 plus=int(ranks[d>0].sum());tail=min(plus,total-plus)
 p=min(1.,2*int(dist[:tail+1].sum())/2**len(d))
 assert p>=2/2**len(d)-1e-15
 return p,plus/2,len(d)
# Validate exact engine on distinct, tied, zero and symmetric cases.
for v in [[1,2,3,4,5,6],[1,2,-3,4,-5,6],[-1,-2,-3,-4,-5,-6]]:
 assert abs(exact(v)[0]-wilcoxon(v,method='exact').pvalue)<1e-12
assert exact([1,1,1,1])[0]==.125 and exact([0,0,0])[0]==1 and exact([1,-1])[0]==1
rng=np.random.default_rng(20260921)
def summarize(g,ci=True):
 v=g.delta.dropna().to_numpy(float);n=len(v);nz=np.round(v,12);nz=nz[nz!=0];pos=int((v>1e-12).sum());neg=int((v< -1e-12).sum())
 out={'n_patients':n,'n_datasets':g.dataset.nunique(),'positive_patients':pos,'negative_patients':neg,'median_delta':float(np.median(v)) if n else np.nan,'mean_delta':float(np.mean(v)) if n else np.nan,'P':np.nan,'P_sign':np.nan,'CI_low':np.nan,'CI_high':np.nan,'minimum_attainable_P':2/2**len(nz) if len(nz) else np.nan}
 if n>=2:
  out['P']=exact(v)[0];out['P_sign']=binomtest(pos,pos+neg,p=.5).pvalue if pos+neg else 1.
  if ci:
   draws=[]
   for _,gg in g.groupby('dataset'):
    a=gg.delta.dropna().to_numpy();draws.append(rng.choice(a,(5000,len(a)),replace=True))
   bs=np.median(np.concatenate(draws,axis=1),axis=1);out['CI_low'],out['CI_high']=np.quantile(bs,[.025,.975])
 out['evidence_position']='exploratory_small_sample' if n<10 else 'discovery_postselection_inference'
 return out
def adjust(df,groups):
 df['FDR']=np.nan
 for _,ix in df.groupby(groups).groups.items():
  p=df.loc[ix,'P'];q=multipletests(p.fillna(1),method='fdr_bh')[1];df.loc[ix,'FDR']=np.where(p.notna(),q,np.nan)
 df['result']=np.select([df.P.isna(),(df.P<.05)&(df.FDR<.05),df.P<.05],['not_evaluable','FDR_significant_discovery','nominally_significant_exploratory'],default='not_significant')
 df['enrichment_supported']=df.P.lt(.05)&df.median_delta.gt(0)
 return df
E=read('within_donor_program_effects.tsv');E['cluster']=E.cluster.astype(str)
eligible=E[E.eligible&E.dataset.isin(COHORTS)].copy();assert not eligible.donor_record.str.contains('UNCONFIRMED|UNKNOWN_PATIENT',case=False).any()
save(eligible,'primary_patient_program_differences.tsv')
master=read('cluster_identity_and_support.tsv');master['cluster']=master.cluster.astype(str)
rows=[];within=[];loo=[]
for r in master.to_dict('records'):
 for pr in PROGS:
  g=eligible[(eligible.layer==r['layer'])&(eligible.cluster==r['cluster'])&(eligible.program==pr)]
  assert not g.donor_record.duplicated().any()
  key={'layer':r['layer'],'cluster':r['cluster'],'program':pr};rows.append({**key,**summarize(g)})
  for ds in COHORTS:within.append({**key,'dataset':ds,**summarize(g[g.dataset==ds])})
  for ds in COHORTS:loo.append({**key,'excluded_dataset':ds,**summarize(g[g.dataset!=ds],False)})
R=adjust(pd.DataFrame(rows),['layer']);W=adjust(pd.DataFrame(within),['layer','dataset']);LO=adjust(pd.DataFrame(loo),['layer','excluded_dataset'])
save(R,'discovery_program_tests.tsv');save(W,'within_cohort_program_tests.tsv');save(LO,'leave_one_cohort_out.tsv')
# Same program tests within only the myCAF-positive cell compartment.
myrows=[];mytests=[];alloc=[]
for layer in ['strict','supported','expanded']:
 d=read(f'{layer}_analysis_labels.tsv.gz');d['cluster']=d.cluster.astype(str)
 dd=d[d.dataset.isin(COHORTS)&d.known_patient&d.myCAF_customer_positive].copy()
 for cl in master.loc[master.layer==layer,'cluster']:
  for donor,g in dd.groupby('analysis_donor'):
   ins=g[g.cluster==cl];outs=g[g.cluster!=cl]
   for pr in PROGS:
    x=ins[pr+'_score'].dropna();y=outs[pr+'_score'].dropna()
    if len(x)>=10 and len(y)>=10:myrows.append({'layer':layer,'cluster':cl,'donor_record':donor,'dataset':g.dataset.iloc[0],'program':pr,'inside_cells':len(x),'outside_cells':len(y),'delta':x.mean()-y.mean()})
 alloc.extend(dd[['sample','patient_id','analysis_donor','dataset']].drop_duplicates().to_dict('records'))
ME=pd.DataFrame(myrows)
for r in master.to_dict('records'):
 for pr in PROGS:
  g=ME[(ME.layer==r['layer'])&(ME.cluster==r['cluster'])&(ME.program==pr)]
  mytests.append({'layer':r['layer'],'cluster':r['cluster'],'program':pr,**summarize(g)})
M=adjust(pd.DataFrame(mytests),['layer']);save(M,'within_myCAF_program_tests.tsv');save(ME,'within_myCAF_patient_differences.tsv');save(pd.DataFrame(alloc).drop_duplicates(),'primary_patient_roster.tsv')
# Candidate genes are a separate exploratory hypothesis family, not mixed with programs.
GE=read('gene_candidate_donor_effects.tsv');GE['cluster']=GE.cluster.astype(str);GE=GE[GE.dataset.isin(COHORTS)]
gene_keys=read('gene_candidate_evidence.tsv')[['layer','cluster','gene_id','symbol']].drop_duplicates();gene_keys['cluster']=gene_keys.cluster.astype(str)
gg=[]
for r in gene_keys.to_dict('records'):
 g=GE[(GE.layer==r['layer'])&(GE.cluster==r['cluster'])&(GE.gene_id==r['gene_id'])]
 gg.append({**r,**summarize(g,False)})
G=adjust(pd.DataFrame(gg),['layer']);save(G,'discovery_marker_tests.tsv')
focus=R[(R.layer=='supported')&R.cluster.isin(['7','9'])].copy();save(focus,'C7_C9_primary_results.tsv')
# A stronger donor-count requirement is a descriptive sensitivity, not parameter tuning.
E20=eligible[(eligible.inside_cells>=20)&(eligible.outside_cells>=20)];r20=[]
for r in master.to_dict('records'):
 for pr in PROGS:r20.append({'layer':r['layer'],'cluster':r['cluster'],'program':pr,**summarize(E20[(E20.layer==r['layer'])&(E20.cluster==r['cluster'])&(E20.program==pr)],False)})
save(adjust(pd.DataFrame(r20),['layer']),'minimum20_cell_sensitivity.tsv')
fig,axs=plt.subplots(1,2,figsize=(11,4.8),sharey=True)
colors={'GSE152048':'#2878a2','GSE162454':'#da8030','GSE318840':'#4d9b66'}
for ax,cl in zip(axs,['7','9']):
 for j,pr in enumerate(PROGS):
  g=eligible[(eligible.layer=='supported')&(eligible.cluster==cl)&(eligible.program==pr)]
  for ds,zz in g.groupby('dataset'):
   ax.scatter(zz.delta,np.full(len(zz),j)+rng.uniform(-.13,.13,len(zz)),s=24,alpha=.8,color=colors[ds],label=ds if j==0 else None)
  rr=focus[(focus.cluster==cl)&(focus.program==pr)].iloc[0];ax.plot([rr.median_delta],[j],marker='|',markersize=15,color='black')
 ax.axvline(0,color='gray',ls='--',lw=.8);ax.set_title('Supported C'+cl);ax.set_xlabel('Within-patient mean logCP10k difference');ax.set_yticks(range(5),['myCAF','iCAF','apCAF','Osteogenic','Mineral']);ax.legend(fontsize=8,loc='best')
fig.suptitle('Discovery cluster minus other cells; each point is one patient',fontsize=12);fig.tight_layout();fig.savefig(O/'figures/patient_program_effects.png',dpi=180);plt.close(fig)
qa={'primary_cohorts':COHORTS,'program_tests':len(R),'supported_family_size':90,'gene_tests':len(G),'no_cell_level_pseudoreplication':True,'unknown_patients_not_pooled':True,'exact_test_unit_checks_pass':True,'external_validation_required':False,'source_files_unchanged':True,'scope':'exploratory post-selection discovery evidence, not new-subtype existence theorem'}
assert len(R)==260 and len(focus)==10
(O/'delivery/QA.json').write_text(json.dumps(qa,indent=2));print(focus.to_string(index=False),flush=True);print('WITHIN_MYCAF',flush=True);print(M[(M.layer=='supported')&M.cluster.isin(['7','9'])].to_string(index=False),flush=True)
