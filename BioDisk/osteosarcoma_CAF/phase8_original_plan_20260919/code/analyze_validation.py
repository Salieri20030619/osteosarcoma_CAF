from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='2'
import json,itertools,datetime
import numpy as np,pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score,confusion_matrix
from statsmodels.stats.multitest import multipletests
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase8_original_plan_20260919'
assert (O/'EXTRACTION_COMPLETE').exists()
panels=json.loads((O/'resources/programs_frozen.json').read_text())
cfg=json.loads((O/'resources/surface_candidates_frozen.json').read_text())
inv=pd.read_csv(O/'tables/analysis_sample_inventory.tsv',sep='\t').set_index('sample')
packs={}
for f in sorted((O/'objects').glob('*_strict_gene_counts.npz')):
 sm=f.name.replace('_strict_gene_counts.npz','');z=np.load(f,allow_pickle=False)
 packs[sm]={k:z[k] for k in z.files}
 # Missing coverage is explicit; arrays contain placeholders only for storage.
assert sum(len(v['raw']) for s,v in packs.items() if s.startswith('GSE152048'))==2047
genes=list(next(iter(packs.values()))['genes'])
customer=sorted(set(panels['Osteogenic_full_customer']+panels['Calcium_regulation_full_customer']))
paired=[];effects=[]
for state in ['0','2']:
 per=[]
 for sm,z in packs.items():
  if not sm.startswith('GSE152048'):continue
  m=z['discovery_states']==state;n=~m
  if m.sum()<20 or n.sum()<20:continue
  c1=z['raw'][m].sum(0)/z['totals'][m].sum()*1e6;c0=z['raw'][n].sum(0)/z['totals'][n].sum()*1e6
  diff=np.log2(c1+1)-np.log2(c0+1)
  for g in customer:
   j=genes.index(g)
   if not z['present'][j]:continue
   row={'state':'H'+state,'sample':sm,'patient_id':inv.at[sm,'patient_id'],'gene':g,'target_cells':int(m.sum()),'other_cells':int(n.sum()),'log2CPM_difference':float(diff[j])}
   paired.append(row);per.append(row)
 per=pd.DataFrame(per)
 for g,df in per.groupby('gene'):
  x=df.log2CPM_difference.to_numpy();n=len(x);p=np.nan
  if n>=3:
   sg=np.array(list(itertools.product([-1.,1.],repeat=n)))
   p=float((np.abs((sg*x).mean(1))>=abs(x.mean())-1e-12).mean())
  effects.append({'state':'H'+state,'gene':g,'patients':n,'median_log2CPM_difference':float(np.median(x)),'positive_patients':int((x>0).sum()),'negative_patients':int((x<0).sum()),'P_value':p})
ef=pd.DataFrame(effects);ef['FDR_BH']=np.nan
for st,ix in ef.groupby('state').groups.items():
 ix=[i for i in ix if pd.notna(ef.at[i,'P_value'])]
 ef.loc[ix,'FDR_BH']=multipletests(ef.loc[ix,'P_value'],method='fdr_bh')[1]
ef.to_csv(O/'tables/full_customer_gene_effects.tsv',sep='\t',index=False)
pd.DataFrame(paired).to_csv(O/'tables/full_customer_patient_effects.tsv',sep='\t',index=False)
# Freeze transfer model before fitting; targets for independent gene checks are excluded from predictors.
heldout_genes=set(customer+cfg['surface_candidates']+cfg['controls'])
features=[g for g in panels['Frozen_H0_top30']+panels['Frozen_H2_top30'] if g not in heldout_genes and all(z['present'][genes.index(g)] for z in packs.values())]
features=list(dict.fromkeys(features));assert len(features)>=15
protocol={'features':features,'excluded_validation_genes':sorted(heldout_genes),'model':'StandardScaler plus L2 multinomial logistic regression C=1; donor-balanced sample weights and balanced class weights','validation':'leave-one-discovery-donor-out','minimum_mean_donor_balanced_accuracy_for_primary_transfer':.6,'confidence_minimum':.7,'interpretation':'candidate-state projection, not new clusters or proof of stable lineage'}
(O/'resources/transfer_protocol.json').write_text(json.dumps(protocol,indent=2))
jj=[genes.index(g) for g in features];X=[];ys=[];donors=[]
for sm,z in packs.items():
 if sm.startswith('GSE152048'):
  X.append(np.log1p(z['raw'][:,jj]/z['totals'][:,None]*1e4));ys.extend(np.where(np.isin(z['discovery_states'],['0','2']),z['discovery_states'],'other'));donors.extend([sm]*len(z['raw']))
X=np.vstack(X);y=np.array(ys);don=np.array(donors)
def fit(x,y,d):
 ss=StandardScaler().fit(x);weights=np.array([1/(d==p).sum() for p in d]);weights*=len(weights)/weights.sum()
 model=LogisticRegression(C=1,class_weight='balanced',max_iter=2000,solver='lbfgs')
 model.fit(ss.transform(x),y,sample_weight=weights);return ss,model
cv=[];pred=np.empty(len(y),dtype=object)
for sm in sorted(set(don)):
 train=don!=sm;test=~train;ss,model=fit(X[train],y[train],don[train])
 pr=model.predict(ss.transform(X[test]));pred[test]=pr
 cv.append({'heldout_donor':sm,'cells':int(test.sum()),'balanced_accuracy':float(balanced_accuracy_score(y[test],pr))})
cv=pd.DataFrame(cv);passed=bool(cv.balanced_accuracy.mean()>=.6)
cv.to_csv(O/'tables/transfer_donor_CV.tsv',sep='\t',index=False)
pd.DataFrame(confusion_matrix(y,pred,labels=['0','2','other']),index=['H0','H2','other'],columns=['H0','H2','other']).to_csv(O/'tables/transfer_CV_confusion.tsv',sep='\t')
ss,model=fit(X,y,don);summary=[];ve=[]
for sm,z in packs.items():
 if sm.startswith('GSE152048'):continue
 assigned={}
 for version,key,tk in [('raw','raw','totals'),('corrected','corrected','corrected_totals')]:
  xm=np.log1p(z[key][:,jj]/np.maximum(z[tk][:,None],1)*1e4);pp=model.predict_proba(ss.transform(xm))
  lab=model.classes_[pp.argmax(1)];confidence=pp.max(1)
  use=(confidence>=.7)&passed
  decision=np.where(use,lab,'unassigned');assigned[version]=decision
  table=pd.DataFrame({'cell_id':z['cell_ids'],'sample':sm,'version':version,'diagnostic_prediction':lab,'confidence':confidence,'accepted_projection':decision,'model_CV_pass':passed})
  table.to_csv(O/'validation'/(sm+'_'+version+'_projection.tsv.gz'),sep='\t',index=False)
  for label in ['0','2','other','unassigned']:
   summary.append({'sample':sm,'dataset':inv.at[sm,'dataset'],'site':inv.at[sm,'site'],'version':version,'projected_state':label,'cells':int((decision==label).sum()),'total':len(lab)})
 # Gene-level validation uses markers not used as classifier features; per donor and dataset.
 for state in ['0','2']:
  m=assigned['raw']==state;n=(assigned['raw']!='unassigned')&~m
  if m.sum()<20 or n.sum()<20:continue
  target=z['raw'][m].sum(0)/z['totals'][m].sum()*1e6;other=z['raw'][n].sum(0)/z['totals'][n].sum()*1e6
  for g in sorted(heldout_genes):
   j=genes.index(g)
   if not z['present'][j]:continue
   ve.append({'sample':sm,'dataset':inv.at[sm,'dataset'],'site':inv.at[sm,'site'],'state':'H'+state,'gene':g,'target_cells':int(m.sum()),'other_cells':int(n.sum()),'log2CPM_difference':float(np.log2(target[j]+1)-np.log2(other[j]+1)),'target_detection_pct':float((z['raw'][m,j]>0).mean()*100),'other_detection_pct':float((z['raw'][n,j]>0).mean()*100),'gene_excluded_from_transfer_features':True})
pd.DataFrame(summary).to_csv(O/'tables/external_projection_counts.tsv',sep='\t',index=False)
pd.DataFrame(ve,columns=['sample','dataset','site','state','gene','target_cells','other_cells','log2CPM_difference','target_detection_pct','other_detection_pct','gene_excluded_from_transfer_features']).to_csv(O/'tables/external_heldout_marker_effects.tsv',sep='\t',index=False)
# Curated protein localization independent of RNA marker performance.
u=json.loads((O/'resources/uniprot_surface_evidence.json').read_text());proteins=[]
for r in u['results']:
 g=r.get('genes',[{}])[0].get('geneName',{}).get('value','')
 loc=[]
 for c in r.get('comments',[]):
  if c.get('commentType')=='SUBCELLULAR LOCATION':
   loc += [s.get('location',{}).get('value','') for s in c.get('subcellularLocations',[])]
 proteins.append({'gene':g,'accession':r['primaryAccession'],'locations':'; '.join(loc),'cell_surface_annotation':any(('Cell membrane' in x or 'Cell surface' in x) for x in loc),'source_url':'https://www.uniprot.org/uniprotkb/'+r['primaryAccession']+'/entry'})
pt=pd.DataFrame(proteins);pt.to_csv(O/'tables/protein_localization.tsv',sep='\t',index=False)
screen=pd.read_csv(O/'tables/RNA_marker_combination_screen.tsv.gz',sep='\t')
rank=[]
for (target,gs),g in screen[(screen.source=='prior')&(screen.detection_count_threshold==1)&(screen.target_cells>=20)].groupby(['target','positive_genes']):
 membrane=all(bool(pt.loc[pt.gene==x,'cell_surface_annotation'].any()) for x in gs.split('+'))
 recall=g.target_coverage.median();purity=g.known_CAF_fraction_gate.fillna(0).median();mal=g.CNV_high_fraction_gate.fillna(0).median()
 rank.append({'target':target,'positive_genes':gs,'donors':g.patient_id.nunique(),'median_RNA_target_coverage':recall,'median_RNA_known_CAF_fraction_gate':purity,'median_RNA_CNV_high_fraction_gate':mal,'worst_RNA_CNV_high_fraction_gate':g.CNV_high_fraction_gate.max(),'surface_localization_supported':membrane,'ranking_score':float(recall*purity*(1-mal)) if membrane else 0,'interpretation':'RNA co-detection prioritization; NOT FACS purity or protein sensitivity'})
rank=pd.DataFrame(rank).sort_values(['target','ranking_score'],ascending=[True,False]);rank.to_csv(O/'tables/surface_combination_ranking.tsv',sep='\t',index=False)
state={'status':'DISCOVERY_COMPLETION_AND_EXTERNAL_VALIDATION_COMPUTED','transfer_features':len(features),'mean_donor_CV_balanced_accuracy':float(cv.balanced_accuracy.mean()),'transfer_CV_pass':passed,'external_samples':len(packs)-6,'full_customer_gene_effect_rows':len(ef),'external_heldout_marker_effect_rows':len(ve),'protein_records':len(pt),'original_labels_changed':False}
(O/'analysis_summary.json').write_text(json.dumps(state,indent=2));(O/'VALIDATION_COMPLETE').write_text('COMPUTATIONAL_VALIDATION_COMPLETE\n')
print(json.dumps(state),flush=True);print(rank[rank.target=='CAF_H0'].head(8).to_string(index=False))
