from pathlib import Path
import os,sys,json,datetime
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
import numpy as np,pandas as pd
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase63_matrix_resistance_20261005';sys.path.insert(0,str(P/'phase58_F1F5_round1_20261004/code'))
from common import tests,fdr,mean_valid
def read(p):return pd.read_csv(p,sep='\t',low_memory=False)
panels={'fibrillar_I_III_V':['COL1A1','COL1A2','COL3A1','COL5A1','COL5A2','COL5A3'],'collagen_IV':['COL4A1','COL4A2','COL4A3','COL4A4','COL4A5','COL4A6'],'cartilage_matrix':['COL2A1','ACAN','COL9A1','COL9A2','COL9A3','COL11A1','COL11A2']}
protocol={'at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'role':'posthoc targeted decomposition after primary crosslink pathway showed COL4A1/COL4A2 dominance','reason':'Distinguish interstitial fibrillar collagen I/III/V from typeIV/basement membrane and cartilage-associated composition; no claim these analyst panels are validated classifiers.','panels':panels,'statistics':'Same patient-site effect and coverage rule; separate posthoc family across all33 states and3 panels within source/matrix/scope/minimum. Primary P/FDR unchanged.','selection':'Every subgroup, not only promising states; complete positive/negative results retained.'}
(O/'resources/POSTHOC_COLLAGEN_STRUCTURE_PROTOCOL.json').write_text(json.dumps(protocol,ensure_ascii=False,indent=2),encoding='utf-8')
ob=read(O/'tables/cell_mechanism_scores_with_confidence.tsv.gz').set_index('cell_id');roster=read(O/'tables/paired_patient_site_roster.tsv.gz');z=np.load(O/'objects/mechanism_counts.npz');counts={k:z[k] for k in ['raw','corrected']};den={k:z[k+'_library'] for k in counts};assay=z['assayed'];gi={g:i for i,g in enumerate(z['genes'])};assert set(sum(panels.values(),[]))<=set(gi);mem={k:np.array([gi[g] for g in gs]) for k,gs in panels.items()}
stats=[];effects=[];siteeffects=[]
for (ds,st,mode,scope,minimum),rg in roster.groupby(['dataset','subgroup_id','matrix','scope','minimum_cells']):
 arr=[];donors=[]
 for donor,dg in rg.groupby('donor_record'):
  sitevals=[]
  for row in dg.itertuples():
   base=(ob.dataset.eq(ds)&ob.donor_record.eq(donor)&ob.site.eq(row.site)).to_numpy()
   if scope=='Fibro':base &= ob.lineage_stratum.eq('Fibro').to_numpy()
   il=np.flatnonzero(base&ob.subgroup_id.eq(st).to_numpy());ir=np.flatnonzero(base&ob.subgroup_id.ne(st).to_numpy());assert len(il)==row.target_cells and len(ir)==row.reference_cells
   x=counts[mode];lib=den[mode];present=assay[np.r_[il,ir]].all(0);dv=np.log2(x[il].sum(0)/lib[il].sum()*1e4+.1)-np.log2(x[ir].sum(0)/lib[ir].sum()*1e4+.1);pv=[]
   for name,jj in mem.items():
    have=jj[present[jj]];ok=len(have)>=min(3,len(jj)) and len(have)/len(jj)>=.8;value=float(dv[have].mean()) if ok else np.nan;pv.append(value);siteeffects.append({'dataset':ds,'subgroup_id':st,'matrix':mode,'scope':scope,'minimum_cells':minimum,'donor_record':donor,'site':row.site,'program':name,'effect':value})
   sitevals.append(pv)
  arr.append(mean_valid(np.array(sitevals),axis=0));donors.append(donor)
 d=np.array(arr);p,n,method=tests(d)
 for j,name in enumerate(panels):
  v=d[:,j];v=v[np.isfinite(v)];stats.append({'dataset':ds,'subgroup_id':st,'matrix':mode,'scope':scope,'minimum_cells':minimum,'program':name,'N':int(n[j]),'mean_effect':float(v.mean()) if len(v) else np.nan,'positive_patients':int((v>0).sum()),'negative_patients':int((v<0).sum()),'P':p[j],'test':method[j],'family':f'{ds}|{mode}|{scope}|min{minimum}|posthoc_structure_panels'})
  for donor,value in zip(donors,d[:,j]):effects.append({'dataset':ds,'subgroup_id':st,'matrix':mode,'scope':scope,'minimum_cells':minimum,'program':name,'donor_record':donor,'effect':value})
s=fdr(pd.DataFrame(stats));s.to_csv(O/'tables/posthoc_collagen_structure_statistics.tsv',sep='\t',index=False);pd.DataFrame(effects).to_csv(O/'tables/posthoc_collagen_structure_patient_effects.tsv',sep='\t',index=False);pd.DataFrame(siteeffects).to_csv(O/'tables/posthoc_collagen_structure_site_effects.tsv.gz',sep='\t',index=False)
print(s[(s.dataset=='SCPCP000017')&(s.matrix=='corrected')&(s.scope=='all')&(s.minimum_cells==1)&s.subgroup_id.isin(['F2_C1','F2_C3','F3_C1','F4_C8'])][['subgroup_id','program','N','mean_effect','positive_patients','P','FDR']].to_string(index=False))
(O/'STRUCTURE_COMPLETE').write_text('SUCCESS\n')
