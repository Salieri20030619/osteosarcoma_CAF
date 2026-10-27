from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
import numpy as np,pandas as pd,json,itertools,hashlib,datetime
from statsmodels.stats.multitest import multipletests
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase56_patient_programs_20261004'
def read(p):return pd.read_csv(p,sep='\t',low_memory=False)
def dump(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str),encoding='utf-8')
signcache={}
def tests(d):
 d=np.asarray(d,dtype=np.float64);n=d.shape[0];p=np.full(d.shape[1],np.nan);method=np.full(d.shape[1],'not_tested_N_lt2',dtype=object)
 if n<2:return p,method
 nz=np.any(d!=0,axis=0);method[~nz]='not_tested_all_zero_effect';ix=np.flatnonzero(nz)
 if n not in signcache:signcache[n]=np.array(list(itertools.product([-1.,1.],repeat=n))) if n<=15 else np.random.default_rng(20261004).choice([-1.,1.],size=(19999,n))
 signs=signcache[n];obs=np.abs(d.mean(0));label='exact_patient_signflip' if n<=15 else 'patient_signflip_MC19999_plus_observed'
 for begin in range(0,len(ix),128):
  j=ix[begin:begin+128];v=d[:,j];hits=(np.abs(signs@v/n)>=obs[j][None,:]-1e-12).sum(0);p[j]=hits/len(signs) if n<=15 else (hits+1)/(len(signs)+1)
 method[ix]=label;return p,method

def bh(p):
 p=np.asarray(p,dtype=float);q=multipletests(np.where(np.isfinite(p),p,1.),method='fdr_bh')[1];q[~np.isfinite(p)]=np.nan;return q

def main():
 assert (O/'PSEUDOBULK_COMPLETE').is_file() and (O/'pseudobulk.exit_code').read_text().strip()=='0'
 m=read(O/'tables/pseudobulk_sample_metadata.tsv');proto=json.loads((O/'resources/PROTOCOL_FROZEN.json').read_text());gs=json.loads((O/'resources/PROGRAMS_FROZEN.json').read_text());sets=gs['programs'];programs=list(sets);scopes=proto['scopes'];rows=[];patientrows=[];siteledger=[];contexts=[];contributions=[]
 for source,ms in m.groupby('source_id',sort=True):
  objs={row.sample_id:np.load(row.object) for row in ms.itertuples(index=False)};first=next(iter(objs.values()));genes=first['genes'].astype(str);groups=first['groups'].astype(str);gix={g:i for i,g in enumerate(genes)};qix={g:i for i,g in enumerate(groups)}
  for v in objs.values():assert np.array_equal(v['genes'].astype(str),genes)
  groupidx={name:[gix[g] for g in members if g in gix] for name,members in sets.items()}
  for normalization in ['raw','corrected']:
   for scope in scopes:
    for reference in ['fixed_background','osteogenic_reference']:
     diffs=[];donors=[];pnames=[]
     for donor,dm in ms.groupby('donor_record',sort=True):
      ds=[];included=[]
      for site,sm in dm.groupby('site',sort=True,dropna=False):
       sid=sm.sample_id.tolist();left=qix[scope];right=qix[reference]
       lc=sum((objs[s][normalization][left] for s in sid),start=np.zeros(len(genes)));rc=sum((objs[s][normalization][right] for s in sid),start=np.zeros(len(genes)))
       ll=sum(float(objs[s][normalization+'_library'][left]) for s in sid);rl=sum(float(objs[s][normalization+'_library'][right]) for s in sid);ln=sum(int(objs[s]['cells'][left]) for s in sid);rn=sum(int(objs[s]['cells'][right]) for s in sid)
       if not ln or not rn or ll<=0 or rl<=0:continue
       effect=np.log2(lc/ll*1e4+.1)-np.log2(rc/rl*1e4+.1);ds.append(effect);included.append(str(site));siteledger.append({'source_id':source,'dataset':ms.dataset.iloc[0],'scope':scope,'reference':reference,'normalization':normalization,'donor_record':donor,'patient_id':dm.patient_id.iloc[0],'site':site,'target_cells':ln,'reference_cells':rn,'target_library':ll,'reference_library':rl,'sample_ids':';'.join(sid)})
      if ds:diffs.append(np.mean(ds,axis=0));donors.append(donor);pnames.append(dm.patient_id.iloc[0])
     n=len(diffs);ctx=f'{source}__{scope}__{reference}__{normalization}';contexts.append({'context_id':ctx,'source_id':source,'dataset':ms.dataset.iloc[0],'scope':scope,'reference':reference,'normalization':normalization,'N':n,'genes':len(genes),'state':'inference_available' if n>=2 else ('direction_only' if n==1 else 'no_matched_patients')})
     if not n:continue
     d=np.vstack(diffs);np.savez_compressed(O/'objects'/f'{ctx}_gene_effects.npz',effects=d,genes=genes,donors=np.asarray(donors),patients=np.asarray(pnames))
     score=np.column_stack([d[:,ix].mean(1) if ix else np.full(n,np.nan) for ix in groupidx.values()]);valid=np.isfinite(score).all(0);pv=np.full(len(programs),np.nan);method=np.full(len(programs),'not_tested_no_measured_genes',dtype=object)
     if valid.any():pv[valid],method[valid]=tests(score[:,valid])
     fdr=bh(pv);np.savez_compressed(O/'objects'/f'{ctx}_program_effects.npz',effects=score,programs=np.asarray(programs),donors=np.asarray(donors),patients=np.asarray(pnames))
     for j,name in enumerate(programs):
      v=score[:,j];good=np.isfinite(v);mn=float(v[good].mean()) if good.any() else np.nan;loo=(v.sum()-v)/(n-1) if n>=2 and good.all() else np.full(n,np.nan)
      rows.append({'context_id':ctx,'source_id':source,'dataset':ms.dataset.iloc[0],'scope':scope,'reference':reference,'normalization':normalization,'program':name,'collection':gs['origin'][name],'defined_genes':len(sets[name]),'measured_genes':len(groupidx[name]),'N':int(good.sum()),'mean_effect':mn,'median_effect':float(np.median(v[good])) if good.any() else np.nan,'positive_patients':int((v>0).sum()),'negative_patients':int((v<0).sum()),'P':pv[j],'FDR':fdr[j],'nominal':bool(pv[j]<.05),'FDR_significant':bool(fdr[j]<.05),'test':method[j],'family':ctx+'|programs','LOO_same_mean_direction_fraction':float(np.mean(np.sign(loo)==np.sign(mn))) if np.isfinite(loo).all() else np.nan,'LOO_min_mean':float(loo.min()) if np.isfinite(loo).all() else np.nan,'LOO_max_mean':float(loo.max()) if np.isfinite(loo).all() else np.nan,'evidence_level':'exploratory_small_N' if n<10 else 'source_scoped_discovery'})
      for i,donor in enumerate(donors):patientrows.append({'context_id':ctx,'source_id':source,'scope':scope,'reference':reference,'normalization':normalization,'donor_record':donor,'patient_id':pnames[i],'program':name,'effect':v[i]})
     print(json.dumps({'context':ctx,'N':n,'nominal_programs':int((pv<.05).sum())}),flush=True)
   # Direct transcript contribution of historical target genes; all sample features enter the library denominator.
   for donor,dm in ms.groupby('donor_record',sort=True):
    for site,sm in dm.groupby('site',dropna=False):
     for scope in scopes:
      for gene in ['COL1A1','LOX','MMP14','RUNX2','ANGPT1','SMAD3']:
       if gene not in gix:continue
       j=gix[gene];sids=sm.sample_id.tolist();target=sum(float(objs[s][normalization][qix[scope],j]) for s in sids);total=sum(float(objs[s][normalization][qix['all_input'],j]) for s in sids)
       contributions.append({'source_id':source,'dataset':ms.dataset.iloc[0],'donor_record':donor,'patient_id':dm.patient_id.iloc[0],'site':site,'scope':scope,'normalization':normalization,'gene':gene,'candidate_counts':target,'all_input_counts':total,'candidate_RNA_share':target/total if total>0 else np.nan})
  for v in objs.values():v.close()
 result=pd.DataFrame(rows);result.to_csv(O/'tables/all_patient_program_statistics.tsv.gz',sep='\t',index=False);result[result.nominal].to_csv(O/'tables/nominal_patient_programs.tsv',sep='\t',index=False);pd.DataFrame(patientrows).to_csv(O/'tables/patient_program_effects.tsv.gz',sep='\t',index=False);pd.DataFrame(siteledger).to_csv(O/'tables/matched_patient_site_roster.tsv',sep='\t',index=False);pd.DataFrame(contexts).to_csv(O/'tables/executed_contexts.tsv',sep='\t',index=False);pd.DataFrame(contributions).to_csv(O/'tables/patient_site_target_gene_RNA_shares.tsv',sep='\t',index=False)
 # Compare scope and correction changes on the exact same patient identifiers, without treating overlapping scopes as replicates.
 pats=pd.DataFrame(patientrows);sens=[]
 comparisons=[('agreed11','agreed11_no_flag','CNV_flag_sensitivity'),('agreed11','all15','candidate_expansion'),('all15','all15_no_flag','expanded_CNV_flag_sensitivity')]
 for (source,ref,norm,prog),g in pats.groupby(['source_id','reference','normalization','program'],sort=False):
  for baseline,other,kind in comparisons:
   a=g[g.scope.eq(baseline)].set_index('donor_record').effect;b=g[g.scope.eq(other)].set_index('donor_record').effect;shared=a.index.intersection(b.index)
   if not len(shared):continue
   aa=a.loc[shared].to_numpy();bb=b.loc[shared].to_numpy();sens.append({'source_id':source,'reference':ref,'normalization':norm,'program':prog,'comparison':kind,'baseline_scope':baseline,'other_scope':other,'N_common':len(shared),'baseline_common_mean':np.mean(aa),'other_common_mean':np.mean(bb),'mean_effect_change':np.mean(bb-aa),'mean_direction_preserved':bool(np.sign(np.mean(aa))==np.sign(np.mean(bb))),'patient_direction_agreement':np.mean(np.sign(aa)==np.sign(bb))})
 pd.DataFrame(sens).to_csv(O/'tables/scope_sensitivity_common_patients.tsv',sep='\t',index=False)
 a=result[result.normalization.eq('raw')];b=result[result.normalization.eq('corrected')];cs=a.merge(b,on=['source_id','scope','reference','program'],suffixes=('_raw','_corrected'),validate='one_to_one');cs['direction_preserved']=np.sign(cs.mean_effect_raw)==np.sign(cs.mean_effect_corrected);cs.to_csv(O/'tables/correction_sensitivity.tsv.gz',sep='\t',index=False)
 # Patient equal-weight contribution summary, keeping sources separate.
 contrib=pd.DataFrame(contributions);patient=contrib.groupby(['source_id','dataset','donor_record','patient_id','scope','normalization','gene'],as_index=False).agg(candidate_RNA_share=('candidate_RNA_share','mean'));patient.to_csv(O/'tables/patient_target_gene_RNA_shares.tsv',sep='\t',index=False)
 csu=patient.groupby(['source_id','dataset','scope','normalization','gene'],as_index=False).agg(patient_records=('donor_record','nunique'),median_RNA_share=('candidate_RNA_share','median'),mean_RNA_share=('candidate_RNA_share','mean'));csu.to_csv(O/'tables/source_target_gene_RNA_shares.tsv',sep='\t',index=False)
 qa={'at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'state':'patient_programs_and_scope_sensitivity_complete','planned_contexts':len(contexts),'evaluable_contexts':sum(x['N']>0 for x in contexts),'programs':len(programs),'rows':len(result),'valid_tests':int(result.P.notna().sum()),'nominal_tests':int(result.nominal.sum()),'FDR_tests':int(result.FDR_significant.sum()),'patients_not_cells_are_inference_units':True,'scope_overlap_not_independent_replication':True,'no_cross_source_patient_pooling':True};dump(O/'delivery/PROGRAM_ACCEPTANCE.json',qa);(O/'PROGRAMS_COMPLETE').write_text('SUCCESS\n');print(json.dumps(qa),flush=True)
if __name__=='__main__':main()