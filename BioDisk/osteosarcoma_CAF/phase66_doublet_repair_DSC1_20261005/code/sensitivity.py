from pathlib import Path
import os,sys,json,datetime
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
import numpy as np,pandas as pd
from current_metadata import load_current_cells
B=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=B/'phase66_doublet_repair_DSC1_20261005'
sys.path.insert(0,str(B/'phase58_F1F5_round1_20261004/code'));from common import tests,fdr,mean_valid
ob=load_current_cells().set_index('cell_id');z=np.load(B/'phase63_matrix_resistance_20261005/objects/mechanism_counts.npz')
assert list(z['cells'])==list(ob.index);genes=z['genes'];gi={g:i for i,g in enumerate(genes)}
cfg=json.loads((O/'resources/PROGRAMS.json').read_text());extra=json.loads((O/'resources/phase65_protocol.json').read_text())
panels=[(p,'full',gs,cfg['roles'][p]) for p,gs in cfg['sets'].items()]
panels += [(p,'full',gs,'posthoc_structure_panels') for p,gs in extra['structure'].items()]
panels += [(p,'drop_dominant_pair',[g for g in gs if g not in extra['drop_dominant'][p]],'reduced_structure_panels') for p,gs in extra['structure'].items()]
mem=[np.array([gi[g] for g in gs]) for _,_,gs,_ in panels]
stats=[];geneout=[];effects=[];rosters=[]
for arm in ['keep_all','exclude_model_positive']:
 keep=np.ones(len(ob),bool) if arm=='keep_all' else ~ob.doublet_positive_any.to_numpy(bool)
 for ds in sorted(ob.dataset.unique()):
  source=ob.dataset.eq(ds).to_numpy()
  donors=list(ob.loc[source].groupby('donor_record',sort=True))
  for st in sorted(ob.loc[source,'subgroup_id'].unique()):
   for mode,scope,minimum in [('raw','all',1),('corrected','all',1),('corrected','all',5),('corrected','Fibro',1)]:
    l=source&ob.subgroup_id.eq(st).to_numpy();r=source&ob.subgroup_id.ne(st).to_numpy()
    if scope=='Fibro':fib=ob.lineage_stratum.eq('Fibro').to_numpy();l &=fib;r &=fib
    if not l.any() or not r.any():continue
    l &= keep;r &= keep;x=z[mode];lib=z[mode+'_library'];gd=[];pdiff=[];names=[]
    for donor,dg in donors:
     sd=[];sp=[]
     for site,sg in dg.groupby('site',dropna=False):
      ix=ob.index.get_indexer(sg.index);il=ix[l[ix]];ir=ix[r[ix]]
      if min(len(il),len(ir))<minimum:continue
      present=z['assayed'][np.r_[il,ir]].all(0)
      dv=np.log2(x[il].sum(0)/lib[il].sum()*1e4+.1)-np.log2(x[ir].sum(0)/lib[ir].sum()*1e4+.1)
      dv[~present]=np.nan;sd.append(dv);pv=[]
      for jj in mem:
       have=jj[present[jj]];ok=len(have)>=3 and len(have)/len(jj)>=.8
       pv.append(float(dv[have].mean()) if ok else np.nan)
      sp.append(pv);rosters.append({'arm':arm,'dataset':ds,'subgroup_id':st,'matrix':mode,'scope':scope,'minimum_cells':minimum,'donor_record':donor,'site':site,'target_cells':len(il),'reference_cells':len(ir)})
     if sd:gd.append(mean_valid(np.array(sd),axis=0));pdiff.append(mean_valid(np.array(sp),axis=0));names.append(donor)
    gd=np.array(gd).reshape((-1,len(genes)));pp=np.array(pdiff).reshape((-1,len(panels)))
    pv,nn,method=tests(pp)
    for j,(p,version,gs,role) in enumerate(panels):
     v=pp[:,j];v=v[np.isfinite(v)]
     stats.append({'arm':arm,'dataset':ds,'subgroup_id':st,'matrix':mode,'scope':scope,'minimum_cells':minimum,'program':p,'version':version,'N':int(nn[j]),'mean_effect':float(v.mean()) if len(v) else np.nan,'positive_patients':int(sum(v>0)),'negative_patients':int(sum(v<0)),'P':pv[j],'test':method[j],'family':f'{arm}|{ds}|{mode}|{scope}|min{minimum}|{role}'})
     for donor,value in zip(names,pp[:,j]):effects.append({'arm':arm,'dataset':ds,'subgroup_id':st,'matrix':mode,'scope':scope,'minimum_cells':minimum,'program':p,'version':version,'donor_record':donor,'effect':value})
    if scope=='all' and minimum==1:
     gp,gn,gm=tests(gd);geneout.append(pd.DataFrame({'arm':arm,'dataset':ds,'subgroup_id':st,'matrix':mode,'gene':genes,'N':gn,'mean_effect':mean_valid(gd,axis=0),'positive_patients':(gd>0).sum(0),'P':gp,'test':gm,'family':f'{arm}|{ds}|{mode}|all_subgroup_targeted_mechanism_genes'}))
  print('DONE',arm,ds,flush=True)
s=fdr(pd.DataFrame(stats));s.to_csv(O/'tables/program_statistics.tsv.gz',sep='\t',index=False)
g=fdr(pd.concat(geneout,ignore_index=True));g.to_csv(O/'tables/gene_statistics.tsv.gz',sep='\t',index=False)
pd.DataFrame(effects).to_csv(O/'tables/patient_program_effects.tsv.gz',sep='\t',index=False);pd.DataFrame(rosters).to_csv(O/'tables/paired_roster.tsv.gz',sep='\t',index=False)
key=['dataset','subgroup_id','matrix','scope','minimum_cells','program','version']
w=s[s.arm=='keep_all'].merge(s[s.arm=='exclude_model_positive'],on=key,suffixes=('_all','_excluded'),validate='one_to_one')
w['effect_change']=w.mean_effect_excluded-w.mean_effect_all
w['direction_flip']=(w.mean_effect_all*w.mean_effect_excluded<0)
w['nominal_change']=w.P_all.lt(.05)!=w.P_excluded.lt(.05)
w['lost_evaluability']=w.P_all.notna()&w.P_excluded.isna()
w.to_csv(O/'tables/program_before_after_comparison.tsv',sep='\t',index=False)
gkey=['dataset','subgroup_id','matrix','gene']
gw=g[g.arm=='keep_all'].merge(g[g.arm=='exclude_model_positive'],on=gkey,suffixes=('_all','_excluded'),validate='one_to_one')
gw['effect_change']=gw.mean_effect_excluded-gw.mean_effect_all;gw['nominal_change']=gw.P_all.lt(.05)!=gw.P_excluded.lt(.05)
gw.to_csv(O/'tables/gene_before_after_comparison.tsv.gz',sep='\t',index=False)
(O/'SENSITIVITY_COMPLETE.json').write_text(json.dumps({'success':True,'program_rows':len(s),'gene_rows':len(g),'retained_cells':8148,'baseline_cells':8174,'at':datetime.datetime.now(datetime.timezone.utc).isoformat()},indent=2))
focus=w[(w.dataset=='SCPCP000017')&(w.matrix=='corrected')&(w.scope=='all')&(w.minimum_cells==1)&w.subgroup_id.isin(['F2_C1','F2_C3'])&w.program.isin(['fibrillar_I_III_V','collagen_IV','collagen_crosslinking'])]
print(focus[['subgroup_id','program','version','N_all','N_excluded','mean_effect_all','mean_effect_excluded','P_all','P_excluded','FDR_all','FDR_excluded']].to_string(index=False),flush=True)
