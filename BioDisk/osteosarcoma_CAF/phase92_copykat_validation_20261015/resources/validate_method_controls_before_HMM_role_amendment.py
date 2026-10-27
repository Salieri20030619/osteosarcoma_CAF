from pathlib import Path
import sys,json,re,datetime
import pandas as pd
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase92_copykat_validation_20261015'
sm=sys.argv[1]
receipt=json.loads((O/'resources'/(sm+'_COMPLETE.json')).read_text(encoding='utf-8'))
assert receipt['exit_code']==0 and receipt['completion_marker'] and (O/'copykat'/sm/'COMPLETE').exists()
c=pd.read_csv(O/'tables'/(sm+'_copykat_full_evidence.tsv.gz'),sep='\t',encoding='utf-8').set_index('cell_id')
assert c.index.is_unique and not (c.reference_used&c.heldout).any()
assert c.patient_id.nunique()==1 and set(c['sample'])=={sm}
log=(O/'logs'/(sm+'_copykat.log')).read_text(encoding='utf-8')
actual_refs=re.findall(r'(\d+) known normal cells found in dataset',log)
assert len(actual_refs)==1,'Native reference count must be explicit and unambiguous'
actual_refs=int(actual_refs[0]);checks=[]
# Include every originally intended heldout stratum, even if all its cells became doublet-positive.
for role,g0 in c[c.heldout].groupby('reference_role'):
 g=g0[g0.heldout_clean];ev=g.native_copykat_call.isin(['diploid','aneuploid']);n=int(ev.sum());norm=int(g.native_copykat_call.eq('diploid').sum())
 checks.append({'reference_role':role,'original_heldouts':len(g0),'new_doublets_excluded':int((~g0.heldout_clean).sum()),'clean_heldouts':len(g),'evaluable':n,'diploid':norm,'not_evaluable':len(g)-n,'diploid_fraction':norm/n if n else None,'pass':bool(n>=10 and norm/n>=.9) if n else False})
copy_ok=actual_refs>=30 and bool(checks) and sum(x['clean_heldouts'] for x in checks)>=20 and all(x['pass'] for x in checks)
hmmdir=(P/'phase91_lung_patient_validation_20261015/cnv/infercnv'/sm) if sm.startswith('GSE') else (P/'phase90_patient_validation_20261015/cnv/infercnv'/sm)
if sm=='NCC__BS7' and (P/'phase90_patient_validation_20261015/cnv_clean/infercnv'/sm/'COMPLETE').exists():hmmdir=P/'phase90_patient_validation_20261015/cnv_clean/infercnv'/sm
hmm={'status':'pending','valid_for_binary_identity':False}
hmm_exit=(P/'phase91_lung_patient_validation_20261015/logs/infercnv.exit_code') if sm.startswith('GSE') else (P/'phase90_patient_validation_20261015/logs'/(sm+'_infercnv.exit_code'))
hmm_exit_ok=hmm_exit.exists() and hmm_exit.read_text(encoding='utf-8').strip()=='0'
if sm=='NCC__BS7' and 'cnv_clean' in str(hmmdir):
 rr=P/'phase90_patient_validation_20261015/resources/BS7_CLEAN_REFERENCE_EXIT.json'
 hmm_exit_ok=rr.exists() and json.loads(rr.read_text(encoding='utf-8')).get('exit_code')==0
if (hmmdir/'COMPLETE').exists() and hmm_exit_ok:
 h=pd.read_csv(hmmdir/'cell_results.tsv',sep='\t',encoding='utf-8').set_index('cell_id');assert h.index.is_unique and set(h.index)==set(c.index)
 cc=c.copy();cc['hmm_call']=h['call'];hc=[]
 for role,g0 in cc[cc.heldout].groupby('reference_role'):
  g=g0[g0.heldout_clean];ev=g.hmm_call.isin(['no_alteration_called','alteration_called']);n=int(ev.sum());norm=int(g.hmm_call.eq('no_alteration_called').sum())
  hc.append({'reference_role':role,'original_heldouts':len(g0),'clean_heldouts':len(g),'evaluable':n,'no_alteration':norm,'no_alteration_fraction':norm/n if n else None,'pass':bool(n>=10 and norm/n>=.9) if n else False})
 ok=bool(hc) and all(x['pass'] for x in hc)
 hmm={'status':'complete_reference_pass' if ok else 'complete_reference_failed','valid_for_binary_identity':ok,'checks':hc,'file':str(hmmdir/'cell_results.tsv'),'raw_calls_preserved':True,'failed_is_not_normal':True}
summary={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'sample':sm,'patient_id':str(c.patient_id.iloc[0]),'full_cells':len(c),'CopyKAT_input_references':int(c.reference_used.sum()),'CopyKAT_software_actual_references':actual_refs,'CopyKAT_reference_checks':checks,'CopyKAT_valid_for_binary_identity':copy_ok,'CopyKAT_calls_including_missing':c.native_copykat_call.fillna('native_prediction_missing').value_counts().to_dict(),'HMM':hmm,'source_patient_N':1,'new_formal_patients':0,'scope':'Method validity only; lineage/ambient/patient admission remain separate. Intended-normal controls are not genomic ground truth or independent biological replicates. Aggregate patients by author-resolved patient_id, not by completed library count.'}
assert sum(summary['CopyKAT_calls_including_missing'].values())==len(c)
(O/'resources'/(sm+'_VALIDITY_RECHECK.json')).write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
print(json.dumps(summary,indent=2))
