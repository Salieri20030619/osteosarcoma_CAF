from pathlib import Path
import pandas as pd,json
p=Path.cwd();out=p/'remediation_20260917/pair_attrition';out.mkdir(exist_ok=True)
a=pd.read_csv(p/'remediation_20260917/state_final_sensitivity/corrected_recluster_cells.tsv.gz',sep='\t',dtype={'harmony_theta2_r0.5':str,'corrected_state':str});rows=[];flows=[];grid=[]
for old,new in [('0','0'),('2','1')]:
 for pat,d in a.groupby('patient_id'):
  o=d['harmony_theta2_r0.5']==old;n=d.corrected_state==new;ot=int(o.sum());nt=int(n.sum());N=len(d);be=ot>=20 and N-ot>=20;af=nt>=20 and N-nt>=20
  reason='remains_eligible' if af else ('insufficient_total_cells_for_two_sides' if N<40 else ('target_below20' if nt<20 else 'other_CAF_below20'))
  rows.append({'original_state':'H'+old,'corrected_cluster':new,'patient':pat,'total_CAF':N,'target_before':ot,'other_before':N-ot,'target_after':nt,'other_after':N-nt,'eligible_before':be,'eligible_after':af,'original_target_retained_in_match':int((o&n).sum()),'moved_out':int((o&~n).sum()),'moved_in':int((n&~o).sum()),'reason':reason})
  for (oc,nc),z in d.groupby(['harmony_theta2_r0.5','corrected_state']):
   if oc==old or nc==new:flows.append({'target_original_state':'H'+old,'patient':pat,'original_cluster':oc,'corrected_cluster':nc,'cells':len(z)})
 for k in [5,10,15,20,25,30]:
  q=[r for r in rows if r['original_state']=='H'+old];before=[r['patient'] for r in q if r['target_before']>=k and r['other_before']>=k];after=[r['patient'] for r in q if r['target_after']>=k and r['other_after']>=k]
  grid.append({'state':'H'+old,'minimum_per_side':k,'eligible_before':len(before),'eligible_after':len(after),'patients_after':';'.join(after),'interpretation':'diagnostic sensitivity only; not a new threshold selection'})
df=pd.DataFrame(rows);assert df[df.original_state=='H0'].eligible_before.sum()==4;assert df[df.original_state=='H2'].eligible_before.sum()==3;assert df.groupby('original_state').eligible_after.sum().eq(2).all()
df.to_csv(out/'patient_pair_attrition.tsv',sep='\t',index=False);pd.DataFrame(flows).to_csv(out/'cluster_cell_flows.tsv',sep='\t',index=False);pd.DataFrame(grid).to_csv(out/'threshold_diagnostic.tsv',sep='\t',index=False)
print(df.to_string(index=False));print(pd.DataFrame(grid).to_string(index=False));(out/'COMPLETE').write_text('PAIR_ATTRITION_AUDIT_COMPLETE\n')
