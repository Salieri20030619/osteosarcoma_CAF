from pathlib import Path
import numpy as np,pandas as pd,json
p=Path.cwd();root=p/'remediation_20260917';out=root/'cnv_comparison';out.mkdir(exist_ok=True)
identity=pd.read_csv(p/'identity/all_cell_identity_evidence.tsv.gz',sep='\t',usecols=['cell_id','sample','caf_strict']).set_index('cell_id');summary=[];roles=[]
for short in ['BC10','BC17','BC2','BC5','BC21','BC22']:
 sample='GSE152048__'+short;alt=root/'cnv_immune_reference/infercnv'/sample
 if not (alt/'COMPLETE').exists():summary.append({'sample':sample,'status':'failed_or_incomplete'});continue
 old=pd.read_csv(p/'infercnv'/sample/'cell_results.tsv',sep='\t');new=pd.read_csv(alt/'cell_results.tsv',sep='\t');frozen=pd.read_csv(p/'inputs'/sample/'cells.tsv',sep='\t').set_index('cell_id').cnv_role
 m=old[['cell_id','call','signal_strength','signal_cutoff']].merge(new[['cell_id','call','signal_strength','signal_cutoff']],on='cell_id',suffixes=('_baseline','_immune'),validate='one_to_one');assert len(m)==len(old)==len(new)
 m['frozen_role']=m.cell_id.map(frozen);m['strict_CAF']=m.cell_id.map(identity.caf_strict).fillna(False);ev=(m.call_baseline!='not_evaluable')&(m.call_immune!='not_evaluable');strict=m[m.strict_CAF]
 m.to_csv(out/(sample+'_cell_comparison.tsv.gz'),sep='\t',index=False)
 for role,d in m[m.frozen_role.str.startswith(('ref_','holdout_'))].groupby('frozen_role'):
  e=d[d.call_immune!='not_evaluable'];roles.append({'sample':sample,'frozen_role':role,'evaluable':len(e),'normal_fraction':float((e.call_immune=='low_signal').mean()) if len(e) else None})
 summary.append({'sample':sample,'status':'complete','jointly_evaluable':int(ev.sum()),'same_call_fraction':float((m.loc[ev,'call_baseline']==m.loc[ev,'call_immune']).mean()),'baseline_strict_CAF':len(strict),'strict_still_low':int((strict.call_immune=='low_signal').sum()),'strict_now_elevated':int((strict.call_immune=='elevated_signal').sum()),'strict_now_ambiguous':int((strict.call_immune=='ambiguous_signal').sum()),'strict_not_evaluable':int((strict.call_immune=='not_evaluable').sum()),'baseline_cutoff':old.signal_cutoff.iloc[0],'immune_cutoff':new.signal_cutoff.iloc[0]})
pd.DataFrame(summary).to_csv(out/'summary.tsv',sep='\t',index=False);pd.DataFrame(roles).to_csv(out/'frozen_reference_strata.tsv',sep='\t',index=False)
(out/'COMPLETE').write_text('CNV_REFERENCE_COMPARISON_COMPLETE\n');print(pd.DataFrame(summary).to_string(index=False))
