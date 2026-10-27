from pathlib import Path
import pandas as pd,json
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase66_doublet_repair_DSC1_20261005')
p=O/'tables/cells_repaired_with_confidence.tsv.gz';d=pd.read_csv(p,sep='\t',low_memory=False)
d['legacy_conservative_fibro']=d.conservative_fibro
mixed=d.mixed_program_flag.astype(str).str.lower().isin(['true','1'])
d['conservative_fibro']=d.lineage_stratum.eq('Fibro')&~mixed&~d.doublet_positive_any&d.lineage_confidence.isin(['高','中'])
d.to_csv(p,sep='\t',index=False)
changed=d[d.conservative_fibro!=d.legacy_conservative_fibro]
changed[['cell_id','subgroup_id','dataset','legacy_conservative_fibro','conservative_fibro','doublet_positive_any']].to_csv(O/'tables/derived_flag_changes.tsv',sep='\t',index=False)
helper=O/'code/current_metadata.py';s=helper.read_text()
s=s.replace("assert d.doublet_assessed_any.notna().all() and d.doublet_positive_any.notna().all()\n return d","""if 'conservative_fibro' in d:
  required={'lineage_stratum','mixed_program_flag','lineage_confidence'}
  if not required<=set(d):raise ValueError('Cannot update conservative_fibro without its defining fields')
  if 'legacy_conservative_fibro' not in d:d['legacy_conservative_fibro']=d['conservative_fibro']
  mixed=d.mixed_program_flag.astype(str).str.lower().isin(['true','1'])
  d['conservative_fibro']=d.lineage_stratum.eq('Fibro')&~mixed&~d.doublet_positive_any&d.lineage_confidence.isin(['高','中'])
 assert d.doublet_assessed_any.notna().all() and d.doublet_positive_any.notna().all()
 return d""")
helper.write_text(s)
audit={'root_cause':'phase55 freeze_annotation_phase55.py line31 LEFT JOIN of targeted phase54 review overlay: non-targeted cells inherit missing doublet metadata, even when source calls existed. phase59 copied these columns; phase62 checked preservation but not source coverage.','fix':'New full8174 source-derived overlay, strict cell-ID join, full assessed coverage assertion, legacy columns, caller details, current entry manifest; derived conservative_fibro recomputed from repaired calls.','derived_flag_changes':len(changed),'lineage_and_subgroup_confidence':'Preserved: phase55 lineage_confidence and candidate_tier derive from cluster review independent of this partial doublet overlay; model risk stored separately.','historical_frozen_files':'Not overwritten; successor version explicitly supersedes doublet fields.','sensitivity_impact_of_derived_fix':'No rerun needed: sensitivity uses direct repaired doublet flag and lineage_stratum, never conservative_fibro.'}
(O/'resources/ROOT_CAUSE_AND_FIX.json').write_text(json.dumps(audit,indent=2));print(json.dumps(audit),flush=True)
