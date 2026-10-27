from pathlib import Path
import pandas as pd
ROOT=Path(__file__).resolve().parents[2]
EVIDENCE=ROOT/'phase66_doublet_repair_DSC1_20261005/tables/current_doublet_evidence.tsv.gz'
CELLS=ROOT/'phase66_doublet_repair_DSC1_20261005/tables/cells_repaired_with_confidence.tsv.gz'
def load_current_cells():
 d=pd.read_csv(CELLS,sep='\t',low_memory=False)
 assert len(d)==8174 and d.cell_id.is_unique
 assert d.doublet_assessed_any.notna().all() and d.doublet_positive_any.notna().all()
 assert d.doublet_assessed_any.all() and d.doublet_positive_any.sum()==26
 return d
def apply_current_doublet_evidence(frame,id_column='cell_id'):
 e=pd.read_csv(EVIDENCE,sep='\t').set_index('cell_id')
 d=frame.copy();ids=d[id_column] if id_column in d else d.index.to_series()
 assert ids.is_unique,'Duplicate cell identifiers'
 missing=set(ids)-set(e.index)
 if missing:raise ValueError(f'No audited doublet evidence for {len(missing)} cells; do not fill missing as False')
 aligned=e.loc[ids]
 for c in aligned:
  if c in ['doublet_assessed_any','doublet_positive_any'] and c in d and 'legacy_'+c not in d:d['legacy_'+c]=d[c]
  d[c]=aligned[c].to_numpy()
 if 'conservative_fibro' in d:
  required={'lineage_stratum','mixed_program_flag','lineage_confidence'}
  if not required<=set(d):raise ValueError('Cannot update conservative_fibro without its defining fields')
  if 'legacy_conservative_fibro' not in d:d['legacy_conservative_fibro']=d['conservative_fibro']
  mixed=d.mixed_program_flag.astype(str).str.lower().isin(['true','1'])
  d['conservative_fibro']=d.lineage_stratum.eq('Fibro')&~mixed&~d.doublet_positive_any&d.lineage_confidence.isin(['高','中'])
 assert d.doublet_assessed_any.notna().all() and d.doublet_positive_any.notna().all()
 return d
