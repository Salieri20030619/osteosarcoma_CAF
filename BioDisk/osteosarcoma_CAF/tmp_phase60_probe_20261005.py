from pathlib import Path
import pandas as pd,json
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');B=P/'phase58_F1F5_round1_20261004';Q=P/'phase59_four_functional_axes_20261004'
def read(p,**kw):return pd.read_csv(p,sep='\t',low_memory=False,**kw)
o=read(Q/'tables/frozen_cells_with_identity.tsv.gz');f=o[o.class_id=='F3'];qc=read(B/'tables/cell_QC_and_class.tsv.gz');f=f.merge(qc[['cell_id','UMI','detected_features','mitochondrial_fraction','ribosomal_fraction']],on='cell_id',validate='one_to_one')
print('COUNTS',f.groupby(['dataset','sample','donor_record','site','round1_state','lineage_stratum'],observed=True).size().to_string())
print('QC',f.groupby(['donor_record','round1_state'],observed=True).agg(cells=('cell_id','size'),UMI=('UMI','median'),features=('detected_features','median'),MT=('mitochondrial_fraction','median'),Ribo=('ribosomal_fraction','median'),contamination=('decontx_contamination','median')).to_string())
print('FLAGS',f[['CNV_evidence_current','mixed_program_flag','doublet_assessed_any','doublet_positive_any']].value_counts(dropna=False).to_string())
print('STATECOLS',read(P/'phase42_F3F4_identity_states_20261001/tables/expanded_F3_internal_cell_states.tsv.gz').columns.tolist())
p=read(B/'tables/corrected_program_statistics.tsv.gz');p=p[(p.dataset=='SCPCP000017')&(p.class_id=='F3')];print('EXISTING_PROGRAM',p[p.program.str.contains('customer|HALLMARK',regex=True)][['program','N','mean_delta','positive_patients','P','FDR']].to_string(index=False))
libs=json.loads((B/'resources/libraries_canonical.json').read_text());print('TRANSLATION',[k for k in libs['reactome'] if any(z in k for z in ['TRANSLATION','RIBOSOMAL','MITOCHONDRIAL_TRANSLATION'])][:50])
