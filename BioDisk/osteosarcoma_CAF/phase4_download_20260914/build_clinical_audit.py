"""SYNTAX_CHECK: python -m py_compile build_clinical_audit.py
RUN_COMMAND: python build_clinical_audit.py
EXPECTED_SUCCESS_MARKER: CLINICAL_SOURCE_AUDIT_COMPLETE
"""
from pathlib import Path
import csv,json,openpyxl
P=Path(__file__).resolve().parent
rows=[]
w=openpyxl.load_workbook(P/'Supplementary Table S1 - Sample characteristics.xlsx',read_only=True,data_only=True)
for sample,source,sex,age,location,ncells in list(w.active.values)[1:7]:
 rows.append(dict(dataset='GSE299023',sample=sample,patient_id='GSE299023:'+sample,source_type_raw='Patient tumor biopsy',source_verified='direct_patient_tumor_non_PDX',sex_raw=sex,age_raw=age,tissue_location_raw=location,sampling_site_verified='Bone',disease_stage_raw='treatment-naive primary tumor (article)',preoperative_treatment_raw='untreated (article)',clinical_group_verified='Unassigned',lung_metastasis_history='not_reported',author_post_QC_cells=ncells,source='PMC12248827 Supplementary Table S1 and Methods 2.1',unresolved='lung_metastasis_history;sampling_date;treatment_timing',note='Author cell counts are not local QC; low-yield samples retained for CAF information assessment'))
w.close()
table=[('BC2','Primary','Femur','4 times (MTX, AP, IFO, MTX)'),('BC3','Primary','Tibia','6 times (MTX, AP, MTX, AP, MTX, MTX)'),('BC5','Primary','Femur','3 times (MTX, IFO, AP)'),('BC6','Primary','Ulna','3 times (MTX, IFO, AP)'),('BC10','Metastasis (Lung)','Femur','2 times (GT, GT)'),('BC11','Recurrent','Femur','3 times (GT, GT, GT)'),('BC16','Primary','Tibia','4 times (IFO, AP, MTX, MTX)'),('BC17','Metastasis (Lung)','Tibia','3 times (GT, GT, GT)'),('BC20','Recurrent','Femur','4 times (AP, IE, AP, IE)'),('BC21','Primary','Femur','4 times (MTX, AP, IFO, AP)'),('BC22','Primary','Femur','4 times (MTX, AP, IFO, MTX)')]
for sample,stage,location,treatment in table:
 lung=stage=='Metastasis (Lung)'
 rows.append(dict(dataset='GSE152048',sample=sample,patient_id='GSE152048:'+sample,source_type_raw='Surgically collected patient tumor',source_verified='direct_patient_tumor_non_PDX',sex_raw='not_extracted',age_raw='not_extracted',tissue_location_raw=location,sampling_site_verified='Lung' if lung else location,disease_stage_raw=stage,preoperative_treatment_raw=treatment,clinical_group_verified='B2' if lung else 'Unassigned',lung_metastasis_history='Yes' if lung else 'not_reported',author_post_QC_cells='',source='PMC7730477 Supplementary Table 1, PDF page 25; Methods and Results',unresolved='sampling_date;complete_treatment_timing'+('' if lung else ';lung_metastasis_history'),note=('Table Location is bone despite lung-metastasis specimen; preserve raw value and use explicit article lung sampling evidence. ' if lung else '')+('BC17 additionally received anti-PD-L1 according to Methods; do not treat as chemo-only. ' if sample=='BC17' else '')))
with (P/'clinical_source_audit.tsv').open('w',encoding='utf-8-sig',newline='') as f:
 v=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');v.writeheader();v.writerows(rows)
assert len(rows)==17 and sum(r['clinical_group_verified']=='B2' for r in rows)==2
print('CLINICAL_SOURCE_AUDIT_COMPLETE',len(rows))
