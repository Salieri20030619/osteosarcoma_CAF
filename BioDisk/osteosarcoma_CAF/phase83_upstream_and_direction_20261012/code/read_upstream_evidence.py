from pathlib import Path
import json
import pandas as pd
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF')
def read(p,**kw):return pd.read_csv(p,sep='\t',low_memory=False,**kw)
for n in ['phase27_identity_review_20260926/code/review.py','phase44_patient_mechanism_20261001/code/audit_anatomical_site.py','phase44_patient_mechanism_20261001/code/stageD_clinical_site_review.py']:
    f=P/n;print('\nCODE',n)
    lines=f.read_text().splitlines()
    if 'review.py'==f.name:print('\n'.join(f'{i+1}: {s}' for i,s in enumerate(lines)))
    else:
        for i,line in enumerate(lines):
            if any(s in line for s in ['Arm','Leg','Thigh','Unknown','UBERON','strict_site','bone','site=']):print(f'{i+1}: {line}')
for n in ['phase27_identity_review_20260926/tables/sample_review_counts.tsv','phase26_client_reanalysis_20260922/tables/all107_sample_outcomes.tsv']:
    d=read(P/n);print('\nTABLE',n,d.columns.tolist());print(d[d.iloc[:,0].astype(str).str.contains('BC17|BC2$|BC11|SCPCS000528|SCPCS000531|SCPCS000532',regex=True)].to_string(index=False))
for n in ['phase27_identity_review_20260926/tables/GSE152048__BC17_review.tsv.gz','phase26_client_reanalysis_20260922/tables/GSE152048__BC17_identity_evidence.tsv.gz']:
    d=read(P/n,nrows=3);print('\nSAMPLE SCHEMA',n,d.columns.tolist());print(d.to_string(index=False,max_colwidth=30))
f=P/'phase27_identity_review_20260926/resources/scpca_clinical_source.csv';d=pd.read_csv(f);print('\nCLINICAL',d.columns.tolist());mask=d.astype(str).apply(lambda c:c.str.contains('SARC0105|SARC0103|SARC0129|SCPCS000528|SCPCS000527|SCPCS000530',regex=True)).any(axis=1);print(d[mask].to_string(index=False))
print('\nSCPCP18 ORIGINAL METADATA FILES')
for f in (P/'phase25_complete_nonPDX_assets_20260921/files/SCPCP000018').rglob('*'):
    if f.is_file() and (f.suffix in ['.csv','.tsv','.json','.txt'] or '528' in str(f)):print(f.relative_to(P),f.stat().st_size)
