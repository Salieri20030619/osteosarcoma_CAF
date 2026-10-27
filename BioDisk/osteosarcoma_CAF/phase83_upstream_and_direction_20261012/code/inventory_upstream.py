from pathlib import Path
import json
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF')
for name in ['phase44_patient_mechanism_20261001/clinical_site_review/resources/CLINICAL_PROTOCOL_FROZEN.json','phase27_identity_review_20260926/PROTOCOL.json']:
    f=P/name;print('\nDOCUMENT',name);print(f.read_text() if f.exists() else 'ABSENT')
for phase in ['phase27_identity_review_20260926','phase26_client_reanalysis_20260922','phase25_complete_nonPDX_assets_20260921','phase15_scpca_download_20260920','phase44_patient_mechanism_20261001','phase28_CAF_clustering_20260926']:
    print('\nPHASE',phase)
    for sub in ['', 'tables','resources','code','delivery']:
        d=P/phase/sub
        if d.exists():print(sub, [f.name for f in d.iterdir()][:70])
print('\nPROJECT METADATA FILES')
for phase in ['phase25_complete_nonPDX_assets_20260921','phase15_scpca_download_20260920','phase16_scpca_supplement_20260920']:
    d=P/phase
    if d.exists():
        for pattern in ['*/*metadata*','*/*manifest*','*/*sample*','*/*project*','*/*/*metadata*','*/*/*manifest*']:
            for f in d.glob(pattern):print(f.relative_to(P),f.stat().st_size if f.is_file() else 'DIR')
