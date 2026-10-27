from pathlib import Path
import pandas as pd,json
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF')
def rd(f,**kw):return pd.read_csv(f,sep='\t',low_memory=False,**kw)
for name in ['rerun_reference_control_behavior.tsv','rerun_candidate_evidence_summary.tsv','job_manifest.tsv']:
    f=P/'phase34_reference_rerun_20260927/tables'/name;print('\nRERUN',name);print(rd(f).to_string(index=False,max_colwidth=40))
f=P/'phase26_client_reanalysis_20260922/code/common.py';lines=f.read_text().splitlines();print('\nPARSER')
hits=[i for i,s in enumerate(lines) if any(k in s for k in ['infercnv_call','state','pred_cnv','alteration_called','no_alteration_called','hmm','reference'])]
idx=sorted(set(j for i in hits for j in range(max(0,i-3),min(len(lines),i+5))))
print('\n'.join(f'{i+1}: {lines[i]}' for i in idx))
for sm in ['GSE152048__BC17','SCPCS000531','SCPCS000532']:
    d=rd(P/'phase27_identity_review_20260926/tables'/f'{sm}_review.tsv.gz')
    ref=d[d.reference.astype(str).str.lower().eq('true')]
    print('\nREFERENCE DETAIL',sm,'n',len(ref));print(ref.groupby(['provisional_lineage','copykat_call','infercnv_call']).size().to_string())
    p=P/'phase26_client_reanalysis_20260922/cnv'/sm
    print('CNV_DIR',sm, [x.name for x in p.iterdir()] if p.exists() else 'ABSENT')
    if p.exists():
        for x in p.glob('*/*'):
            if x.is_file() and any(k in x.name for k in ['pred_cnv','summary','calls','cell_groups','annotation']):print(x.relative_to(P),x.stat().st_size)
