
from pathlib import Path
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='2'
import pandas as pd,numpy as np,json,hashlib
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'patient_expansion_download_investigation_20261013_v1'
def read(f,**kw):return pd.read_csv(f,sep='\t',low_memory=False,**kw)
def truth(s):return s.astype(str).str.lower().isin(['true','1','1.0'])
m=read(P/'phase26_client_reanalysis_20260922/tables/analysis_source_manifest.tsv');clinical=read(P/'phase81_allCAF_bone_lung_20261012/tables/clinical_sample_manifest.tsv').set_index('sample');final=read(P/'phase73_F134_bone_lung_20261009/strict_v2/tables/expanded_fixed_cells_with_sites.tsv.gz');rows=[];reasons=[];candidate_cells=[];refs=[]
for r in m[m.dataset.eq('SCPCP000023')].itertuples():
 f=P/'phase27_identity_review_20260926/tables'/f'{r.sample}_review.tsv.gz';d=read(f,keep_default_na=False);fib=d.provisional_lineage.eq('Fibroblast');low=d.CNV_evidence.isin(['single_method_no_alteration','concordant_no_alteration']);doublet=truth(d.known_doublet);candidate=truth(d.CAF_candidate);assert (candidate==(fib&low&~doublet)).all()
 g=final[final['sample'].eq(r.sample)];cc=clinical.loc[r.sample]
 row=dict(sample=r.sample,donor_record=cc.donor_record,strict_site=cc.strict_site_group,full_called_cells=len(d),provisional_fibroblast_cells=int(fib.sum()),known_doublet_fibro=int((fib&doublet).sum()),stage27_CAF_candidates=int(candidate.sum()),final_candidates=len(g),final_Fibro=int(g.review_lineage_code.eq('Fibro').sum()),original_count_object_exists=Path(r.counts_path).exists())
 for state in sorted(d.CNV_evidence.unique()):row['fibro_'+state]=int((fib&~doublet&d.CNV_evidence.eq(state)).sum())
 rows.append(row)
 z=d[fib].groupby(['copykat_call','infercnv_call','CNV_evidence','review_class']).size().reset_index(name='cells');z['sample']=r.sample;z['donor_record']=cc.donor_record;z['strict_site']=cc.strict_site_group;reasons.append(z)
 re=d[truth(d.reference)];rr=re.groupby(['copykat_call','infercnv_call','CNV_evidence']).size().reset_index(name='cells');rr['sample']=r.sample;refs.append(rr)
 if candidate.any():candidate_cells.append(d[candidate].assign(strict_site=cc.strict_site_group))
ret=pd.DataFrame(rows);ret.to_csv(O/'tables/S23_all_samples_retention_waterfall.tsv',sep='\t',index=False);pd.concat(reasons).to_csv(O/'tables/S23_Fibro_CNV_exclusions.tsv',sep='\t',index=False);pd.concat(refs).to_csv(O/'tables/S23_reference_CNV_calls.tsv',sep='\t',index=False)
cs=pd.concat(candidate_cells);cs.to_csv(O/'tables/S23_stage27_candidate_cells.tsv.gz',sep='\t',index=False)
# Original downstream reviewed annotation for those cells, without relabeling.
cols=['original_cell_id','cell_id','sample','dataset','review_lineage_code','lineage_confidence','malignancy_assessment','phase54_review_decision','candidate_tier','CNV_evidence_current','confirmed_CAF_identity']
parts=[]
for chunk in read(P/'phase55_annotation_freeze_20261004/tables/all_cell_annotation_v2.tsv.gz',usecols=cols,chunksize=100000):
 chunk=chunk[chunk.dataset.eq('SCPCP000023')]
 if len(chunk):parts.append(chunk)
ann=pd.concat(parts);ann.to_csv(O/'tables/S23_frozen_downstream_annotation.tsv.gz',sep='\t',index=False)
print('S23_LUNG_WATERFALL');print(ret[ret.strict_site.eq('Lung')].fillna(0).to_string(index=False))
print('S23_ANNOTATIONS',ann.groupby(['review_lineage_code','candidate_tier']).size().to_string())
res={'status':'complete','samples':len(ret),'stage27_candidates':int(ret.stage27_CAF_candidates.sum()),'final_candidates':int(ret.final_candidates.sum()),'final_Fibro':int(ret.final_Fibro.sum()),'lung_stage27_patients':int(ret[(ret.strict_site=='Lung')&ret.stage27_CAF_candidates.gt(0)].donor_record.nunique()),'lung_final_Fibro_patients':int(ret[(ret.strict_site=='Lung')&ret.final_Fibro.gt(0)].donor_record.nunique()),'newly_admitted_patients':0,'upstream_gate_reproduced':True,'baseline_changed':False}
(O/'resources/S23_RETENTION_COMPLETE.json').write_text(json.dumps(res,indent=2));print(json.dumps(res))
