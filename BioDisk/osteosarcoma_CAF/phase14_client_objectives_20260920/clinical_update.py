from pathlib import Path
import json,csv,zipfile,io,xml.etree.ElementTree as E
import pandas as pd,numpy as np
O=Path(__file__).resolve().parent;T=O/'tables'
ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
def doc_rows(archive,name):
    with zipfile.ZipFile(O/archive) as z:
        with zipfile.ZipFile(io.BytesIO(z.read(name))) as d:r=E.fromstring(d.read('word/document.xml'))
    return [[''.join(t.text or '' for t in c.findall('.//w:t',ns)) for c in row.findall('w:tc',ns)] for row in r.findall('.//w:tr',ns)]
old=doc_rows('GSE162454_original_supplementary.zip','Table_1.docx');new=doc_rows('SEC16B_supplementary.zip','FBA2-7-e70025-s003.docx')
base={r[0]:r for r in old[1:]};upd={r[1]:r for r in new if len(r)>5 and r[1] in base}
assert len(base)==len(upd)==6
gsm=['GSM4952363','GSM4952364','GSM4952365','GSM5155198','GSM5155199','GSM5155200'];mapped=[]
for i,sample in enumerate(gsm,1):
    pid='OS'+str(i);a=base[pid];b=upd[pid];assert a[1:3]==b[2:4]
    mapped.append(dict(sample='GSE162454__'+sample,paper_patient_id=pid,sex=b[2],age=b[3],location=b[4],metastasis_any=b[5],neoadjuvant_chemotherapy=a[3],mapping_evidence='matching OS1-6 IDs and all six sex/age profiles across original and followup papers; GEO Patient1-6 numbered mapping',mapping_level='cross_publication_reconciled_not_new_author_confirmation',metastatic_organ='not_specified_in_TableS1',metastasis_timing='not_specified_in_TableS1'))
m=pd.DataFrame(mapped);m.to_csv(T/'01_GSE162454_cross_publication_clinical.tsv',sep='\t',index=False)
c=pd.read_csv(T/'01_clinical_sample_audit.tsv',sep='\t');c=c.merge(m,on='sample',how='left');c['updated_comparison_role']=c.clinical_group
c.loc[c.metastasis_any.eq('No'),'updated_comparison_role']='A_candidate_reported_no_metastasis_timing_unknown'
c.loc[c.metastasis_any.eq('Yes'),'updated_comparison_role']='metastatic_bone_primary_lung_and_timing_unconfirmed'
c.to_csv(T/'01_clinical_sample_audit_updated.tsv',sep='\t',index=False)
e=pd.read_csv(T/'02_patient_gene_expression.tsv.gz',sep='\t');case='GSE162454__GSM4952365';control='GSE162454__GSM5155200'
a=e[e['sample'].eq(case)].set_index(['set','expression','gene']);b=e[e['sample'].eq(control)].set_index(['set','expression','gene']);joined=a[['cells','CP10k','detection_fraction']].join(b[['cells','CP10k','detection_fraction']],lsuffix='_metastasis',rsuffix='_no_metastasis')
joined['delta_log2CP10k']=np.log2(joined.CP10k_metastasis+.1)-np.log2(joined.CP10k_no_metastasis+.1);joined['detection_difference']=joined.detection_fraction_metastasis-joined.detection_fraction_no_metastasis
joined['interpretation']='one patient each; cross-publication any-metastasis labels; not verified lung-specific B1; descriptive only'
joined.reset_index().to_csv(T/'02_GSE162454_any_metastasis_descriptive.tsv',sep='\t',index=False)
d=pd.read_csv(T/'03_surface_screen_donor.tsv.gz',sep='\t');a=d[d.donor_record.eq('GSE162454|GSE162454__Patient_3')].set_index(['markers','count_threshold']);b=d[d.donor_record.eq('GSE162454|GSE162454__Patient_6')].set_index(['markers','count_threshold']);keep=['target_cells','target_coverage','known_target_fraction','known_CAF_fraction','CNV_risk_fraction','mural_fraction','gate_cells'];q=a[keep].join(b[keep],lsuffix='_metastasis',rsuffix='_no_metastasis');q['target_coverage_difference']=q.target_coverage_metastasis-q.target_coverage_no_metastasis;q['interpretation']='preparatory RNA contrast, one patient each; neither protein purity nor metastasis-specific target established';q.reset_index().to_csv(T/'03_GSE162454_marker_pair_contrast.tsv',sep='\t',index=False)
print(c[c.dataset.eq('GSE162454')][['sample','metastasis_any','unified_strict','updated_comparison_role']].to_string(index=False))
print(joined.reset_index().query("set=='customer_myCAF' and expression=='corrected'").sort_values('delta_log2CP10k',ascending=False).head(12).to_string(index=False))
print(q.reset_index().query('count_threshold==1 and target_coverage_metastasis>=0.2').sort_values('target_coverage_difference',ascending=False).head(8).to_string(index=False))
