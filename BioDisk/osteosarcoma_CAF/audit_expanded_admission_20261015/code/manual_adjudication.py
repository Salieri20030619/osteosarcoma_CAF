from pathlib import Path
import pandas as pd, json, datetime
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF')
O=P/'audit_expanded_admission_20261015'
c=pd.read_csv(O/'tables/all_cells_review.tsv.gz',sep='\t',dtype={'cluster':str},low_memory=False)
cl=pd.read_csv(O/'tables/cluster_packets_for_manual_review.tsv',sep='\t',dtype={'cluster':str},keep_default_na=False)
competitors={'T_NK','B_Plasma','Myeloid','DC','Endo','Epi','Mast','Erythroid','Mural','Myogenic_pending','OC_like','Mixed'}
extra={
 ('NCC__BS2','0'):'T-lineage IL7R/LTB/CD40LG/CCR7',
 ('NCC__BS2','7'):'Cytotoxic T-lineage CD8A/CD8B/KLRK1/GZMB',
 ('NCC__BS2','8'):'Treg FOXP3/IL2RA/CTLA4/IKZF2',
 ('NCC__BS5','24'):'Lymphatic endothelial CCL21/PROX1/MMRN1',
 ('NCC__BS7','0'):'Skeletal myogenic DMD/RYR1/TRDN/MYOM3/MLIP; full top50 checked',
 ('NCC__BS7','4'):'Skeletal myogenic MYL1/TNNI2/TNNC2/TTN/TNNT3',
 ('NCC__BS14','7'):'Osteoclast RUFY4/SIGLEC15/ATP6V0D2'
}
matrix_note={
 ('NCC__BS1','5'):'Matrix-mural overlap COL3A1/COL12A1/COL6A3/FAP/TNC with EDNRA/PDGFRB; retain unresolved, do not equate with pure mural or confirmed CAF',
 ('NCC__BS2','5'):'DCN/LUM/FBLN1/COL6A1/COL12A1 fibro-matrix',
 ('NCC__BS4','1'):'Contractile-mesenchymal TPM2/MYL9/PRRX1/BGN/THY1/TAGLN; identity remains operational',
 ('NCC__BS5','20'):'SFRP2/ASPN/FNDC1/CXCL14/COL12A1 matrix state with osteochondral overlap',
 ('NCC__BS7','12'):'MFAP5/SFRP2/CFD/APOD/THBS4 stromal fibroblast program',
 ('NCC__BS7','22'):'GPC6/SFRP4/COL11A1/CDH11/POSTN/COL12A1 matrix program',
 ('NCC__BS9','1'):'DCN/CXCL5/MMP1/MMP3/IL6 inflammatory matrix, reference failure retained',
 ('NCC__BS9','4'):'PODN/FNDC1/OGN/ASPN/COL14A1/COL11A1 matrix, reference failure retained',
 ('NCC__BS12','0'):'COL14A1/COL8A1/C1R with TAGLN/NDUFA4L2 overlap; genomic context unresolved',
 ('NCC__BS14','2'):'POSTN/ADAMTSL1/AEBP1/VCAN/TAGLN matrix-contractile overlap',
 ('CXG__EXT_OS_2','16'):'FAP/ASPN/FNDC1/ADAMTS12/ITGBL1 matrix with genomic-context risk'
}
ann=[]
for r in cl[cl.native_entry.gt(0)].to_dict('records'):
 key=(r['sample'],str(r['cluster']))
 excluded=r['original_lineage'] in competitors
 reason='Explicit competing lineage supported by cluster marker packet: '+r['original_lineage'] if excluded else 'Retain identity-uncertain/osteogenic/cycling/stress state for expanded exploration; no confirmed CAF claim'
 if key in extra:excluded=True;reason=extra[key]
 if key in matrix_note:
  reason=matrix_note[key]
  if key==('NCC__BS1','5'):excluded=False
 r['manual_cluster_decision']='competing_lineage_background' if excluded else 'retain_expanded_exploratory_candidate'
 r['manual_reason']=reason
 r['matrix_program_highlight']=key in matrix_note
 r['review_evidence']='Existing cluster top15; top50 additionally inspected for BS7:0/19, BS1:5, BS5:20, BS4:1; per-cell Phase27 scores and CNV/ambient evidence'
 ann.append(r)
a=pd.DataFrame(ann)
a.to_csv(O/'tables/manual_cluster_adjudication.tsv',sep='\t',index=False,encoding='utf-8')
n=c[c.compatible_native_entry].merge(a[['sample','cluster','manual_cluster_decision','manual_reason','matrix_program_highlight']],on=['sample','cluster'],validate='many_to_one')
assert len(n)==2289 and n.manual_cluster_decision.notna().all()
n['expanded_review_retained']=n.manual_cluster_decision.eq('retain_expanded_exploratory_candidate')
n['reference_pass_retained']=n.expanded_review_retained & n.reference_QC_pass
n['reference_failed_retained']=n.expanded_review_retained & ~n.reference_QC_pass
n['formal_admission']=False
n.to_csv(O/'tables/final_cell_review_decisions.tsv.gz',sep='\t',index=False,encoding='utf-8')
s=n.groupby(['sample','dataset','patient_id'],sort=False).agg(native_entry=('cell_id','size'),retained_exploratory=('expanded_review_retained','sum'),reference_pass_retained=('reference_pass_retained','sum'),reference_failed_retained=('reference_failed_retained','sum')).reset_index()
s['competing_background']=s.native_entry-s.retained_exploratory
s['Bone_Lung_eligible_new_patients']=0
s.to_csv(O/'tables/final_sample_decisions.tsv',sep='\t',index=False,encoding='utf-8')
pt=s.groupby(['dataset','patient_id'],sort=False).agg(libraries=('sample','size'),retained_exploratory=('retained_exploratory','sum'),reference_pass_retained=('reference_pass_retained','sum'),reference_failed_retained=('reference_failed_retained','sum')).reset_index()
pt.to_csv(O/'tables/final_patient_capacity.tsv',sep='\t',index=False,encoding='utf-8')
tot={'review_scope':'Expanded layer prioritized; reviewed existing Phase88 retained cells; no source data or baseline modification','raw_expression_candidates':29923,'native_entry_audit':2289,'retained_exploratory':int(s.retained_exploratory.sum()),'competing_background':int(s.competing_background.sum()),'reference_pass_retained':int(s.reference_pass_retained.sum()),'reference_failed_retained':int(s.reference_failed_retained.sum()),'retained_source_donor_labels':int(pt.retained_exploratory.gt(0).sum()),'reference_pass_source_donor_labels':int(pt.reference_pass_retained.gt(0).sum()),'signal_only_native_unresolved':728,'retained_corrected_Fib_ge2':int((n.expanded_review_retained & n.corrected_det_Fib.ge(2)).sum()),'retained_corrected_strict_lineage_support':int((n.expanded_review_retained & n.corrected_lineage_supported).sum()),'strict179_overlap_with_retained':int((n.expanded_review_retained & n.supported_CAF_candidate).sum()),'strict_not_union_added':1,'reviewed_native_candidate_clusters':len(a),'new_formal_Bone_Lung_patients':0,'formal_baseline_unchanged':True,'time_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
(O/'REVIEW_COUNTS.json').write_text(json.dumps(tot,indent=2)+'\n',encoding='utf-8')
fig,ax=plt.subplots(figsize=(11,6.5));yy=range(len(s))
ax.barh(yy,s.reference_pass_retained,label='Exploratory retained / reference pass',color='#257c79')
ax.barh(yy,s.reference_failed_retained,left=s.reference_pass_retained,label='Exploratory retained / reference failed',color='#d99538')
ax.barh(yy,s.competing_background,left=s.retained_exploratory,label='Competing-lineage background',color='#afb4bd')
for i,r in s.iterrows():ax.text(r.native_entry+4,i,str(r.retained_exploratory),va='center',fontsize=8)
ax.set_yticks(list(yy),s['sample']);ax.invert_yaxis();ax.set_xlim(0,s.native_entry.max()*1.16);ax.set_xlabel('Cells; right labels = exploratory retained')
ax.set_title('Expanded-layer review: candidate cells, not admitted patients')
ax.legend(loc='lower right',fontsize=8);fig.tight_layout();fig.savefig(O/'figures/final_review_capacity.png',dpi=170);plt.close(fig)
print(s.to_string(index=False))
print(json.dumps(tot))
