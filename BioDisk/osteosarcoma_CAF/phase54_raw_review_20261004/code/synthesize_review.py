from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
import json,numpy as np,pandas as pd
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase54_raw_review_20261004';S=P/'phase53_de_novo_reannotation_20261002';V=P/'review_intake_20261004/received_review'
def read(p,**kw):return pd.read_csv(p,sep='\t',low_memory=False,**kw)
def write(df,n):df.to_csv(O/'tables'/n,sep='\t',index=False)
for stage,marker in [('raw_review','RAW_REVIEW_COMPLETE'),('cnv_review','CNV_REVIEW_COMPLETE'),('doublet_metadata','DOUBLET_METADATA_COMPLETE')]:assert (O/(stage+'.exit_code')).read_text().strip()=='0' and (O/marker).is_file()
d=pd.concat([read(f) for f in sorted((O/'cells').glob('*_target_cells.tsv.gz'))],ignore_index=True);assert d.cell_id.is_unique
bm=read(O/'tables/target_cell_doublet_evidence.tsv.gz').set_index('cell_id');assert set(d.cell_id)==set(bm.index);bm=bm.loc[d.cell_id]
methodcols=['input_scDblFinder_class','phase31_scDblFinder_class','input_predicted_doublet']
valid=[];positive=[]
for col in methodcols:
 vals=bm[col].astype(str).str.lower() if col in bm else pd.Series('missing',index=bm.index)
 valid.append(vals.isin(['singlet','doublet','true','false','1','0']).to_numpy());positive.append(vals.isin(['doublet','true','1']).to_numpy())
d['doublet_assessed_any']=np.any(valid,axis=0);d['doublet_positive_any']=np.any(positive,axis=0)
d['doublet_evidence_group']=np.select([d.doublet_positive_any,d.doublet_assessed_any],['model_positive','assessed_no_positive'],default='unassessed')
structure=[]
for source in d.source_id.unique():
 x=read(S/'blind_review'/source/'cell_structure.tsv.gz',usecols=['cell_id','r05','r10','r20']);structure.append(x)
ss=pd.concat(structure);d=d.merge(ss[['cell_id','r05','r20']],on='cell_id',validate='one_to_one')
primary={'D01_R10_K002':('PF4','C1QA'),'D01_R10_K034':('IFITM5','C1QA'),'D04_R10_K007':('CD163','DCN'),'D06_R10_K022':('RUNX2','CD163'),'D07_R10_K027':('PECAM1','SFTPB'),'D09_R10_K010':('C1QA','COL3A1')}
fm=read(O/'tables/frozen_pair_feature_mapping.tsv');pm=read(O/'tables/fixed_panel_feature_mapping.tsv')
for cl,genes in primary.items():
 for gene in genes:
  frozen=fm.loc[(fm.cluster_id==cl)&(fm.gene==gene),'feature_id'].iloc[0];fixed=pm.loc[(pm.source_id==cl[:3])&(pm.gene==gene),'feature_id'].iloc[0];assert frozen==fixed
mixsub=[]
for cl,(ga,gb) in primary.items():
 x=d[d.r10.eq(cl)].copy()
 x['raw_both']=(x['raw_'+ga]>0)&(x['raw_'+gb]>0);x['corrected_both_ge1']=(x['corrected_'+ga]>=1)&(x['corrected_'+gb]>=1)
 for dim in ['patient_id','r20','doublet_evidence_group']:
  for name,g in x.groupby(dim):mixsub.append({'cluster_id':cl,'gene_A':ga,'gene_B':gb,'stratum_type':dim,'stratum':name,'cells':len(g),'raw_both':int(g.raw_both.sum()),'raw_fraction':g.raw_both.mean(),'corrected_both_ge1':int(g.corrected_both_ge1.sum()),'corrected_fraction_ge1':g.corrected_both_ge1.mean(),'doublet_assessed':int(g.doublet_assessed_any.sum()),'doublet_positive':int(g.doublet_positive_any.sum()),'contamination_median':g.decontx_contamination.median(),'raw_library_median':g.raw_library.median()})
write(pd.DataFrame(mixsub),'priority_mixed_stratified_review.tsv')
summary=[];patientrows=[];subrows=[]
for cl,g in d.groupby('r10'):
 z={'cluster_id':cl,'source_id':g.source_id.iloc[0],'cells':len(g),'patients':g.patient_id.nunique(),'samples':g.sample_id.nunique(),'max_patient_fraction':g.patient_id.value_counts(normalize=True).max(),'largest_patient':g.patient_id.value_counts().index[0],'contamination_median':g.decontx_contamination.median(),'raw_library_median':g.raw_library.median(),'corrected_library_median':g.corrected_library.median(),'doublet_assessed_any':int(g.doublet_assessed_any.sum()),'doublet_positive_any':int(g.doublet_positive_any.sum()),'doublet_positive_fraction_assessed':g.doublet_positive_any.sum()/max(g.doublet_assessed_any.sum(),1)}
 for gene in ['DCN','LUM','COL3A1','PDGFRA','RUNX2','IFITM5','RGS5','CSPG4']:
  z[gene+'_raw_detection']=g['raw_'+gene].gt(0).mean();z[gene+'_corrected_ge1']=g['corrected_'+gene].ge(1).mean();z[gene+'_count_retention']=g['corrected_'+gene].sum()/g['raw_'+gene].sum() if g['raw_'+gene].sum()>0 else np.nan
 summary.append(z)
 cols=[c for c in g if c.endswith('_mean_logCP10k')]
 for patient,v in g.groupby('patient_id'):
  q={'cluster_id':cl,'patient_id':patient,'cells':len(v),'samples':v.sample_id.nunique(),'contamination_median':v.decontx_contamination.median()}
  q.update({c:v[c].median() for c in cols});patientrows.append(q)
 for sub,v in g.groupby('r20'):
  q={'cluster_id':cl,'r20':sub,'cells':len(v),'fraction':len(v)/len(g),'patients':v.patient_id.nunique(),'doublet_positive_fraction':v.doublet_positive_any.mean(),'contamination_median':v.decontx_contamination.median()}
  q.update({c:v[c].median() for c in cols});subrows.append(q)
write(pd.DataFrame(summary),'target_cluster_raw_review_summary.tsv');write(pd.DataFrame(patientrows),'target_patient_panel_summary.tsv');write(pd.DataFrame(subrows),'target_r20_panel_summary.tsv')
write(d[['cell_id','original_cell_id','source_id','sample_id','patient_id','r10','r05','r20','raw_library','corrected_library','decontx_contamination','doublet_assessed_any','doublet_positive_any','doublet_evidence_group']],'target_cell_review_overlay.tsv.gz')
# Pooled same-sample osteogenic comparator avoids choosing the best matching cluster.
pr=read(O/'tables/CNV_chromosome_profiles_recomputed.tsv.gz');det=read(V/'blind_submission/cluster_review_details.tsv');lm=det.set_index('cluster_id').lineage_code;pr['lineage_code']=pr.cluster_id.map(lm);pooled=[]
for (sid,chrom),g in pr[pr.lineage_code.eq('Osteo')].groupby(['sample_id','chromosome']):
 pooled.append({'sample_id':sid,'chromosome':chrom,'pooled_osteo_cells':g.cells.sum(),'pooled_osteo_signal':np.average(g.ref_centered_signal,weights=g.cells)})
z=pr.merge(pd.DataFrame(pooled),on=['sample_id','chromosome'],how='left',validate='many_to_one');compar=[]
for (sid,cl),g in z[z.cluster_id.isin(set(det.loc[det.lineage_code.isin(['Fibro','Mesench']),'cluster_id']))].groupby(['sample_id','cluster_id']):
 eligible=bool(g.adequate_reference.iloc[0] and g.cells.iloc[0]>=20 and g.pooled_osteo_cells.iloc[0]>=20);finite=g.ref_centered_signal.notna()&g.pooled_osteo_signal.notna();a=g.loc[finite,'ref_centered_signal'];b=g.loc[finite,'pooled_osteo_signal'];r=float(np.corrcoef(a,b)[0,1]) if eligible and len(a)>2 and a.std()>0 and b.std()>0 else np.nan
 compar.append({'source_id':g.source_id.iloc[0],'sample_id':sid,'patient_id':g.patient_id.iloc[0],'cluster_id':cl,'target_copykat_cells':g.cells.iloc[0],'pooled_osteo_cells':g.pooled_osteo_cells.iloc[0],'eligible':eligible,'chromosomes':int(finite.sum()),'pooled_osteo_chromosome_r':r,'no_comparator_selection_by_correlation':True})
write(pd.DataFrame(compar),'CNV_pooled_osteogenic_similarity.tsv')
# Plots are descriptive, not inference or diagnosis.
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
F=O/'figures';F.mkdir(exist_ok=True);plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42})
ps=read(O/'tables/cluster_pair_codetection.tsv');rows=[]
for cl,(a,b) in primary.items():rows.append(ps[(ps.cluster_id==cl)&(ps.gene_A==a)&(ps.gene_B==b)].iloc[0])
q=pd.DataFrame(rows);fig,ax=plt.subplots(figsize=(10,5));y=np.arange(len(q));ax.barh(y-.18,q.raw_both_fraction*100,.34,label='Raw joint detection > 0',color='#286c8e');ax.barh(y+.18,q.corrected_both_ge1_fraction*100,.34,label='DecontX: both >= 1 count',color='#d28b37');ax.set_yticks(y,[f'{r.cluster_id}\n{r.gene_A} + {r.gene_B} (n={r.cells})' for r in q.itertuples()]);ax.set_xlim(0,105);ax.set_xlabel('Percentage of barcodes within cluster');ax.invert_yaxis();ax.legend(loc='lower right');ax.set_title('Mixed transcription signals: raw counts and correction sensitivity');fig.tight_layout();fig.savefig(F/'mixed_codetection.png',dpi=180);fig.savefig(F/'mixed_codetection.pdf');plt.close(fig)
x=read(O/'tables/CNV_patient_recomputed.tsv');x=x[x.cluster_id.eq('D01_R10_K017')].dropna(subset=['fold']).sort_values('patient_id');fig,ax=plt.subplots(figsize=(8,4));ax.bar(x.patient_id,x.fold,color='#625092');ax.axhline(1,color='gray',lw=1);ax.axhline(2,color='gray',ls='--',lw=1);ax.set_ylabel('Native CNV amplitude / sample reference');ax.set_title('D01_K017: patient-level RNA CNV amplitude');fig.text(.01,.01,'Reference is expression-defined; ratios and guide lines do not diagnose malignancy.',fontsize=8);fig.tight_layout(rect=[0,.05,1,1]);fig.savefig(F/'D01_K017_patient_CNV.png',dpi=180);fig.savefig(F/'D01_K017_patient_CNV.pdf');plt.close(fig)
(O/'SYNTHESIS_COMPLETE').write_text('SUCCESS\n');print(pd.DataFrame(summary).loc[pd.DataFrame(summary).cluster_id.isin(list(primary)+['D01_R10_K017','D02_R10_K012','D05_R10_K018','D07_R10_K024','D07_R10_K041'])].to_json(orient='records'))