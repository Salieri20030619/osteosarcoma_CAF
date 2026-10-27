from pathlib import Path
import json,re,io,zipfile,gzip
import numpy as np,pandas as pd
R=Path(r'C:\Users\31937\CAF_final_report_20260926');O=Path(r'C:\Users\31937\CAF_report_reaudit_20261010');D=json.loads((R/'inputs/report_data.json').read_text(encoding='utf-8'));T={(x['page'],x['table']):x['rows'] for x in json.loads((O/'pdf_tables.json').read_text(encoding='utf-8'))};rows=[]
def ck(p,t,r,c,value,src):
 actual=T[p,t][r][c].replace('\n',' ').strip();raw=actual.replace(',','').replace('−','-').replace('%','').strip()
 if pd.isna(value):ok=raw in ['NA','nan',''];delta=None
 else:
  try:
   num=float(raw);dec=len(raw.split('.')[-1]) if '.' in raw else None;tol=.50001*10**(-dec) if dec is not None else 1e-9;delta=num-float(value);ok=abs(delta)<=tol
  except ValueError:ok=actual.replace(' ','')==str(value).replace(' ','');delta=None
 rows.append(dict(page=p,table=t,row=r,column=c,printed=actual,source_value=value,status='PASS' if ok else 'MISMATCH',difference=delta,source=src))
def one(df,**kw):
 for k,v in kw.items():df=df[df[k].eq(v)]
 assert len(df)==1,(kw,len(df));return df.iloc[0]
def rowstats(p,t,r,x,fields):
 for col,key in fields.items():ck(p,t,r,col,x[key],str(x.get('family','historical_table'))+'|'+str(x.get('feature',x.get('symbol',''))))
Z=zipfile.ZipFile(R/'output/CAF_statistics_and_evidence_20260927.zip')
def zread(n):return pd.read_csv(io.BytesIO(Z.read(n)),sep='\t',compression='gzip' if n.endswith('.gz') else None,low_memory=False)
cp=zread('S02/classification_and_programs/class_program_statistics.tsv');cs=zread('S02/classification_and_programs/class_surface_gene_statistics.tsv')
# pp4 and13 frozen counts directly from archived cell records.
for layer,col in [('supported',2),('expanded',3)]:
 cells=zread(f'S02/classification_and_programs/{layer}_classified_cells.tsv.gz')
 for i,cl in enumerate(['F1','F2','F3','F4','F5'],1):ck(4,1,i,col,int(cells.class_id.eq(cl).sum()),'archived classified_cells')
 print('CNV_COUNTS',layer,len(cells),cells.copykat_call.value_counts().to_dict(),cells.infercnv_call.value_counts().to_dict())
 for c,v in [(1,int(cells.copykat_call.eq('diploid').sum())),(2,int(cells.infercnv_call.eq('no_alteration_called').sum())),(3,int((cells.copykat_call.eq('diploid') & cells.infercnv_call.eq('no_alteration_called')).sum()))]:
  if True:ck(13,3,1 if layer=='supported' else 2,c,v,'archived classified_cells')
for i,(layer,feature) in enumerate([('supported','program_myCAF_customer'),('expanded','program_myCAF_customer'),('expanded','program_mCAF_matrix')],1):
 x=one(cp,layer=layer,pool='SCPCP000017',class_id='F2',feature=feature);rowstats(5,1,i,x,{1:'n_independent_donors',2:'median_delta',4:'P',5:'FDR'});ck(5,1,i,3,f'{int(x.positive_donors)}/{int(x.n_independent_donors)}','class_program_statistics')
for i,(layer,sym) in enumerate([('supported','PDGFRB'),('expanded','PDGFRB'),('supported','MCAM'),('expanded','MCAM')],1):
 x=one(cs,layer=layer,pool='SCPCP000017',class_id='F2',symbol=sym);rowstats(5,2,i,x,{1:'n_independent_donors',2:'median_delta',4:'P',5:'FDR'});ck(5,2,i,3,f'{int(x.positive_donors)}/{int(x.n_independent_donors)}','class_surface_gene_statistics')
for i,row in enumerate(T[6,1][1:],1):
 layer='expanded' if row[0].startswith('E') else 'supported';cluster=int(row[0][-1]);x=one(pd.DataFrame(D[layer+'_myCAF_counts']),cluster=cluster);ck(6,1,i,1,x.cells,'internal cluster counts')
g=pd.read_csv(R/'focus_audit_20260928/myCAF_customer_gene_statistics.tsv',sep='\t');sel=pd.DataFrame(D['selected_comparisons'])
for table,syms in [(1,['ACTA2','TAGLN','MYH11','MYL9','TPM2','CNN1']),(2,['MCAM','PDGFRB','JAG1','THY1','FAP'])]:
 for i,sym in enumerate(syms,1):
  pool=g if table==1 else sel;x=one(pool,layer='expanded',study='SCPCP000017',analysis='myCAF_internal',target='subgroup_1',contrast='subgroup_vs_other_myCAF',symbol=sym)
  rowstats(7,table,i,x,{1:'n_first',2:'median_effect',4:'P',5:'FDR'});ck(7,table,i,3,f'{int(x.positive_donors)}/{int(x.n_first)}','internal gene statistics')
pr=pd.read_csv(R/'focus_audit_20260928/myCAF_main_program_statistics.tsv',sep='\t')
for i,row in enumerate(T[8,1][1:],1):
 layer='expanded' if row[0].startswith('E') else 'supported';target='subgroup_'+row[0][-1]
 for feature,ec,pc in [('program_Osteogenic_customer',2,3),('program_Mineral_customer',4,5)]:
  x=one(pr,layer=layer,study='SCPCP000017',target=target,feature=feature);ck(8,1,i,1,x.n_first,'internal program statistics');ck(8,1,i,ec,x.median_effect,'internal program statistics');ck(8,1,i,pc,x.P,'internal program statistics')
for i,(layer,target) in enumerate([('supported','all_CAF_candidates'),('expanded','all_CAF_candidates'),('supported','F3')],1):
 x=one(pd.DataFrame(D['strict_site_all_surface']),layer=layer,target=target,symbol='GPC3',paired=False);rowstats(9,1,i,x,{2:'median_effect',3:'P',4:'FDR'});ck(9,1,i,1,f'{int(x.n_bone)}/{int(x.n_lung)}','strict anatomical comparison')
l=pd.read_csv(R/'lung_resident_review_20260928/lung_resident_source_summary.tsv',sep='\t')
for i,row in enumerate(T[10,1][1:],1):
 x=one(l,dataset=row[0],layer='supported' if row[1]=='支持' else 'expanded',scope='all_candidates');ck(10,1,i,2,f'{int(x.evaluable_patients)}/{int(x.source_patients)}','lung resident summary')
 for col,key in [(3,'candidate_cells'),(4,'AdvF_cells'),(5,'AlvF_cells')]:ck(10,1,i,col,x[key],'lung resident summary')
 ck(10,1,i,6,x.patient_mean_fraction*100,'lung resident summary')
SZ=zipfile.ZipFile(R/'external_validation_results/delivery/CAF_external_validation_complete_20260927.zip')
def sr(n):return pd.read_csv(io.BytesIO(SZ.read('tables/'+n)),sep='\t')
sp=sr('independent_R_site_statistics_QA.tsv');cc=sr('marker_component_patient_level_tests.tsv')
print('SPATIAL_CONTRASTS',sp.contrast.unique())
for i,sym in enumerate(['GPC3','GPRC5A','PDGFRB','MCAM'],1):
 x=one(sp,feature=sym,contrast='independent_lung_minus_primary');rowstats(12,1,i,x,{1:'mean_effect',2:'P',3:'FDR'})
 x=one(sp,feature=sym,contrast='paired_lung_minus_primary');txt=T[12,1][i][4];nums=re.findall(r'\d+(?:\.\d+)?',txt);rows.append(dict(page=12,table=1,row=i,column=4,printed=txt,source_value=f'{x.positive_donors}/4; P={x.P}',status='PASS' if int(nums[0])==x.positive_donors and np.isclose(float(nums[-1]),x.P) else 'MISMATCH',source='paired spatial statistics'))
for i,(sym,group) in enumerate([('GPC3','fibroblast'),('MCAM','endothelial'),('MCAM','mural'),('PDGFRB','mural'),('PDGFRB','fibroblast')],1):
 x=one(cc,marker=sym,component_group=group);rowstats(13,1,i,x,{1:'median_rho',3:'P',4:'FDR'});ck(13,1,i,2,f'{int(x.positive_patients)}/{int(x.n_patients)}','marker_component_patient_level_tests')
ext=pd.read_csv(R/'external_singlecell_validation_20260928/tables/patient_level_marker_program_statistics.tsv',sep='\t')
for i,(gene,layer) in enumerate([('PDGFRB','expanded'),('PDGFRB','supported'),('MCAM','expanded'),('MCAM','supported')],1):
 x=one(ext,grouping=layer+'_reference_like',tissue='primary',left='F2',right='rest_of_F2',endpoint=gene);rowstats(14,1,i,x,{2:'mean_effect',4:'P',5:'FDR'});ck(14,1,i,3,f'{int(x.positive_patients)}/{int(x.N_patients)}','external singlecell patient statistics')
for i,(layer,cl,ft) in enumerate([('expanded','F4','Osteogenic_customer'),('expanded','F4','Mineral_customer'),('expanded','F2','Osteogenic_customer'),('expanded','F2','Mineral_customer'),('supported','F3','mCAF_matrix')],1):
 x=one(cp,layer=layer,pool='SCPCP000017',class_id=cl,feature='program_'+ft);rowstats(19,1,i,x,{1:'n_independent_donors',2:'median_delta',4:'P',5:'FDR'});ck(19,1,i,3,f'{int(x.positive_donors)}/{int(x.n_independent_donors)}','class_program_statistics')
for i,target in enumerate(['all_candidates','F4','F5'],1):
 x=one(pd.DataFrame(D['clinical_proportion_nominal']),layer='expanded',study='SCPCP000023',target=target);ck(19,2,i,1,f'{int(x.n_first)}/{int(x.n_second)}','treatment proportion stats');ck(19,2,i,2,x.mean_first*100,'treatment proportion stats');ck(19,2,i,3,x.mean_second*100,'treatment proportion stats');ck(19,2,i,4,x.P,'treatment proportion stats');ck(19,2,i,5,x.FDR,'treatment proportion stats')
for i,x in enumerate(D['datasets'],1):
 j=next(j for j,row in enumerate(T[22,1]) if row[0]==x['dataset'])
 for col,key in [(1,'sample_records'),(2,'donor_records'),(3,'full_called_cells')]:ck(22,1,j,col,x[key],'source metadata')
cnt=pd.read_csv(R/'external_singlecell_validation_20260928/tables/independent_count_QA.tsv',sep='\t')
for i,row in enumerate(T[22,2][1:],1):ck(22,2,i,3,one(cnt,sample=row[0]).cells,'independent external count QA')
out=pd.DataFrame(rows);out.to_csv(O/'numeric_reconciliation.tsv',sep='\t',index=False);print('NUMERIC_CHECKS',len(out));print(out.status.value_counts().to_dict());print(out[out.status!='PASS'].to_string(index=False))




