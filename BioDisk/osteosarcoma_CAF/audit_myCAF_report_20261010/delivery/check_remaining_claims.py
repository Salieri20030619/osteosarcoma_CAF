from pathlib import Path
import pandas as pd,numpy as np,json,re,zipfile,io
R=Path(r'C:\Users\31937\CAF_final_report_20260926');O=Path(r'C:\Users\31937\CAF_report_reaudit_20261010');T={(x['page'],x['table']):x['rows'] for x in json.loads((O/'pdf_tables.json').read_text(encoding='utf-8'))};D=json.loads((R/'inputs/report_data.json').read_text(encoding='utf-8'));checks=[]
def add(page,claim,ok,evidence):checks.append({'page':page,'claim':claim,'status':'PASS' if ok else 'CHECK','evidence':str(evidence)})
l=pd.read_csv(R/'lung_resident_review_20260928/lung_resident_source_summary.tsv',sep='\t')
for i,row in enumerate(T[11,1][1:],1):
 for col,layer in [(1,'supported'),(2,'expanded')]:
  a=l[(l.dataset==row[0])&(l.scope=='F2')&(l.layer==layer)].iloc[0];s=row[col];nums=re.findall(r'\d+(?:\.\d+)?',s)
  ok=('NA' in s and a.evaluable_patients==0) if 'NA' in s else abs(float(nums[0])-100*a.patient_mean_fraction)<=.050001 and int(nums[1])==a.evaluable_patients and int(nums[2])==a.candidate_cells
  add(11,f'{row[0]} {layer} F2比例',ok,a.to_dict())
b=pd.read_csv(R/'external_singlecell_validation_20260928/tables/RNA_pair_codetection.tsv',sep='\t')
for i,(pair,group) in enumerate([('PDGFRB__MCAM','Fibroblast'),('PDGFRB__MCAM','Mural'),('GPC3__THY1','Fibroblast')],1):
 for col,num in [(1,7),(2,8)]:
  sub=b[(b.tissue=='primary')&(b.patient==f'GSE237070_patient_{num}')&(b.pair==pair)&(b.group==group)&(b.grouping=='lineage_surface_markers_held_out')]
  a=sub.iloc[0];nums=re.findall(r'\d+(?:\.\d+)?',T[13,2][i][col].replace(',',''));add(13,f'{pair} {group}患者{num}',int(nums[0])==a.both_detected and int(nums[1])==a.cells and abs(float(nums[2])-100*a.fraction_both)<=.005001,a.to_dict())
mapd=pd.read_csv(R/'external_singlecell_validation_20260928/tables/reference_mapping_dependence.tsv',sep='\t');pr=mapd[mapd['sample'].isin(['GSM7594318','GSM7594321'])];add(14,'两参考F2人数39和3050',pr.supported_F2_like.sum()==39 and pr.expanded_F2_like.sum()==3050,pr.to_dict('records'))
for i,layer in enumerate(['expanded','supported','supported_myCAF','expanded_myCAF'],1):
 a=pd.DataFrame(D['sensitivity']);a=a[a.layer==layer];nums=[float(z) for z in re.findall(r'\d+\.\d+',T[15,1][i][1])];add(15,layer+' ARI范围',abs(nums[0]-a.ARI.min())<=.0005001 and abs(nums[1]-a.ARI.max())<=.0005001,{'min':a.ARI.min(),'max':a.ARI.max(),'baseline_clusters':a.baseline_clusters.unique().tolist()})
g=pd.read_csv(R/'focus_audit_20260928/myCAF_customer_gene_statistics.tsv',sep='\t');genes=D['panels']['Osteogenic_customer']+D['panels']['Mineral_customer'];x=g[(g.study=='SCPCP000017')&g.symbol.isin(genes)];add(8,'31基因全部内部比较无有效P<.05',x.P.lt(.05).sum()==0,{'rows':len(x),'tested':int(x.P.notna().sum()),'NA':int(x.P.isna().sum()),'evaluation':x.evaluation_status.value_counts().to_dict()})
sp=pd.read_csv(R/'lung_resident_review_20260928/Gong_patient_statistics.tsv',sep='\t');print('GONG',sp.columns.tolist());add(11,'108项Gong比较无P<.05',len(sp)==108 and sp.P.lt(.05).sum()==0,{'rows':len(sp),'valid':int(sp.P.notna().sum()),'min_P':sp.P.min()})
add(22,'总样本和细胞数',sum(v['sample_records'] for v in D['datasets'])==107 and sum(v['full_called_cells'] for v in D['datasets'])==876648,D['scope'])
z=zipfile.ZipFile(R/'output/CAF_statistics_and_evidence_20260927.zip');d=pd.read_csv(io.BytesIO(z.read('S06/patient_effect_review/new_doublet_calls_by_class.tsv')),sep='\t');print('DOUBLET',d.to_dict('records'))
# The recorded primary six-gene result must match historical P and median before new sensitivity families are interpreted.
a=pd.read_csv(O/'six_gene_recheck/exact_six_gene_patient_statistics.tsv',sep='\t');x=a[(a.matrix=='raw')&(a.scope=='frozen_all')&(a.target=='F2')].iloc[0];old=next(v for v in D['main_program_statistics'] if v['layer']=='expanded' and v['class_id']=='F2' and v['feature']=='program_myCAF_customer');add(5,'原始6基因从完整矩阵重算复现',x.N==old['n_independent_donors'] and abs(x.P-old['P'])<1e-9 and abs(x.median_effect-old['median_delta'])<1e-6,{'old':old,'new':x.to_dict()})
pd.DataFrame(checks).to_csv(O/'supplementary_claim_checks.tsv',sep='\t',index=False);print('SUPPLEMENT_CHECKS',len(checks),pd.DataFrame(checks).status.value_counts().to_dict())
