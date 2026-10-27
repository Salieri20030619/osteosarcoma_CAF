from pathlib import Path
import json,zipfile,io,pandas as pd,numpy as np,itertools
from statsmodels.stats.multitest import multipletests
r=Path(r'C:\Users\31937\CAF_final_report_20260926');o=Path(r'C:\Users\31937\CAF_report_reaudit_20261010');z=zipfile.ZipFile(r/'output/CAF_statistics_and_evidence_20260927.zip')
for layer in ['supported','expanded']:
 d=pd.read_csv(io.BytesIO(z.read(f'S02/classification_and_programs/{layer}_classified_cells.tsv.gz')),sep='\t',compression='gzip');print('JOINT',layer,int((d.copykat_call.eq('diploid')&d.infercnv_call.eq('no_alteration_called')).sum()))
cp=pd.read_csv(io.BytesIO(z.read('S02/classification_and_programs/class_program_statistics.tsv')),sep='\t');checks=[]
for layer in ['supported','expanded']:
 for cl in cp[cp.pool.eq('SCPCP000017')&cp.layer.eq(layer)].class_id.unique():
  f=f'S02/classification_and_programs/{layer}_SCPCP000017_{cl}_patient_program_effects.tsv';a=pd.read_csv(io.BytesIO(z.read(f)),sep='\t');rows=cp[(cp.layer==layer)&(cp.pool=='SCPCP000017')&(cp.class_id==cl)]
  for _,v in rows.iterrows():
   d=a[v.feature].to_numpy();n=len(d);signs=np.array(list(itertools.product([-1.,1.],repeat=n)));pv=float((abs(signs@d/n)>=abs(d.mean())-1e-12).mean());checks.append({'layer':layer,'class':cl,'feature':v.feature,'N':n,'P_saved':v.P,'P_recomputed':pv,'pass':abs(pv-v.P)<1e-10,'median_recomputed':np.median(d),'median_saved':v.median_delta})
for fam,g in cp.groupby('family'):
 q=multipletests(g.P.fillna(1),method='fdr_bh')[1];ix=g.P.notna().to_numpy();assert np.allclose(q[ix],g.FDR.to_numpy()[ix])
pd.DataFrame(checks).to_csv(o/'independent_main_program_recalculation.tsv',sep='\t',index=False);print('MAIN_PROGRAM_RECALC',len(checks),sum(x['pass'] for x in checks),'BH_ALL_FAMILIES_PASS')
for file in ['RNA_pair_codetection.tsv','reference_mapping_dependence.tsv']:
 d=pd.read_csv(r/'external_singlecell_validation_20260928/tables'/file,sep='\t');print(file,d.columns.tolist());print(d.head(2).to_dict('records'))
