from pathlib import Path
import pandas as pd,numpy as np,itertools,json,zipfile,hashlib,shutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
B=Path('C:/Users/31937/CAF_C3_matrix_strength_20261010');R=B/'results';D=B/'delivery';T=D/'tables';T.mkdir(parents=True,exist_ok=True);(D/'resources').mkdir(exist_ok=True);(D/'code').mkdir(exist_ok=True);(D/'figures').mkdir(exist_ok=True)
for p in R.iterdir():shutil.copy2(p,T/p.name)
shutil.copy2(B/'analyze_functional_strength.py',D/'code/analyze_functional_strength.py')
s=pd.read_csv(R/'functional_patient_statistics.tsv',sep='\t');e=pd.read_csv(R/'patient_effects.tsv',sep='\t');rs=pd.read_csv(R/'residual_processing_statistics.tsv',sep='\t');re=pd.read_csv(R/'residual_patient_effects.tsv',sep='\t');cs=pd.read_csv(R/'cell_functional_scores.tsv.gz',sep='\t');checks=0
for stats,eff,keys in [(s,e,['matrix','min_cells','reference','panel']),(rs,re,['matrix','min_cells','reference'])]:
 for _,r in stats.iterrows():
  q=eff.copy()
  for k in keys:q=q[q[k]==r[k]]
  v=q.effect.to_numpy();assert len(v)==r.N and q.donor_record.nunique()==r.N
  assert np.isclose(v.mean(),r.mean_effect) if len(v) else pd.isna(r.mean_effect)
  assert sum(v>0)==r.positive
  if len(v)>=2 and np.any(v!=0):
   p=sum(abs(np.mean(v*np.asarray(a)))>=abs(v.mean())-1e-12 for a in itertools.product([-1,1],repeat=len(v)))/(2**len(v));assert np.isclose(p,r.P)
  else:assert pd.isna(r.P)
  checks+=1
 for fam,g in stats[stats.P.notna()].groupby('family'):
  o=np.argsort(g.P.to_numpy());v=g.P.to_numpy()[o];adj=np.minimum(1,np.minimum.accumulate((v*len(v)/np.arange(1,len(v)+1))[::-1])[::-1]);assert np.allclose(g.FDR.to_numpy()[o],adj)
# Reconstruct every unadjusted patient effect from exported cell scores, independently.
for _,r in e.iterrows():
 d=cs[(cs.matrix==r.matrix)&(cs.donor_record==r.donor_record)];vals=[]
 for site,z in d.groupby('site'):
  l=z[z.subgroup_id=='F2_C3'];q=z[z.subgroup_id!='F2_C3']
  if r.reference=='other_F2':q=q[q.class_id=='F2']
  elif r.reference!='other_Fibro':q=q[q.subgroup_id==r.reference]
  if min(len(l),len(q))>=r.min_cells:vals.append(l[r.panel].mean()-q[r.panel].mean())
 assert len(vals) and np.isclose(np.mean(vals),r.effect)
# Reconstruct held-out regressions and residual contrasts from scores and model coefficients.
models=pd.read_csv(R/'LODO_models.tsv',sep='\t')
for _,r in re.iterrows():
 d=cs[(cs.matrix==r.matrix)&(cs.donor_record==r.donor_record)].copy();m=models[(models.matrix==r.matrix)&(models.heldout_donor==r.donor_record)].iloc[0];d['residual']=d.processing2-m.intercept-m.beta_COL5_program*d.Published_COL5A1_MSC4-m.beta_library*np.log1p(d.library);vals=[]
 for site,z in d.groupby('site'):
  l=z[z.subgroup_id=='F2_C3'];q=z[z.subgroup_id!='F2_C3']
  if r.reference=='other_F2':q=q[q.class_id=='F2']
  elif r.reference!='other_Fibro':q=q[q.subgroup_id==r.reference]
  if min(len(l),len(q))>=r.min_cells:vals.append(l.residual.mean()-q.residual.mean())
 assert len(vals) and np.isclose(np.mean(vals),r.effect)
A=Path('C:/Users/31937/CAF_all_clusters_20261006/analysis')
pg=pd.read_csv(A/'tables/all_source_program_statistics.tsv',sep='\t');pg[(pg.dataset=='SCPCP000017')&pg.subgroup_id.isin(['F2_C3','F2_C1','F2_C0'])].to_csv(T/'historical_C3_C1_C0_program_statistics.tsv',sep='\t',index=False)
g=pd.read_csv(A/'tables/nominal_upregulated_all_sources_with_sensitivity.tsv',sep='\t');g[g.subgroup_id=='F2_C3'].to_csv(T/'historical_C3_nominal_genes.tsv',sep='\t',index=False)
for name in ['PROTOCOL.json','PROGRAMS.json']:shutil.copy2(A/'resources'/name,D/'resources'/('historical_'+name))
z=next(Path('C:/Users/31937/Downloads').glob('CAF*GSEA*.zip'))
with zipfile.ZipFile(z) as a:
 for tail in ['GSEA_all_33_subgroups_477_records.tsv','GSEA_full_original_families.tsv.gz','GSEA_PROTOCOL_phase62.json']:
  name=next(n for n in a.namelist() if n.endswith(tail));dest=(D/'resources' if tail.endswith('.json') else T)/tail;dest.write_bytes(a.read(name))
 gd=pd.read_csv(T/'GSEA_all_33_subgroups_477_records.tsv',sep='\t');gd[gd.contrast.isin(['F2_C3','F2_C1'])].to_csv(T/'historical_C3_C1_GSEA.tsv',sep='\t',index=False)
# Verify preserved model BH using all original hypotheses, not the selected 477 rows.
gf=pd.read_csv(T/'GSEA_full_original_families.tsv.gz',sep='\t');familykeys=['context','library'];gn=0
for key,q in gf.groupby(familykeys):
 if not ('F2_C3' in str(key) or 'F2_C1' in str(key)):continue
 v=q.GSEA_model_P.fillna(1).to_numpy();o=np.argsort(v);adj=np.minimum(1,np.minimum.accumulate((v[o]*len(v)/np.arange(1,len(v)+1))[::-1])[::-1]);assert np.allclose(adj,q.model_BH_FDR.to_numpy()[o]);gn+=1
# Patient effects figure, dots are patients, not cells.
fig,axs=plt.subplots(1,2,figsize=(11.6,5.5),gridspec_kw={'width_ratios':[1.5,1]},constrained_layout=True)
items=[('Published_COL5A1_MSC4','Published COL5A1 program'),('C3_assembly5','Collagen / processing core'),('processing2','PCOLCE + LOXL2'),('residual','Processing after adjustment')]
for y,(panel,label) in enumerate(items):
 q=(re if panel=='residual' else e);q=q[(q.matrix=='corrected')&(q.min_cells==1)&(q.reference=='other_Fibro')]
 if panel!='residual':q=q[q.panel==panel]
 vals=q.effect.to_numpy();axs[0].scatter(vals,np.repeat(y,len(vals))+np.linspace(-.085,.085,len(vals)),s=35,color='#287fa1',alpha=.85);axs[0].scatter([vals.mean()],[y],marker='D',s=66,color='#b54235',zorder=4)
axs[0].set_yticks(range(4),[x[1] for x in items]);axs[0].invert_yaxis();axs[0].axvline(0,color='#888',lw=1);axs[0].set_xlabel('C3 minus other Fibro: mean log1p(CP10k) score');axs[0].set_title('A. Patient-matched effects (N = 7)',loc='left',fontweight='bold');axs[0].grid(axis='x',alpha=.16)
q=e[(e.matrix=='corrected')&(e.min_cells==1)&(e.reference=='F2_C1')]
for y,(panel,label) in enumerate(items[1:3]):
 v=q[q.panel==panel].effect.to_numpy();axs[1].scatter(v,np.repeat(y,len(v))+np.linspace(-.06,.06,len(v)),s=38,color='#287fa1');axs[1].scatter([v.mean()],[y],marker='D',s=66,color='#b54235')
axs[1].set_yticks([0,1],['Collagen core','Processing']);axs[1].set_ylim(1.5,-.5);axs[1].axvline(0,color='#888',lw=1);axs[1].set_xlabel('C3 minus C1: score difference');axs[1].set_title('B. Direct C1 comparison (N = 4)',loc='left',fontweight='bold');axs[1].grid(axis='x',alpha=.16)
fig.suptitle('C3 collagen-processing strength: patient-level comparisons',fontsize=14,fontweight='bold');fig.savefig(D/'figures/patient_functional_effects.png',dpi=180);fig.savefig(D/'figures/patient_functional_effects.pdf');plt.close(fig)
qa={'status':'passed','remote_analysis_exit_code':0,'statistical_records_recomputed':checks,'functional_patient_effects_reconstructed':len(e),'residual_patient_effects_reconstructed':len(re),'all_new_BH_families_verified':True,'original_C3_C1_GSEA_families_verified':gn,'cells_per_matrix':int(len(cs)/2),'C3_patients':int(cs[(cs.matrix=='corrected')&(cs.subgroup_id=='F2_C3')].donor_record.nunique()),'source_GSEA_zip':str(z),'source_GSEA_zip_sha256':hashlib.sha256(z.read_bytes()).hexdigest(),'new_vs_old_FDR_note':'phase76 newly defined three functional panels and seven comparisons; historical families remain unchanged, not a replacement for phase75 FDR'}
(D/'QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(qa,ensure_ascii=False))
