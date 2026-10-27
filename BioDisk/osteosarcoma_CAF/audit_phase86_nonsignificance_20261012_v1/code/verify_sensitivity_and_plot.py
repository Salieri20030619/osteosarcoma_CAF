
from pathlib import Path
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='2'
import pandas as pd,numpy as np,json
from scipy.stats import permutation_test
from statsmodels.stats.multitest import multipletests
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'audit_phase86_nonsignificance_20261012_v1';B=P/'phase86_functional_bone_lung_20261012'
f=pd.read_csv(B/'tables/all_functional_bone_lung_statistics.tsv',sep='\t');q=pd.read_csv(O/'tables/19_all800_equal_sample_weight_sensitivity.tsv',sep='\t');a=pd.read_csv(O/'tables/17_equal_sample_patient_endpoints.tsv',sep='\t');original=pd.read_csv(B/'tables/patient_site_programs.tsv',sep='\t')
def ids(v):return [] if pd.isna(v) else v.split(';')
def exact(x,y):return permutation_test((np.array(x),np.array(y)),lambda a,b:abs(b.mean()-a.mean()),n_resamples=np.inf,alternative='greater',vectorized=False).pvalue
errs=[]
for r in q[q.nominal_sensitivity].itertuples():
 sel=(f.scope==r.scope)&(f.dataset==r.dataset)&(f.matrix==r.matrix)&(f.program==r.program)&(f.minimum_cells==r.minimum_cells)&(f.design==r.design)&(f.endpoint==r.endpoint);v=f[sel].iloc[0]
 t=a[(a.scope==r.scope)&(a.dataset==r.dataset)&(a.matrix==r.matrix)&(a.program==r.program)];x=t[t.site=='Bone'].set_index('donor_record').loc[ids(v.bone_patients),r.endpoint].to_numpy();y=t[t.site=='Lung'].set_index('donor_record').loc[ids(v.lung_patients),r.endpoint].to_numpy();pv=exact(x,y);errs.append(abs(pv-r.equal_sample_P));assert abs(pv-r.equal_sample_P)<1e-10
for fam,t in q.groupby('family'):
 t=t[t.equal_sample_P.notna()]
 if len(t):assert np.max(abs(multipletests(t.equal_sample_P,method='fdr_bh')[1]-t.equal_sample_FDR))<1e-10
# Perturb only the identified donor's within-patient sample weight to diagnose its leverage.
rows=[]
for matrix in ['raw','corrected']:
 r=f[(f.scope=='Fibro_reviewed')&(f.dataset=='SCPCP000017')&(f.matrix==matrix)&(f.program=='contraction')&(f.minimum_cells==1)&(f.design=='independent')&(f.endpoint=='mean_score')].iloc[0]
 s=(a.scope=='Fibro_reviewed')&(a.dataset=='SCPCP000017')&(a.matrix==matrix)&(a.program=='contraction');g=a[s]
 x=g[g.site=='Bone'].set_index('donor_record').loc[ids(r.bone_patients),'mean_score'];y=g[g.site=='Lung'].set_index('donor_record').loc[ids(r.lung_patients),'mean_score'];donor='SCPCP000017|CCG1106_007'
 oo=original[(original.scope=='Fibro_reviewed')&(original.dataset=='SCPCP000017')&(original.matrix==matrix)&(original.program=='contraction')&(original.donor_record==donor)&(original.site=='Lung')].iloc[0].mean_score
 for name,yy in [('all_samples_equal',y),('only_CCG1106_007_restored_to_original_cell_weight',y.copy()),('CCG1106_007_omitted_as_patient_influence_check',y.drop(donor))]:
  if name.startswith('only_'):yy.loc[donor]=oo
  rows.append(dict(matrix=matrix,program='contraction',diagnostic=name,bone_N=len(x),lung_N=len(yy),effect=yy.mean()-x.mean(),P_diagnostic=exact(x.to_numpy(),yy.to_numpy()),note='Mechanism-of-sensitivity check only; not a new candidate-screening family or replacement primary test'))
l=pd.DataFrame(rows);l.to_csv(O/'tables/24_contraction_weight_leverage.tsv',sep='\t',index=False)
print(l.to_string(index=False))
qa={'status':'PASS','new_nominal_tests_independently_checked_with_scipy':len(errs),'max_P_error':max(errs),'all_equal_weight_FDR_families_checked':True,'original_analysis_unchanged':True}
(O/'SENSITIVITY_QA.json').write_text(json.dumps(qa,indent=2));print(json.dumps(qa))
# Publication-style static audit figure, with one row per independent biological donor.
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
(O/'figures').mkdir(exist_ok=True)
r=f[(f.scope=='Fibro_reviewed')&(f.dataset=='SCPCP000017')&(f.matrix=='raw')&(f.program=='contraction')&(f.minimum_cells==1)&(f.design=='independent')&(f.endpoint=='mean_score')].iloc[0]
fig,ax=plt.subplots(1,2,figsize=(11,7),gridspec_kw={'width_ratios':[1.4,1]})
g=original[(original.scope=='Fibro_reviewed')&(original.dataset=='SCPCP000017')&(original.matrix=='raw')&(original.program=='contraction')];t=a[(a.scope=='Fibro_reviewed')&(a.dataset=='SCPCP000017')&(a.matrix=='raw')&(a.program=='contraction')];roster=[('Bone',d) for d in ids(r.bone_patients)]+[('Lung',d) for d in ids(r.lung_patients)]
for i,(site,d) in enumerate(roster):
 v=g[(g.site==site)&(g.donor_record==d)].iloc[0].mean_score;w=t[(t.site==site)&(t.donor_record==d)].iloc[0].mean_score;col='#2878B5' if site=='Bone' else '#C8553D'
 ax[0].plot([v,w],[i,i],color=col,lw=2);ax[0].scatter(v,i,color=col,s=35,marker='o');ax[0].scatter(w,i,facecolors='white',edgecolors=col,s=55,marker='s')
ax[0].set_yticks(range(len(roster)),[site+' '+d.split('|')[1] for site,d in roster],fontsize=8);ax[0].set_xlabel('Contraction score, mean log1p(CP10k)');ax[0].set_title('S17 Fibro: original vs equal-sample weighting\nFilled circle: original; open square: sensitivity',fontsize=10);ax[0].grid(axis='x',alpha=.2)
ax[1].bar(['Not testable','P floor >= 0.05','P < 0.05 attainable'],[300,250,250],color=['#AAA','#E2AC54','#4C899D'])
for i,n in enumerate([300,250,250]):ax[1].text(i,n+5,str(n),ha='center')
ax[1].set_ylim(0,350);ax[1].set_ylabel('Registered comparisons (overlapping analyses)');ax[1].tick_params(axis='x',labelrotation=25,labelsize=8);ax[1].set_title('Why 800 registered comparisons\nprovide much less independent evidence',fontsize=10)
fig.tight_layout();fig.savefig(O/'figures/diagnostic_overview.png',dpi=170);fig.savefig(O/'figures/diagnostic_overview.pdf');plt.close(fig)
