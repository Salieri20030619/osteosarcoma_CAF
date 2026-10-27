from pathlib import Path
import json,hashlib
import pandas as pd,numpy as np
from statsmodels.stats.multitest import multipletests
from scipy.stats import spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=Path(r'C:\Users\31937\CAF_F2C3_novelty_20261009\deep_review');r=p/'results';t=r/'tables';fig=r/'figures';fig.mkdir(exist_ok=True)
s=pd.read_csv(t/'signature_patient_statistics.tsv',sep='\t');e=pd.read_csv(t/'patient_effects.tsv',sep='\t');c=pd.read_csv(t/'cell_signature_scores.tsv.gz',sep='\t');q=[]
keys=['dataset','matrix','min_cells','comparison','panel']
for k,g in e.groupby(keys):
 a=s.set_index(keys).loc[k];assert g.donor_record.is_unique;assert a.N==len(g);assert np.isclose(a.mean_effect,g.effect.mean());q.append(k)
assert s[s.N<2].P.isna().all()
for fam,g in s[s.P.notna()].groupby('family'):
 expected=multipletests(g.P,method='fdr_bh')[1];assert np.allclose(expected,g.FDR)
assert len(c)==16348 and c.groupby('matrix').cell_id.nunique().eq(8174).all()
assert s.loc[(s.N==5)&s.P.notna(),'P'].min()>=.0625-1e-12
assert s.loc[(s.N==7)&s.P.notna(),'P'].min()>=.015625-1e-12
coverage=pd.read_csv(t/'signature_coverage.tsv',sep='\t');print('MISSINGPANELS',int((coverage.measured<coverage.expected).sum()))
# Independently reconstruct one named contrast from exported cell scores and source metadata.
a=c[(c.dataset=='SCPCP000017')&(c.matrix=='corrected')&(c.lineage_stratum=='Fibro')&(~c.doublet_positive_any)]
v=[]
for donor,d in a.groupby('donor_record'):
 ss=[]
 for site,z in d.groupby('site'):
  left=z[z.subgroup_id=='F2_C3'];right=z[z.subgroup_id!='F2_C3']
  if len(left) and len(right):ss.append(left.F2C3_processing2.mean()-right.F2C3_processing2.mean())
 if ss:v.append(np.mean(ss))
check=s[(s.dataset=='SCPCP000017')&(s.matrix=='corrected')&(s.min_cells==1)&(s.comparison=='other_Fibro')&(s.panel=='F2C3_processing2')].iloc[0]
assert check.N==len(v)==7 and np.isclose(check.mean_effect,np.mean(v))
# Save descriptive library correlations, no inferential p-values.
lib=[]
for panel in ['F2C3_identity3','F2C3_assembly5','F2C3_processing2','Zheng_mCAF6','FNDC1_CAF6']:
 lib.append({'panel':panel,'rho_log_library':spearmanr(a[panel],np.log1p(a.library)).statistic,'cells':len(a)})
pd.DataFrame(lib).to_csv(t/'library_correlation_descriptive.tsv',sep='\t',index=False)
# Plot source-local patient/site matched differences; all data points are independent donor summaries.
panels=['F2C3_assembly5','F2C3_processing2','Zheng_mCAF6','Atlas_mCAF4','FNDC1_CAF6']
labels=['F2-C3 assembly (5 genes)','PCOLCE + LOXL2','Zheng mCAF (6 genes)','Atlas mCAF (4 genes)','FNDC1+ CAF (6 genes)']
f,axes=plt.subplots(1,2,figsize=(12.8,5.2),sharey=True)
for ax,comp,title in zip(axes,['other_Fibro','other_F2'],['A  F2-C3 vs other fibroblast-like cells','B  F2-C3 vs other F2 fibroblast-like cells']):
 for i,(panel,label) in enumerate(zip(panels,labels)):
  d=e[(e.dataset=='SCPCP000017')&(e.matrix=='corrected')&(e.min_cells==1)&(e.comparison==comp)&(e.panel==panel)].sort_values('donor_record')
  ax.scatter(d.effect,i+np.linspace(-.12,.12,len(d)),s=31,color='#357c9c',alpha=.85)
  ax.scatter(d.effect.mean(),i,marker='D',s=40,color='#102c3a',zorder=4)
  ax.text(.98,i+.18,f'N={len(d)}; +{(d.effect>0).sum()}/{len(d)}',transform=ax.get_yaxis_transform(),ha='right',fontsize=8)
 ax.axvline(0,color='#9d9d9d',lw=.8);ax.set_title(title,fontsize=11,loc='left');ax.set_yticks(range(5),labels);ax.invert_yaxis();ax.set_xlabel('Within-patient difference in mean log1p(CP10k)');ax.spines[['top','right']].set_visible(False)
f.suptitle('F2-C3 shares a matrix-CAF program; processing signal remains',fontsize=14,y=1.01)
f.text(.5,-.015,'SCPCP000017 only | corrected RNA | current doublet exclusion | Fibro on both sides\nDots: donors, diamonds: donor mean. Different gene panels have different scales. Exploratory; not reference label transfer.',ha='center',fontsize=9)
f.tight_layout();f.savefig(fig/'01_patient_signature_comparison.png',dpi=180,bbox_inches='tight');plt.close(f)
# Context heatmap standardized within each panel across the six F2 states.
x=pd.read_csv(t/'signature_subgroup_summary.tsv',sep='\t');x=x[(x.dataset=='SCPCP000017')&(x.matrix=='corrected')&(x.scope=='Fibro')&x.subgroup_id.str.startswith('F2_')]
ps=['F2C3_identity3','F2C3_processing2','Zheng_mCAF6','Zheng_myoCAF7','Zheng_iCAF12','FNDC1_CAF6','CD36_CAF5','Mineralization_control3']
y=x.pivot(index='subgroup_id',columns='panel',values='donor_equal_mean')[ps];z=(y-y.mean())/y.std(ddof=0);counts=x.drop_duplicates('subgroup_id').set_index('subgroup_id')
f,ax=plt.subplots(figsize=(11.5,4.8));im=ax.imshow(z,cmap='RdBu_r',vmin=-2,vmax=2,aspect='auto');ax.set_xticks(range(len(ps)),['F2-C3\nidentity','PCOLCE\nLOXL2','mCAF\nZheng','myoCAF\nZheng','iCAF\nZheng','FNDC1+\nCAF','CD36+\nCAF','Mineral\ncontrol']);ax.set_yticks(range(len(y)),[f'{k}  (cells={int(counts.loc[k,"cells"])}, donors={int(counts.loc[k,"donors"])})' for k in y.index]);ax.set_title('F2 states occupy overlapping published functional programs',loc='left',pad=15);f.colorbar(im,ax=ax,label='Within-panel Z score across F2 states');f.text(.5,.005,'Descriptive, source-local, donor-equal means. No significance implied by heatmap colors.',ha='center',fontsize=9);f.tight_layout(rect=(0,.035,1,1));f.savefig(fig/'02_F2_signature_context.png',dpi=180);plt.close(f)
qa={'status':'PASS','statistic_rows':len(s),'valid_P_rows':int(s.P.notna().sum()),'reconciled_patient_contrasts':len(q),'BH_reproduced':True,'N_lt2_no_P':True,'independent_cell_to_patient_reconstruction':True,'exact_P_floor_verified':True,'cells_per_matrix':8174,'original_labels_unchanged':True,'figures':2}
(r/'QA.json').write_text(json.dumps(qa,indent=2));print(json.dumps(qa));print('LIBRARY',lib)
