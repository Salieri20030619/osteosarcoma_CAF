from pathlib import Path
import json,pandas as pd,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');D=P/'phase44_patient_mechanism_20261001';L=P/'phase45_LR_expression_20261001'
g=pd.read_csv(D/'delivery/main_selected_gene_patient_stability.tsv',sep='\t');fig,axes=plt.subplots(1,2,figsize=(10,4.2))
for ax,(layer,cls,genes,color,title) in zip(axes,[('supported','F3',['COL1A1','LOX','MMP14'],'#176D71','Sup F3'),('expanded','F4',['RUNX2','ANGPT1','SMAD3'],'#7846A0','Exp F4')]):
 z=g[(g.layer==layer)&g.contrast.eq(cls+'_vs_rest')].set_index('symbol').loc[genes]
 yy=np.arange(len(z));value=z.mean_delta.to_numpy();lower=z.patient_leaveout_mean_min.to_numpy();upper=z.patient_leaveout_mean_max.to_numpy()
 ax.errorbar(value,yy,xerr=np.vstack([value-lower,upper-value]),fmt='o',color=color,ecolor=color,capsize=4,markersize=7)
 ax.axvline(0,color='#777777',lw=.8);ax.set_yticks(yy,genes);ax.invert_yaxis();ax.set_xlabel('Mean delta');ax.set_title(title+' (N='+str(int(z.N.iloc[0]))+')');ax.spines[['top','right']].set_visible(False);ax.set_xlim(min(0,lower.min())-.2,upper.max()+.5)
fig.suptitle('Patient leave-out direction');fig.text(.5,.015,'Points: full mean; ranges: leave-one-patient-out means; delta: log2(CP10k+0.1)',ha='center',fontsize=9);fig.tight_layout(rect=[0,.06,1,.93]);fig.savefig(D/'figures/selected_gene_patient_leaveout.png',dpi=180);fig.savefig(D/'figures/selected_gene_patient_leaveout.pdf');plt.close(fig)
# Preserve full clinical model signals and their sensitivity side by side, not select the convenient reference.
sens=pd.read_csv(D/'clinical_model_sensitivity/tables/all_adjusted_HC3_t_reference_sensitivity.tsv',sep='\t');sens[sens.primary_P.lt(.05)].to_csv(D/'delivery/clinical_nominal_primary_and_t_sensitivity.tsv',sep='\t',index=False)
tf=pd.read_csv(D/'delivery/main_all_nominal_TF_response.tsv',sep='\t');tf[(tf.layer.eq('supported')&tf.contrast.eq('F3_vs_rest'))|(tf.regulator.isin(['CBFB','TGIF1','LYL1','NRL','CRX','ALX4']))][['layer','contrast','regulator','N','defined_targets','measured_targets','mean_delta','P','FDR']].sort_values('P').head(25).to_csv(D/'delivery/TF_regulon_coverage_review.tsv',sep='\t',index=False)
print('TF_COVERAGE',pd.read_csv(D/'delivery/TF_regulon_coverage_review.tsv',sep='\t').to_string(index=False))
print('LR_ACCEPTED', (L/'QA_RESUME_STATUS.json').read_text());print('CORE_ACCEPTED',(D/'delivery/CORE_ACCEPTED_STATUS.json').read_text())