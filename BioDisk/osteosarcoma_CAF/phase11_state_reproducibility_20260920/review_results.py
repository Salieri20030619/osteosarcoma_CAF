from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

O=Path(__file__).resolve().parent;T=O/'tables';F=O/'figures';F.mkdir(exist_ok=True)
read=lambda f:pd.read_csv(T/f,sep='\t')
e=read('patient_program_effects.tsv');elig=read('patient_comparison_eligibility.tsv')
tr=read('all_cluster_transitions.tsv');m=read('best_cluster_transitions.tsv');s=read('stress_sensitivity_metrics.tsv')
assert int((O/'exit_code').read_text())==0 and int((O/'stress_exit_code').read_text())==0
assert len(elig)==14*11 and set(elig.cluster)==set(range(11))
assert ((elig.in_cells+elig.out_cells)==elig.groupby('donor_record').in_cells.transform('sum')).all()
assert elig.groupby('cluster').in_cells.sum().sum()==3976
assert e[['raw_delta','corrected_delta','coverage']].notna().all().all()
assert (e.coverage.between(0,1)).all()
for (model,c),g in tr.groupby(['model','cluster']):
    assert g.cells.sum()==g.remaining_primary_cells.iloc[0]
    assert np.isclose(g.recall.sum(),1)
assert (m.Jaccard.dropna().between(0,1)).all()
for name in s.model:
    z=read(name+'_labels.tsv.gz');assert len(z)==3976 and z.cell_id.nunique()==3976

programs=['Stress_core','Contractile','Matrix_remodeling','Inflammatory','Antigen_presentation','Cell_cycle']
q=e[(e.cluster==0)&e.primary20&e.dataset.ne('GSE250015')&e.program.isin(programs)]
donors=sorted(q.donor_record.unique());fig,axs=plt.subplots(1,2,figsize=(13,4.8))
v=max(abs(q.raw_delta).max(),abs(q.corrected_delta).max())
for ax,kind in zip(axs,['raw','corrected']):
    matrix=q.pivot(index='donor_record',columns='program',values=kind+'_delta').reindex(index=donors,columns=programs)
    im=ax.imshow(matrix,cmap='RdBu_r',vmin=-v,vmax=v,aspect='auto')
    for i in range(len(donors)):
        for j in range(len(programs)):ax.text(j,i,f'{matrix.iloc[i,j]:.2f}',ha='center',va='center',fontsize=8)
    ax.set_xticks(range(len(programs)),[p.replace('_','\n') for p in programs],rotation=35,ha='right',fontsize=8)
    ax.set_yticks(range(len(donors)),[d.split('|')[-1].replace('GSE162454__','') for d in donors],fontsize=8)
    ax.set_title('U0 within-donor effect | '+kind)
fig.subplots_adjust(bottom=.3,wspace=.4,right=.91);cax=fig.add_axes([.93,.3,.015,.55]);fig.colorbar(im,cax=cax,label='Mean log1p(CP10k) difference')
fig.savefig(F/'01_U0_patient_programs.png',dpi=180);plt.close(fig)

models=['hvg1500_pc30','corrected_counts','exclude_stress','exclude_cycle','exclude_stress_cycle']
q=m[m.model.isin(models)];fig,axs=plt.subplots(1,2,figsize=(13,5))
for ax,metric in zip(axs,['Jaccard','recall']):
    matrix=q.pivot(index='cluster',columns='model',values=metric).reindex(index=range(11),columns=models)
    im=ax.imshow(matrix,cmap='viridis',vmin=0,vmax=1,aspect='auto')
    for i in range(11):
        for j in range(len(models)):ax.text(j,i,f'{matrix.iloc[i,j]:.2f}',ha='center',va='center',fontsize=8,color='white' if matrix.iloc[i,j]<.5 else 'black')
    ax.set_yticks(range(11),['U'+str(i) for i in range(11)])
    ax.set_xticks(range(len(models)),[n.replace('_','\n') for n in models],fontsize=8)
    ax.set_title('Best-match '+metric)
fig.tight_layout();fig.savefig(F/'02_boundary_sensitivity.png',dpi=180);plt.close(fig)
qa={'all_checks_passed':True,'strict_cells':3976,'clusters':11,'donor_records':14,'transition_count_conservation':True,'targeted_refits':len(s),'biological_pvalues':False,'original_labels_modified':False}
(O/'REVIEW_QA.json').write_text(json.dumps(qa,indent=2))
print(s.to_string(index=False))
print(m[m.model.str.startswith('exclude_')].to_string(index=False))
print(json.dumps(qa))
