from pathlib import Path
import json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=Path(__file__).resolve().parent;T=O/'tables';F=O/'figures';F.mkdir(exist_ok=True)
def read(n):
    try:return pd.read_csv(T/n,sep='\t')
    except pd.errors.EmptyDataError:return pd.DataFrame()
assert int((O/'exit_code').read_text())==0
fold=read('fold_summary.tsv');attr=read('assignment_attrition.tsv');cand=read('heldout_candidate_results.tsv');features=read('training_selected_features.tsv')
assert len(fold)==4 and fold.test_cells.sum()==3577
sets=[]
for study in fold.heldout_study:
    tr=read(study+'_training_labels.tsv.gz');te=read(study+'_heldout_assignments.tsv.gz')
    assert not set(tr.cell_id)&set(te.cell_id)
    assert not set(tr.donor_record)&set(te.donor_record)
    assert te.dataset.eq(study).all() and tr.dataset.ne(study).all()
    assert len(tr)+len(te)==3577
    assert (te.accepted70==((te.probability>=.7)&(te.distance<=te.radius95))).all()
    assert te.cell_id.is_unique and tr.cell_id.is_unique
    sets.extend(te.cell_id)
assert len(sets)==len(set(sets))==3577
assert (attr.accepted+attr.rejected==attr.total).all()
for _,g in attr.groupby(['fold','donor_record']):assert g.sort_values('threshold').accepted.diff().dropna().le(0).all()
fig,axs=plt.subplots(1,2,figsize=(12,4.5))
names=fold.heldout_study.tolist();x=np.arange(len(names))
for i,cut in enumerate([.5,.7,.9]):
    ss=attr[attr.threshold.eq(cut)].groupby('fold')[['accepted','total']].sum().reindex(names)
    axs[0].bar(x+(i-1)*.24,ss.accepted/ss.total,width=.24,label=f'p >= {cut:.2f}')
axs[0].set_xticks(x,names,rotation=20);axs[0].set_ylim(0,1);axs[0].set_ylabel('Assigned fraction (distance gate also applied)');axs[0].set_title('Held-out study assignment');axs[0].legend(fontsize=8)
axs[1].bar(names,fold.candidate_program_pairs,color='#467c9b');axs[1].set_xticks(x,names,rotation=20);axs[1].set_ylabel('Training-selected cluster/program pairs');axs[1].set_title('Candidates selected before held-out evaluation')
fig.tight_layout();fig.savefig(F/'holdout_overview.png',dpi=180);plt.close(fig)
qa={'all_checks_passed':True,'folds':4,'heldout_cells_once_each':3577,'train_test_records_disjoint':True,'accepted_flags_recomputed':True,'rejection_sensitivity_monotonic':True}
(O/'REVIEW_QA.json').write_text(json.dumps(qa,indent=2))
print(fold.to_string(index=False));print(attr.to_string(index=False))
print('CANDIDATES');print(cand.to_string(index=False) if len(cand) else 'No training-eligible candidate program pairs')
