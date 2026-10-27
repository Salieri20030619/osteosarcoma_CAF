from pathlib import Path
import json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=Path(__file__).resolve().parent;T=O/'tables';F=O/'figures';F.mkdir(exist_ok=True)
r=pd.read_csv(T/'donor_program_coherence.tsv',sep='\t');alt=pd.read_csv(T/'alternative_splits.tsv.gz',sep='\t');parts=json.loads((O/'disjoint_gene_panels.json').read_text());other=json.loads((O/'alternative_splits.json').read_text())
assert int((O/'exit_code').read_text())==0
assert len(r)==14*17*3 and not r.duplicated(['donor_record','program','expression']).any()
for p,s in parts.items():assert not set(s['A'])&set(s['B'])
for p,ss in other['partitions'].items():
    for a,b in ss:assert not set(a)&set(b) and set(a+b)==set(parts[p]['A']+parts[p]['B'])
assert r.partial_rho.dropna().between(-1,1).all()
assert alt.partial_rho.dropna().between(-1,1).all()
assert ((r[r.evaluable].cells>=50)&(r[r.evaluable].coverage>=.8)).all()
programs=['Stress_core','Matrix_remodeling','Contractile','Cell_cycle','Inflammatory','Antigen_presentation']
q=r[r.dataset.ne('GSE250015')&r.program.isin(programs)];donors=sorted(q.donor_record.unique());fig,axs=plt.subplots(1,3,figsize=(16,7))
for ax,kind in zip(axs,['raw','corrected','corrected_count1']):
    z=q[q.expression.eq(kind)].copy();z.loc[~z.evaluable,'partial_rho']=np.nan
    m=z.pivot(index='donor_record',columns='program',values='partial_rho').reindex(index=donors,columns=programs)
    cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad('#dddddd');im=ax.imshow(m,cmap=cmap,vmin=-.8,vmax=.8,aspect='auto')
    for i in range(len(donors)):
        for j in range(len(programs)):
            v=m.iloc[i,j];ax.text(j,i,'NA' if pd.isna(v) else f'{v:.2f}',ha='center',va='center',fontsize=7,color='white' if pd.notna(v) and abs(v)>.6 else 'black')
    ax.set_xticks(range(6),[p.replace('_','\n') for p in programs],fontsize=8)
    ax.set_yticks(range(len(donors)),[d.split('|')[-1].replace('GSE152048:','').replace('GSE162454__','').replace('GSE318840__Osteosarcoma_','FLEX:') for d in donors],fontsize=8)
    ax.set_title(kind.replace('_',' '))
fig.suptitle('Disjoint gene-panel coherence | within-donor partial Spearman rho\nGray: insufficient cells, detection, or coverage; not zero correlation',fontsize=12)
fig.tight_layout(rect=[0,0,1,.94]);fig.savefig(F/'program_coherence.png',dpi=180);plt.close(fig)
paired=pd.read_csv(T/'paired_expression_coherence.tsv',sep='\t');out=[]
for p,g in paired[paired.all_three_evaluable].groupby('program'):
    out.append(dict(program=p,same_evaluable_records=len(g),raw_median=g.raw.median(),corrected_median=g.corrected.median(),count1_median=g.corrected_count1.median(),all_three_positive=int(((g.raw>0)&(g.corrected>0)&(g.corrected_count1>0)).sum())))
pd.DataFrame(out).to_csv(T/'same_record_expression_summary.tsv',sep='\t',index=False)
qa={'all_checks_passed':True,'disjoint_panels_verified':True,'donor_program_rows':len(r),'alternative_split_rows':len(alt),'correlations_bounded':True,'missing_not_zero':True,'no_P_values':True}
(O/'REVIEW_QA.json').write_text(json.dumps(qa,indent=2))
print(pd.DataFrame(out).query('program in @programs').to_string(index=False))
for p in ['Stress_core','Matrix_remodeling']:
    z=q[(q.program==p)&q.expression.eq('corrected')&q.evaluable]
    print(p,z.assign(FLEX=z.technology.str.contains('FLEX')).groupby('FLEX').partial_rho.agg(['count','median','min','max']).to_string())
