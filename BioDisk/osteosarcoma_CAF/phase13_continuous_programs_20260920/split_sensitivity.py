from pathlib import Path
import itertools,json
import numpy as np,pandas as pd,anndata as ad
from scipy.stats import rankdata
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');B=P/'phase10_batch_gates_20260919';V=P/'phase11_state_reproducibility_20260920';O=P/'phase13_continuous_programs_20260920'
programs=['Stress_core','Cell_cycle','Contractile','Inflammatory','Matrix_remodeling','Antigen_presentation']
z=np.load(V/'resources/target_gene_expression.npz');genes=z['genes'].astype(str).tolist();ids=z['cell_ids'].astype(str)
a=ad.read_h5ad(B/'objects/unified_strict_integrated.h5ad',backed='r');ob=a.obs.loc[ids].copy();a.file.close()
raw=z['raw'];cor=z['corrected'].copy();n=np.expm1(cor)*ob.corrected_library_counts.to_numpy(float)[:,None]/1e4;cor[np.isfinite(n)&(n<1e-8)]=0;guard=cor.copy();guard[np.isfinite(n)&(n<1)]=0
panels=json.loads((B/'resources/programs.json').read_text());rng=np.random.default_rng(20260920);partitions={};rows=[]
for p in programs:
    gs=sorted(set(panels[p]));parts=[]
    for aa in itertools.combinations(gs,len(gs)//2):
        bb=tuple(g for g in gs if g not in aa)
        if len(aa)==len(bb) and aa>bb:continue
        parts.append((list(aa),list(bb)))
    if len(parts)>25:parts=[parts[i] for i in sorted(rng.choice(len(parts),25,replace=False))]
    partitions[p]=parts
    for donor,ii in ob.groupby('donor_record',observed=True).indices.items():
        meta=ob.iloc[ii];ds=str(meta.dataset.iloc[0])
        if ds=='GSE250015':continue
        cc=np.column_stack([np.ones(len(ii)),rankdata(np.log1p(meta.full_library_counts)),rankdata(np.log1p(meta.n_genes))]);proj=np.linalg.pinv(cc)
        for kind,mat in [('raw',raw),('corrected',cor),('corrected_count1',guard)]:
            x=mat[ii];present={g for g in gs if np.isfinite(x[:,genes.index(g)]).all()};cv=len(present)/len(gs)
            for k,(aa,bb) in enumerate(parts):
                aa=[g for g in aa if g in present];bb=[g for g in bb if g in present]
                active=lambda gg:sum((x[:,genes.index(g)]>0).sum()>=10 and np.std(x[:,genes.index(g)])>1e-12 for g in gg)
                ok=len(ii)>=50 and cv>=.8 and min(len(aa),len(bb),active(aa),active(bb))>=2
                value=np.nan
                if ok:
                    vx=rankdata(x[:,[genes.index(g) for g in aa]].mean(1));vy=rankdata(x[:,[genes.index(g) for g in bb]].mean(1));rx=vx-cc@(proj@vx);ry=vy-cc@(proj@vy)
                    if min(np.std(rx),np.std(ry))>1e-10:value=float(np.corrcoef(rx,ry)[0,1])
                rows.append(dict(program=p,donor_record=donor,dataset=ds,technology=str(meta.technology.iloc[0]),expression=kind,split=k,evaluable=ok,partial_rho=value))
r=pd.DataFrame(rows);r.to_csv(O/'tables/alternative_splits.tsv.gz',sep='\t',index=False)
summ=[]
for keys,g in r.groupby(['program','donor_record','dataset','technology','expression']):
    q=g[g.evaluable&g.partial_rho.notna()]
    summ.append(dict(zip(['program','donor_record','dataset','technology','expression'],keys),total_splits=len(g),evaluable_splits=len(q),positive_splits=int((q.partial_rho>0).sum()),median_partial_rho=q.partial_rho.median(),min_partial_rho=q.partial_rho.min(),max_partial_rho=q.partial_rho.max()))
pd.DataFrame(summ).to_csv(O/'tables/alternative_split_summary.tsv',sep='\t',index=False)
(O/'resources/alternative_splits.json').write_text(json.dumps({'reason':'verify fixed split dependence after observing initial results; no split selection by outcome','seed':20260920,'rule':'all unique balanced partitions if <=25; otherwise 25 sampled without replacement','partitions':partitions},indent=2))
(O/'SPLIT_SENSITIVITY_COMPLETE').write_text('SPLIT_SENSITIVITY_COMPLETE\n')
print('SPLIT_SENSITIVITY_COMPLETE',len(rows))
