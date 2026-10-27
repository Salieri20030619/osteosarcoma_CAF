"""SYNTAX_CHECK: python -m py_compile continuous.py
RUN_COMMAND: python -u continuous.py
EXPECTED_SUCCESS_MARKER: ANALYSIS_COMPLETE and exit_code=0
"""
from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
import json,hashlib
import numpy as np,pandas as pd,anndata as ad
from scipy.stats import rankdata
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');B=P/'phase10_batch_gates_20260919';V=P/'phase11_state_reproducibility_20260920';O=P/'phase13_continuous_programs_20260920'
for d in ['tables','resources','logs']: (O/d).mkdir(parents=True,exist_ok=True)
primary=['Stress_core','Cell_cycle','Contractile','Inflammatory','Matrix_remodeling','Antigen_presentation']
protocol={'design':'label-free within-donor coherence of fixed panels; no reclustering or U labels used for evaluation',
 'panels':'all 17 existing non-single-marker panels; six primary programs listed separately',
 'primary_programs':primary,'split':'alphabetically sorted gene symbols, alternating into two disjoint halves; frozen before expression evaluation',
 'unit':'donor record; exclude GSE250015 from main summaries; study and technology kept separate',
 'eligibility':'at least50 cells, coverage>=80%, each half >=2 measured genes and >=2 genes detected in >=10 cells and with variable expression',
 'score':'mean gene log1p(CP10k) separately for each half; source full library normalization',
 'correlation':'Spearman half-score rho; partial Spearman uses rank residuals controlling rank(log1p full_library_counts) and rank(log1p n_genes)',
 'ambient':'raw vs DecontX; corrected inferred counts <1e-8 set to zero for numerical residual handling',
 'sensitivity':'corrected inferred counts <1 set to zero; one-gene omission with each surviving half at least2 measured genes',
 'inference':'no cell-level or patient-level P values; cell correlation describes within-patient coherence; no independent subtype claim',
 'effect_display':'rho reported continuously; descriptive flags rho>=0.2 and positive are not significance or subtype acceptance criteria'}
(O/'PROTOCOL.json').write_text(json.dumps(protocol,indent=2))
z=np.load(V/'resources/target_gene_expression.npz');ids=z['cell_ids'].astype(str);genes=z['genes'].astype(str).tolist()
a=ad.read_h5ad(B/'objects/unified_strict_integrated.h5ad',backed='r');ob=a.obs.loc[ids].copy();a.file.close()
assert len(ob)==3976 and ob.index.is_unique
raw=z['raw'].copy();corrected=z['corrected'].copy();den=ob.corrected_library_counts.to_numpy(float)
inferred=np.expm1(corrected)*den[:,None]/1e4
corrected[np.isfinite(inferred)&(inferred<1e-8)]=0
guard=corrected.copy();guard[np.isfinite(inferred)&(inferred<1)]=0
panels=json.loads((B/'resources/programs.json').read_text());panels={k:sorted(set(v)) for k,v in panels.items() if not k.endswith('_marker')}
splits={p:{'A':gs[::2],'B':gs[1::2]} for p,gs in panels.items()}
(O/'resources/disjoint_gene_panels.json').write_text(json.dumps(splits,indent=2))
def rho(x,y):
    if np.std(x)<1e-12 or np.std(y)<1e-12:return np.nan
    return float(np.corrcoef(rankdata(x),rankdata(y))[0,1])
def partial(x,y,cov):
    xx=rankdata(x);yy=rankdata(y);c=np.column_stack([np.ones(len(x)),*[rankdata(cov[:,j]) for j in range(cov.shape[1])]])
    rx=xx-c@np.linalg.lstsq(c,xx,rcond=None)[0];ry=yy-c@np.linalg.lstsq(c,yy,rcond=None)[0]
    if np.std(rx)<1e-10 or np.std(ry)<1e-10:return np.nan
    return float(np.corrcoef(rx,ry)[0,1])
rows=[];omits=[];pairrows=[]
for donor,idx in ob.groupby('donor_record',observed=True).indices.items():
    meta=ob.iloc[idx];ds=str(meta.dataset.iloc[0]);tech=str(meta.technology.iloc[0]);cov=np.column_stack([np.log1p(meta.full_library_counts.to_numpy(float)),np.log1p(meta.n_genes.to_numpy(float))])
    for p,split in splits.items():
        for kind,mat in [('raw',raw),('corrected',corrected),('corrected_count1',guard)]:
            x=mat[idx];present={g for g in panels[p] if np.isfinite(x[:,genes.index(g)]).all()};coverage=len(present)/len(panels[p])
            aa=[g for g in split['A'] if g in present];bb=[g for g in split['B'] if g in present]
            active=lambda gs:sum(int((x[:,genes.index(g)]>0).sum()>=10 and np.std(x[:,genes.index(g)])>1e-12) for g in gs)
            na=active(aa);nb=active(bb);ok=len(idx)>=50 and coverage>=.8 and min(len(aa),len(bb),na,nb)>=2
            xs=x[:,[genes.index(g) for g in aa]].mean(1) if aa else np.zeros(len(idx));ys=x[:,[genes.index(g) for g in bb]].mean(1) if bb else np.zeros(len(idx))
            rec=dict(donor_record=donor,dataset=ds,technology=tech,program=p,primary_program=p in primary,expression=kind,cells=len(idx),coverage=coverage,measured_A=len(aa),measured_B=len(bb),active_A=na,active_B=nb,evaluable=ok,rho=rho(xs,ys) if aa and bb else np.nan,partial_rho=partial(xs,ys,cov) if aa and bb else np.nan)
            rows.append(rec)
            if not ok:continue
            for gene in aa+bb:
                aaa=[g for g in aa if g!=gene];bbb=[g for g in bb if g!=gene]
                if min(len(aaa),len(bbb))<2:continue
                vx=x[:,[genes.index(g) for g in aaa]].mean(1);vy=x[:,[genes.index(g) for g in bbb]].mean(1)
                omits.append(dict(donor_record=donor,dataset=ds,technology=tech,program=p,expression=kind,omitted_gene=gene,partial_rho=partial(vx,vy,cov)))
            if kind in ['raw','corrected']:
                gg=aa+bb
                for i,g in enumerate(gg):
                    for h in gg[i+1:]:
                        vx=x[:,genes.index(g)];vy=x[:,genes.index(h)]
                        if (vx>0).sum()<10 or (vy>0).sum()<10:continue
                        pairrows.append(dict(donor_record=donor,dataset=ds,program=p,expression=kind,gene1=g,gene2=h,partial_rho=partial(vx,vy,cov)))
    print('DONOR_COMPLETE',donor,len(idx),flush=True)
def save(d,n):d.to_csv(O/'tables'/n,sep='\t',index=False)
r=pd.DataFrame(rows);save(r,'donor_program_coherence.tsv');loo=pd.DataFrame(omits);save(loo,'leave_gene_out_coherence.tsv.gz');save(pd.DataFrame(pairrows),'gene_pair_coherence.tsv.gz')
summ=[]
for scope in ['all_records','exclude_GSE250015']+sorted(ob.dataset.unique())+['technology:'+t for t in sorted(ob.technology.unique())]:
    sub=r if scope=='all_records' else r[r.dataset.ne('GSE250015')] if scope=='exclude_GSE250015' else r[r.technology.eq(scope.split(':',1)[1])&r.dataset.ne('GSE250015')] if scope.startswith('technology:') else r[r.dataset.eq(scope)]
    for (p,kind),g in sub.groupby(['program','expression']):
        q=g[g.evaluable&g.partial_rho.notna()]
        summ.append(dict(scope=scope,program=p,expression=kind,available_records=len(g),evaluable_records=len(q),studies=q.dataset.nunique(),positive_records=int((q.partial_rho>0).sum()),rho_ge02_records=int((q.partial_rho>=.2).sum()),median_partial_rho=q.partial_rho.median(),min_partial_rho=q.partial_rho.min(),max_partial_rho=q.partial_rho.max()))
save(pd.DataFrame(summ),'coherence_summary.tsv')
# Paired expression-sensitivity summaries use only the exact same evaluable donor set.
wide=r.pivot(index=['donor_record','dataset','technology','program'],columns='expression',values=['partial_rho','evaluable']);paired=[]
for ix,q in wide.iterrows():
    if ix[1]=='GSE250015':continue
    valid=all(bool(q[('evaluable',k)]) and np.isfinite(q[('partial_rho',k)]) for k in ['raw','corrected','corrected_count1'])
    vals={k:float(q[('partial_rho',k)]) for k in ['raw','corrected','corrected_count1']}
    paired.append(dict(donor_record=ix[0],dataset=ix[1],technology=ix[2],program=ix[3],all_three_evaluable=valid,**vals))
save(pd.DataFrame(paired),'paired_expression_coherence.tsv')
assert len(r)==len(ob.donor_record.unique())*len(panels)*3
assert not any(set(s['A'])&set(s['B']) for s in splits.values())
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
(O/'resources/input_sha256.json').write_text(json.dumps({str(p):sha(p) for p in [V/'resources/target_gene_expression.npz',B/'objects/unified_strict_integrated.h5ad',B/'resources/programs.json']},indent=2))
(O/'STATUS.json').write_text(json.dumps({'status':'CONTINUOUS_PROGRAMS_COMPUTED_REVIEW_PENDING','cells':3976,'donor_records':14,'main_excludes_GSE250015':True,'panels':17,'independent_subtype_validation':False,'P_values_computed':False},indent=2))
(O/'ANALYSIS_COMPLETE').write_text('ANALYSIS_COMPLETE\n')
