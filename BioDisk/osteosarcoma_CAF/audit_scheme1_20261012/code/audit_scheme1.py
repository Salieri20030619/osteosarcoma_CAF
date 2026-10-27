from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']: os.environ[k]='2'
import json,hashlib,itertools,math,traceback
from datetime import datetime,timezone
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse

P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF')
S=P/'phase73_F134_bone_lung_20261009/strict_v2'
E=P/'phase80_scheme1_anatomical_20261012'
O=P/'audit_scheme1_20261012'
for d in ['tables','QA','delivery','resources']: (O/d).mkdir(parents=True,exist_ok=True)
def read(f,**kw):return pd.read_csv(f,sep='\t',low_memory=False,**kw)
def save(d,n):d.to_csv(O/'tables'/n,sep='\t',index=False)
def dump(x,p):p.write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str)+'\n',encoding='utf-8')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def bh(p):
    p=np.asarray(p,float);ix=np.argsort(p,kind='stable');q=np.minimum.accumulate((p[ix]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]
    out=np.empty(len(p));out[ix]=np.minimum(q,1);return out
def exact_p(x,y,paired):
    if paired:
        ok=np.isfinite(x)&np.isfinite(y);v=y[ok]-x[ok]
        if len(v)<2 or not np.any(v):return np.nan
        obs=abs(v.mean());null=[abs(np.mean(v*np.array(sign))) for sign in itertools.product([-1,1],repeat=len(v))]
    else:
        x=x[np.isfinite(x)];y=y[np.isfinite(y)]
        if min(len(x),len(y))<2:return np.nan
        v=np.r_[x,y]
        if np.all(v==v[0]):return np.nan
        obs=abs(y.mean()-x.mean());null=[]
        for idx in itertools.combinations(range(len(v)),len(y)):
            mask=np.zeros(len(v),bool);mask[list(idx)]=True
            null.append(abs(v[mask].mean()-v[~mask].mean()))
    return float(np.mean(np.asarray(null)>=obs-1e-12))

eligible=read(E/'tables/scheme1_all_eligible_comparisons.tsv')
clinical=read(S/'resources/clinical_sample_manifest.tsv',keep_default_na=False)
cells={layer:read(S/'tables'/f'{layer}_fixed_cells_with_sites.tsv.gz',keep_default_na=False) for layer in ['supported','expanded']}
tracked=[S/'code/analyze_bone_lung_v2.py',S/'resources/clinical_sample_manifest.tsv',P/'CAF_SCHEME1_CURRENT.json']
hashes={str(p):sha(p) for p in tracked}
dump(hashes,O/'resources/BEFORE_SHA256.json')
source_stats=[];site_errors=[];target_meta={}
for layer,ob in cells.items():
    cm=clinical.set_index('sample')
    mismatch=ob.strict_site.ne(ob['sample'].map(cm.strict_site_group))
    site_errors.append(dict(layer=layer,cell_site_mismatches=int(mismatch.sum()),cells=len(ob)))
    for (ds,site),sg in clinical.groupby(['dataset','strict_site_group']):
        z=ob[ob['sample'].isin(sg['sample'])]
        if site not in ['Bone','Lung']:continue
        source_stats.append(dict(layer=layer,dataset=ds,site=site,source_samples=len(sg),source_patient_records=sg.donor_record.nunique(),full_called_cells=sg.full_called_cells.sum(),candidate_cells=len(z),patients_with_candidates=z.donor_record.nunique(),F1_cells=int(z.analysis_class.eq('F1').sum()),F3_cells=int(z.analysis_class.eq('F3').sum()),F4_cells=int(z.analysis_class.eq('F4').sum())))
    target_meta[layer]=ob
save(pd.DataFrame(source_stats),'source_site_retention.tsv');save(pd.DataFrame(site_errors),'grouping_checks.tsv')

families=[];contexts=[];probes=[];precision=[];effect_checks=[]
gene_paths=sorted((S/'tables').glob('*__all_gene_statistics.tsv.gz'))
gene_columns=['context_id','family','sensitivity','layer','dataset','target','level','design','feature','symbol','P','FDR','mean_lung_minus_bone','bone_N','lung_N','test']
all_eligible_stats={}
for path in gene_paths:
    df=read(path,usecols=gene_columns);df=df[df.sensitivity.eq('all_cells')].copy()
    for fam,g in df.groupby('family'):
        good=g.P.notna(); p=g.P.fillna(1).to_numpy();fullq=bh(p)
        err=float(np.max(np.abs(fullq[good.to_numpy()]-g.loc[good,'FDR'].to_numpy()))) if good.any() else 0
        qvalid=bh(g.loc[good,'P']) if good.any() else np.array([])
        families.append(dict(family=fam,all_hypothesis_rows=len(g),valid_P_rows=int(good.sum()),untestable_rows=int((~good).sum()),untestable_fraction=float((~good).mean()),original_min_FDR=float(g.FDR.min()) if good.any() else None,valid_only_diagnostic_min_FDR=float(qvalid.min()) if good.any() else None,original_FDR_hits=int(g.FDR.lt(.05).sum()),valid_only_diagnostic_FDR_hits=int((qvalid<.05).sum()),BH_reproduction_max_error=err))
        df.loc[g.loc[good].index,'valid_only_family_q_diagnostic']=qvalid
    for ctx,g in df.groupby('context_id'):
        if ctx not in set(eligible.context_id):continue
        g=g.copy();all_eligible_stats[ctx]=g
        qcontrast=bh(g.loc[g.P.notna(),'P'])
        row=eligible[eligible.context_id.eq(ctx)].iloc[0]
        nx,ny=int(row.bone_N),int(row.lung_N)
        assignments=2**nx if row.design=='paired' else math.comb(nx+ny,ny)
        pbound=(2/assignments if row.design=='paired' or nx==ny else 1/assignments)
        contexts.append(dict(context_id=ctx,layer=row.layer,dataset=row.dataset,target=row.target,design=row.design,bone_N=nx,lung_N=ny,hypothesis_rows=len(g),valid_P_rows=int(g.P.notna().sum()),nominal_rows=int(g.P.lt(.05).sum()),min_observed_P=float(g.P.min()),permutation_assignments=assignments,theoretical_min_P_lower_bound=pbound,original_min_FDR=float(g.FDR.min()),diagnostic_valid_family_min_q=float(g.valid_only_family_q_diagnostic.min()),diagnostic_valid_context_min_q=float(qcontrast.min()),diagnostic_context_FDR_hits=int((qcontrast<.05).sum())))
        with np.load(S/'objects'/f'{ctx}.npz',allow_pickle=True) as z:
            features=z['features'].astype(str);symbols=z['symbols'].astype(str)
            x=z['bone_expression'];y=z['lung_expression']
        gm=g.set_index('feature').loc[features]
        with np.errstate(invalid='ignore'):
            e=np.nanmean(y,axis=0)-np.nanmean(x,axis=0)
        ef=np.abs(e-gm.mean_lung_minus_bone.to_numpy())
        effect_checks.append(dict(context_id=ctx,genes=len(features),finite_effects=int(np.isfinite(ef).sum()),max_effect_error=float(np.nanmax(ef))))
        # Independent exact enumeration for lowest-P, evenly spaced, and predefined marker features.
        candidates=list(g[g.P.notna()].sort_values(['P','feature']).head(4).feature)
        valid_index=np.where(gm.P.notna().to_numpy())[0]
        candidates+=list(features[valid_index[np.linspace(0,len(valid_index)-1,min(10,len(valid_index)),dtype=int)]])
        candidates+=list(features[np.isin(symbols,['ACTA2','COL1A1','GPC3','SPP1','CD74','IL6'])])
        for f in dict.fromkeys(candidates):
            i=int(np.where(features==f)[0][0]);p=exact_p(x[:,i],y[:,i],row.design=='paired');old=gm.loc[f,'P']
            probes.append(dict(context_id=ctx,feature=f,symbol=symbols[i],stored_P=old,independent_P=p,error=abs(p-old) if np.isfinite(p) and np.isfinite(old) else None,NA_agrees=bool(np.isnan(p)==pd.isna(old))))
    del df
save(pd.DataFrame(families),'FDR_family_diagnostics.tsv');save(pd.DataFrame(contexts),'test_resolution_and_FDR.tsv');save(pd.DataFrame(probes),'independent_exact_P_checks.tsv');save(pd.DataFrame(effect_checks),'all_gene_effect_checks.tsv')
print('GENE_STATISTICS_AUDITED',len(probes),flush=True)

# Check both prevalence denominators and sample-to-site labels without changing any labels.
pv=read(S/'tables/patient_site_proportions.tsv.gz');pv=pv[pv.sensitivity.eq('all_cells')]
proportion_checks=[]
for r in pv.itertuples():
    z=target_meta[r.layer];samples=r.samples.split(';');z=z[z['sample'].isin(samples)]
    num=int(z.analysis_class.eq(r.target).sum()) if r.level=='main' else int(z.analysis_internal.eq(r.target).sum())
    den=int(clinical[clinical['sample'].isin(samples)].full_called_cells.sum()) if r.denominator=='full_called' else (len(z) if r.denominator=='layer' else int(z.analysis_class.eq(r.parent).sum()))
    fraction=num/den if den else np.nan
    proportion_checks.append(dict(layer=r.layer,dataset=r.dataset,target=r.target,denominator=r.denominator,donor_record=r.donor_record,site=r.site,numerator_matches=num==r.numerator,denominator_matches=den==r.denominator_cells,fraction_error=abs(fraction-r.fraction) if den else 0))
pc=pd.DataFrame(proportion_checks);save(pc,'proportion_denominator_checks.tsv')
prop=read(S/'tables/all_proportion_statistics.tsv');prop=prop[prop.sensitivity.eq('all_cells')].copy()
prop['mean_difference_percentage_points']=prop.mean_lung_minus_bone*100
save(prop[prop.P.notna()].sort_values('P'),'valid_proportion_effects_and_P.tsv')

# Reconstruct the saved pseudobulk from raw counts, checking every feature in all 10 eligible contexts.
aa=ad.read_h5ad(P/'phase41_full_transcriptome_setup_20261001/objects/all_candidate_full_genes.h5ad')
present=sparse.load_npz(P/'phase41_full_transcriptome_setup_20261001/objects/sample_feature_presence.npz').tocsr()
sr=read(P/'phase41_full_transcriptome_setup_20261001/tables/feature_presence_sample_index.tsv').set_index('sample').presence_row.to_dict()
rawchecks=[];quality=[];featuremiss=[]
for row in eligible.itertuples():
    ob=target_meta[row.layer]
    ob=ob[ob.dataset.eq(row.dataset)&(ob.analysis_class.eq(row.target) if row.level=='main' else ob.analysis_internal.eq(row.target))]
    with np.load(S/'objects'/f'{row.context_id}.npz',allow_pickle=True) as z:
        feats=z['features'].astype(str);gi=aa.var_names.get_indexer(feats)
        for side in ['bone','lung']:
            mat=z[f'{side}_expression'];pat=z[f'{side}_patients'].astype(str)
            for i,donor in enumerate(pat):
                q=ob[ob.donor_record.eq(donor)&ob.strict_site.eq(side.title())]
                ix=aa.obs_names.get_indexer(q.cell_id);assert (ix>=0).all()
                raw=np.asarray(aa.X[ix][:,gi].sum(axis=0)).ravel()
                library=float(aa.X[ix].sum());stored_lib=q.full_library_counts.astype(float).sum()
                samples=q['sample'].unique();valid=np.asarray(present[[sr[s] for s in samples]][:,gi].sum(axis=0)).ravel()==len(samples)
                v=np.log2(raw/library*1e4+.1);v[~valid]=np.nan
                diff=np.abs(v-mat[i]);na_ok=np.array_equal(np.isnan(v),np.isnan(mat[i]))
                rawchecks.append(dict(context_id=row.context_id,site=side,donor_record=donor,cells=len(q),library_UMIs=library,stored_library_error=abs(library-stored_lib),max_expression_error=float(np.nanmax(diff)),missingness_matches=na_ok,measured_features=int(valid.sum())))
                # Complete-library and detection metrics describe precision; no cell threshold is applied.
                counts=q.UMI.astype(float)
                quality.append(dict(context_id=row.context_id,layer=row.layer,dataset=row.dataset,target=row.target,design=row.design,site=side,donor_record=donor,cells=len(q),library_UMIs=library,median_cell_UMIs=float(counts.median()),expressed_features=int(((raw>0)&valid).sum()),measured_features=int(valid.sum()),zero_fraction=float(np.mean(raw[valid]==0)),lineage_composition=json.dumps(q.review_lineage_code.value_counts().to_dict()),CNV_composition=json.dumps(q.CNV_evidence_current.value_counts().to_dict())))
save(pd.DataFrame(rawchecks),'raw_count_reconstruction_checks.tsv');save(pd.DataFrame(quality),'patient_expression_precision.tsv')
print('RAW_INPUTS_RECONSTRUCTED',len(rawchecks),flush=True)

# Program and proportion P are audited from their patient values using a separate enumerator.
programs=read(S/'tables/all_program_statistics.tsv.gz');programs=programs[programs.sensitivity.eq('all_cells')]
program_probes=[];program_families=[]
for fam,g in programs.groupby('family'):
    good=g.P.notna();qfull=bh(g.P.fillna(1));qvalid=bh(g.loc[good,'P']) if good.any() else np.array([])
    program_families.append(dict(family=fam,all_rows=len(g),valid_P_rows=int(good.sum()),original_min_FDR=float(g.FDR.min()) if good.any() else None,valid_only_diagnostic_min_FDR=float(qvalid.min()) if good.any() else None,valid_only_diagnostic_FDR_hits=int((qvalid<.05).sum()),BH_error=float(np.max(np.abs(qfull[good]-g.loc[good,'FDR']))) if good.any() else 0))
for row in eligible.itertuples():
    sub=programs[programs.context_id.eq(row.context_id)].set_index('program')
    with np.load(S/'objects'/f'{row.context_id}.npz',allow_pickle=True) as z:
        names=z['programs'].astype(str);x=z['bone_programs'];y=z['lung_programs']
    chosen=list(sub[sub.P.notna()].sort_values('P').head(5).index)+list(sub.index[sub.index.str.startswith('customer|')])
    for term in dict.fromkeys(chosen):
        idx=int(np.where(names==term)[0][0]);pp=exact_p(x[:,idx],y[:,idx],row.design=='paired');old=sub.loc[term,'P']
        program_probes.append(dict(context_id=row.context_id,program=term,stored_P=old,independent_P=pp,error=abs(pp-old) if np.isfinite(pp) and np.isfinite(old) else None,NA_agrees=np.isnan(pp)==pd.isna(old)))
save(pd.DataFrame(program_probes),'independent_program_P_checks.tsv');save(pd.DataFrame(program_families),'program_FDR_diagnostics.tsv')

propprobes=[]
for r in prop[prop.P.notna()].itertuples():
    v=pv[(pv.layer==r.layer)&(pv.dataset==r.dataset)&(pv.target==r.target)&(pv.denominator==r.denominator)].dropna(subset=['fraction'])
    b=v[v.site.eq('Bone')].set_index('donor_record').fraction;l=v[v.site.eq('Lung')].set_index('donor_record').fraction
    paired=sorted(set(b.index)&set(l.index))
    if r.design=='paired':
        x=b.loc[paired].to_numpy();y=l.loc[paired].to_numpy();perm=2**len(x)
    else:
        x=b.loc[sorted(set(b.index)-set(paired))].to_numpy();y=l.loc[sorted(set(l.index)-set(paired))].to_numpy();perm=math.comb(len(x)+len(y),len(y))
    if perm>65536:continue
    pp=exact_p(x,y,r.design=='paired')
    propprobes.append(dict(layer=r.layer,dataset=r.dataset,target=r.target,denominator=r.denominator,design=r.design,bone_N=len(x),lung_N=len(y),stored_P=r.P,independent_P=pp,error=abs(pp-r.P)))
save(pd.DataFrame(propprobes),'independent_proportion_P_checks.tsv')

# Describe lineage mixing and unequal capture; this audit does not change candidate identity.
lineage=[]
for row in eligible.itertuples():
    ob=target_meta[row.layer]
    donors=read(E/'tables/scheme1_patient_cell_coverage.tsv')
    donors=donors[donors.context_id.eq(row.context_id)]
    ob=ob[ob.dataset.eq(row.dataset)&(ob.analysis_class.eq(row.target) if row.level=='main' else ob.analysis_internal.eq(row.target))]
    for r in donors.itertuples():
        z=ob[ob.donor_record.eq(r.donor_record)&ob.strict_site.eq(r.site.title())]
        for code,n in z.review_lineage_code.value_counts().items():
            lineage.append(dict(context_id=row.context_id,site=r.site,donor_record=r.donor_record,review_lineage=code,cells=n,fraction=n/len(z)))
save(pd.DataFrame(lineage),'lineage_composition_by_patient.tsv')
quality=pd.DataFrame(quality)
save(quality.groupby(['context_id','site']).agg(patients=('donor_record','nunique'),cells=('cells','sum'),min_cells=('cells','min'),median_cells=('cells','median'),patients_with_one_cell=('cells',lambda s:int((s==1).sum())),median_library_UMIs=('library_UMIs','median'),median_expressed_features=('expressed_features','median'),median_zero_fraction=('zero_fraction','median')).reset_index(),'precision_summary.tsv')

checks=dict(audited_at=datetime.now(timezone.utc).isoformat(),audit_only=True,analysis_inputs_unchanged=all(sha(Path(f))==v for f,v in hashes.items()),
    site_assignment_mismatches=sum(r['cell_site_mismatches'] for r in site_errors),
    raw_patient_reconstructions=len(rawchecks),max_raw_expression_error=max(r['max_expression_error'] for r in rawchecks),
    max_library_error=max(r['stored_library_error'] for r in rawchecks),all_missingness_matches=all(r['missingness_matches'] for r in rawchecks),
    gene_P_probes=len(probes),max_gene_P_error=float(pd.DataFrame(probes).error.max()),gene_P_missingness_match=bool(pd.DataFrame(probes).NA_agrees.all()),
    max_gene_effect_error=max(r['max_effect_error'] for r in effect_checks),max_gene_BH_error=max(r['BH_reproduction_max_error'] for r in families),
    program_P_probes=len(program_probes),max_program_P_error=float(pd.DataFrame(program_probes).error.max()),max_program_BH_error=max(r['BH_error'] for r in program_families),
    proportion_denominator_records=len(pc),all_proportion_numerators_match=bool(pc.numerator_matches.all()),all_proportion_denominators_match=bool(pc.denominator_matches.all()),
    max_proportion_fraction_error=float(pc.fraction_error.max()),proportion_exact_P_probes=len(propprobes),max_proportion_P_error=float(pd.DataFrame(propprobes).error.max()),
    diagnostic_valid_family_gene_FDR_hits=sum(r['valid_only_diagnostic_FDR_hits'] for r in families),diagnostic_per_context_gene_FDR_hits=sum(r['diagnostic_context_FDR_hits'] for r in contexts),
    diagnostic_valid_family_program_FDR_hits=sum(r['valid_only_diagnostic_FDR_hits'] for r in program_families),original_inferential_results_overwritten=False)
checks['implementation_checks_passed']=checks['analysis_inputs_unchanged'] and checks['site_assignment_mismatches']==0 and checks['all_missingness_matches'] and checks['all_proportion_numerators_match'] and checks['all_proportion_denominators_match'] and all(checks[k]<1e-9 for k in ['max_raw_expression_error','max_library_error','max_gene_P_error','max_gene_effect_error','max_gene_BH_error','max_program_P_error','max_program_BH_error','max_proportion_fraction_error','max_proportion_P_error'])
dump(checks,O/'QA/AUDIT_CHECKS.json')
(O/'AUDIT_COMPUTATION_COMPLETE').write_text('Diagnostic computations complete; interpretation review pending.\n')
print(json.dumps(checks,ensure_ascii=False,indent=2),flush=True)
print('CONTEXT_DIAGNOSTICS\n'+pd.DataFrame(contexts).to_string(index=False),flush=True)
print('PRECISION\n'+read(O/'tables/precision_summary.tsv').to_string(index=False),flush=True)
print('FAMILY_DIAGNOSTICS\n'+pd.DataFrame(families).query('valid_P_rows > 0').to_string(index=False),flush=True)
print('PROPORTION_MIN_P\n'+prop[prop.P.notna()].sort_values('P').head(8)[['layer','dataset','target','denominator','bone_N','lung_N','mean_difference_percentage_points','P','FDR']].to_string(index=False),flush=True)
