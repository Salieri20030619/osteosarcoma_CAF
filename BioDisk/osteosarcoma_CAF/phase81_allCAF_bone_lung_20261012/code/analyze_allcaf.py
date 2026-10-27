from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='2'
import json,hashlib,importlib.util,traceback
from datetime import datetime,timezone
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF')
O=P/'phase81_allCAF_bone_lung_20261012'
S=P/'phase73_F134_bone_lung_20261009/strict_v2'
A=P/'phase41_full_transcriptome_setup_20261001'
for d in ['code','tables','objects','figures','resources','QA','delivery','logs']:(O/d).mkdir(parents=True,exist_ok=True)
def read(f,**kw):return pd.read_csv(f,sep='\t',low_memory=False,**kw)
def save(d,f):d.to_csv(f,sep='\t',index=False,na_rep='NA')
def dump(x,f):f.write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str)+'\n',encoding='utf-8')
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
def now():return datetime.now(timezone.utc).isoformat()
def truth(s):return s.astype(str).str.lower().isin(['true','1','1.0'])
spec=importlib.util.spec_from_file_location('audited_stats',S/'code/analyze_bone_lung_v2.py');core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core);core.SEED=20261012
def stats(x,y,paired):
    z=core.stats(x,y,paired)
    if paired:
        both=np.isfinite(x)&np.isfinite(y);d=np.where(both,y-x,np.nan)
        z['mean_lung_minus_bone']=core.mean(d)
        with np.errstate(invalid='ignore'):z['median_lung_minus_bone']=np.nanmedian(d,axis=0)
    return z
def adjust(df):
    df=df.copy();df['FDR']=np.nan;df['legacy_NA_as_one_FDR_diagnostic']=np.nan;df['FDR_family_valid_tests']=0;df['family_total_rows']=0
    for fam,g in df.groupby('family',sort=False):
        idx=g.index[g.P.notna()];df.loc[g.index,'FDR_family_valid_tests']=len(idx);df.loc[g.index,'family_total_rows']=len(g)
        if len(idx):
            df.loc[idx,'FDR']=multipletests(df.loc[idx,'P'],method='fdr_bh')[1]
            full=multipletests(g.P.fillna(1),method='fdr_bh')[1]
            df.loc[g.index,'legacy_NA_as_one_FDR_diagnostic']=np.where(g.P.notna(),full,np.nan)
    df['nominal']=df.P.lt(.05);df['FDR_significant']=df.FDR.lt(.05);return df

def main():
    assert not (O/'COMPUTATION_COMPLETE.json').exists()
    protocol=dict(frozen_at=now(),user_request='先用所有CAF细胞中的两转移灶做比对',
        interpretation='Under previously approved scheme1: all CAF candidates pooled across F classes, anatomical Bone vs Lung regardless of metastatic history; not two selected individual samples',
        supported_layer=1224,expanded_layer=8174,layers_separate_not_independent=True,
        labels='Original candidate IDs, classes and identity confidence retained; all F1-F5 pooled within each layer for testing',
        candidate_identity='Expanded pool includes unresolved/osteogenic expression; not equivalent to confirmed pure fibroblasts',
        clinical_sites='Phase73 strict_v2 accepted Bone/Lung anatomical adjudication; other or mixed sites separately retained',
        replication='Source-specific patient-site pseudobulk; paired donors removed from independent contrast of the same endpoint',
        gene_expression='Raw sum per patient/site; full-library log2(CP10k+0.1)',
        program_expression='Mean per-cell log1p(CP10k), average measured frozen panel genes, patient/site mean',
        testing='Audited exact patient label permutation/sign-flip; predefined MC65536+observed when exact assignments exceed65536 (paired exact if N<=15); seed20261012',
        insufficient_replication='Less than2 per independent arm or less than2 paired donors: effect only, no P',
        nominal_selection='Target raw P<0.05, retain up/down, do not require FDR',
        FDR='BH among valid P values within layer/source/design/sensitivity/hypothesis type; each source and paired design separate',
        FDR_change='Invalid tests remain NA and do not count as P=1 in BH. This family definition is frozen before analysis. Legacy NA-as-one q retained as diagnostic.',
        proportion='All CAF candidates / full provider-called cells; zero numerators retained, zero denominator undefined',
        sensitivity='Exclude26 repaired model-positive cells where present; keep original labels and design',
        ORA='Nominal gene-list up/down, canonical valid-gene background, per-library/direction BH; modeled enrichment not patient-effect P',
        no_new_cell_entry_cutoffs=True,no_surface_or_flow_filter=True)
    dump(protocol,O/'resources/PROTOCOL_FROZEN.json')
    paths=dict(clinical=S/'resources/clinical_sample_manifest.tsv',supported=S/'tables/supported_fixed_cells_with_sites.tsv.gz',expanded=S/'tables/expanded_fixed_cells_with_sites.tsv.gz',counts=A/'objects/all_candidate_full_genes.h5ad',presence=A/'objects/sample_feature_presence.npz',canonical=P/'phase58_F1F5_round1_20261004/resources/source_canonical_feature_map.tsv',programs=P/'phase58_F1F5_round1_20261004/resources/programs_frozen.json',libraries=P/'phase58_F1F5_round1_20261004/resources/libraries_canonical.json',stats_code=S/'code/analyze_bone_lung_v2.py')
    before={k:dict(path=str(v),sha256=sha(v),bytes=v.stat().st_size) for k,v in paths.items()}
    dump(before,O/'resources/INPUT_MANIFEST.json')
    cm=read(paths['clinical'],keep_default_na=False).set_index('sample');save(cm.reset_index(),O/'tables/clinical_sample_manifest.tsv')
    layers={k:read(paths[k],keep_default_na=False).set_index('cell_id',drop=False) for k in ['supported','expanded']}
    assert len(layers['supported'])==1224 and len(layers['expanded'])==8174
    assert set(layers['supported'].index)<=set(layers['expanded'].index)
    for layer,ob in layers.items():
        assert (ob.strict_site==ob['sample'].map(cm.strict_site_group)).all()
        save(ob[['cell_id','sample','dataset','donor_record','strict_site','analysis_class','review_lineage_code','CNV_evidence_current','doublet_positive_any']],O/'tables'/f'{layer}_allCAF_membership.tsv.gz')
    a=ad.read_h5ad(paths['counts']);assert a.shape==(8174,63398)
    assert np.allclose(np.asarray(a.X.sum(1)).ravel(),layers['expanded'].loc[a.obs_names,'full_library_counts'].astype(float))
    features=a.var_names.astype(str).to_numpy();symbols=a.var.symbol.astype(str).to_numpy()
    presence=sparse.load_npz(paths['presence']).tocsr();sidx=read(A/'tables/feature_presence_sample_index.tsv').set_index('sample').presence_row.to_dict()
    canonical=read(paths['canonical'],keep_default_na=False)
    panels=json.loads(paths['programs'].read_text())['sets'];names=list(panels);libs=json.loads(paths['libraries'].read_text())
    logx=a.X.astype(float).tocsr();logx=sparse.diags(1e4/np.asarray(a.X.sum(1)).ravel())@logx;logx.data=np.log1p(logx.data)
    source_mats={};coverage=[]
    for ds,sg in cm.groupby('dataset'):
        gi=np.flatnonzero(np.asarray(presence[[sidx[s] for s in sg.index]].sum(0)).ravel()>0)
        pos={features[j]:i for i,j in enumerate(gi)};can=canonical[canonical.dataset.eq(ds)];cp={r.symbol:pos[r.feature] for r in can.itertuples() if r.feature in pos}
        rr=[];cc=[]
        for i,name in enumerate(names):
            m=[cp[s] for s in panels[name] if s in cp];rr+=m;cc += [i]*len(m)
            coverage.append(dict(dataset=ds,program=name,defined=len(panels[name]),measured=len(m)))
        source_mats[ds]=(gi,sparse.csr_matrix((np.ones(len(rr)),(rr,cc)),shape=(len(gi),len(names))))
    save(pd.DataFrame(coverage),O/'tables/program_coverage.tsv')
    props=[];prop_values=[];contexts=[];roster=[];all_genes=[];all_programs=[];composition=[]
    for layer,original in layers.items():
      for sensitivity in ['all_cells','exclude_doublets']:
        ob=original if sensitivity=='all_cells' else original[~truth(original.doublet_positive_any)]
        for ds,ms in cm[cm.strict_site_group.isin(['Bone','Lung'])].groupby('dataset'):
            frac={'Bone':{},'Lung':{}};expr={'Bone':{},'Lung':{}};program={'Bone':{},'Lung':{}}
            gi,pm=source_mats[ds];fe=features[gi];sy=symbols[gi]
            for (donor,site),sg in ms.groupby(['donor_record','strict_site_group']):
                q=ob[ob['sample'].isin(sg.index)];den=int(sg.full_called_cells.astype(int).sum())
                if den:frac[site][donor]=len(q)/den
                prop_values.append(dict(layer=layer,sensitivity=sensitivity,dataset=ds,donor_record=donor,site=site,cells=len(q),full_called_cells=den,fraction=len(q)/den if den else np.nan))
                if not len(q):continue
                ix=a.obs_names.get_indexer(q.index);assert (ix>=0).all()
                raw=np.asarray(a.X[ix][:,gi].sum(0)).ravel();library=float(q.full_library_counts.astype(float).sum())
                ss=q['sample'].unique();valid=np.asarray(presence[[sidx[s] for s in ss]][:,gi].sum(0)).ravel()==len(ss)
                val=np.log2(raw/library*1e4+.1);val[~valid]=np.nan;expr[site][donor]=val
                lm=np.asarray(logx[ix][:,gi].mean(0)).ravel();nm=np.asarray(valid.astype(float)@pm).ravel();nv=np.asarray(np.where(valid,lm,0)@pm).ravel();pv=np.full(len(names),np.nan);np.divide(nv,nm,out=pv,where=nm>0);program[site][donor]=pv
                roster.append(dict(layer=layer,sensitivity=sensitivity,dataset=ds,donor_record=donor,site=site,cells=len(q),library_UMIs=library,samples=';'.join(ss)))
                for cls,n in q.analysis_class.value_counts().items():composition.append(dict(layer=layer,sensitivity=sensitivity,dataset=ds,donor_record=donor,site=site,class_id=cls,cells=n,fraction=n/len(q)))
            for endpoint,values in [('proportion',frac),('expression',expr)]:
                paired=sorted(set(values['Bone'])&set(values['Lung']))
                for design in ['paired','independent']:
                    ids1=paired if design=='paired' else sorted(set(values['Bone'])-set(paired));ids2=paired if design=='paired' else sorted(set(values['Lung'])-set(paired))
                    key='__'.join([sensitivity,layer,ds,'All_CAF',design]);meta=dict(sensitivity=sensitivity,layer=layer,dataset=ds,target='All_CAF',design=design,context_id=key)
                    if endpoint=='proportion':
                        x=np.array([values['Bone'][d] for d in ids1]).reshape(-1,1);y=np.array([values['Lung'][d] for d in ids2]).reshape(-1,1)
                        st=stats(x,y,design=='paired');st['bone_mean_fraction']=core.mean(x);st['lung_mean_fraction']=core.mean(y)
                        for k,v in meta.items():st[k]=v
                        st['family']=key+'|proportion';props.append(st)
                        continue
                    contexts.append(dict(**meta,bone_N=len(ids1),lung_N=len(ids2),bone_patients=';'.join(ids1),lung_patients=';'.join(ids2),status='testable' if min(len(ids1),len(ids2))>=2 else ('direction_only' if ids1 and ids2 else 'group_absent')))
                    if not ids1 or not ids2:continue
                    x=np.vstack([values['Bone'][d] for d in ids1]);y=np.vstack([values['Lung'][d] for d in ids2]);bx=np.vstack([program['Bone'][d] for d in ids1]);ly=np.vstack([program['Lung'][d] for d in ids2])
                    gd=stats(x,y,design=='paired');gd['feature']=fe;gd['symbol']=sy
                    pg=stats(bx,ly,design=='paired');pg['program']=names
                    for df,kind in [(gd,'genes'),(pg,'programs')]:
                        for k,v in meta.items():df[k]=v
                        df['family']=key+'|'+kind
                    all_genes.append(gd);all_programs.append(pg)
                    np.savez_compressed(O/'objects'/f'{key}.npz',features=fe.astype(str),symbols=sy.astype(str),bone_patients=np.array(ids1,dtype=str),lung_patients=np.array(ids2,dtype=str),bone_expression=x,lung_expression=y,bone_programs=bx,lung_programs=ly,programs=np.array(names,dtype=str))
                    print('COMPARED',key,len(ids1),len(ids2),'valid',int(gd.P.notna().sum()),flush=True)
    genes=adjust(pd.concat(all_genes,ignore_index=True));programs=adjust(pd.concat(all_programs,ignore_index=True));pr=adjust(pd.concat(props,ignore_index=True))
    for df,name in [(genes,'genes'),(programs,'programs'),(pr,'proportions')]:
        save(df,O/'tables'/f'all_{name}_statistics.tsv.gz');save(df[df.nominal],O/'tables'/f'nominal_{name}.tsv.gz')
    save(pd.DataFrame(contexts),O/'tables/comparison_registry.tsv');save(pd.DataFrame(roster),O/'tables/patient_expression_roster.tsv');save(pd.DataFrame(prop_values),O/'tables/patient_CAF_fractions.tsv');save(pd.DataFrame(composition),O/'tables/patient_Fclass_composition.tsv')
    summary=[]
    for key,g in genes.groupby('context_id'):
        row=next(r for r in contexts if r['context_id']==key);pg=programs[programs.context_id.eq(key)]
        summary.append(dict(**row,valid_genes=int(g.P.notna().sum()),nominal_up=int((g.nominal&g.mean_lung_minus_bone.gt(0)).sum()),nominal_down=int((g.nominal&g.mean_lung_minus_bone.lt(0)).sum()),FDR_genes=int(g.FDR_significant.sum()),nominal_programs=int(pg.nominal.sum()),FDR_programs=int(pg.FDR_significant.sum()),legacy_diagnostic_gene_FDR_hits=int(g.legacy_NA_as_one_FDR_diagnostic.lt(.05).sum()),min_gene_P=g.P.min(),min_gene_FDR=g.FDR.min()))
    save(pd.DataFrame(summary),O/'tables/comparison_results_summary.tsv')
    ora=[]
    for key,g in genes[genes.sensitivity.eq('all_cells')].groupby('context_id'):
        can=canonical[canonical.dataset.eq(g.dataset.iloc[0])][['feature','symbol']].drop_duplicates('feature');v=g.drop(columns='symbol').merge(can,on='feature',validate='one_to_one');v=v[v.P.notna()];background=set(v.symbol)
        for direction in ['up','down']:
            sel=set(v.loc[v.P.lt(.05)&(v.mean_lung_minus_bone.gt(0) if direction=='up' else v.mean_lung_minus_bone.lt(0)),'symbol'])
            if not sel:continue
            for library,terms in libs.items():
                for term,members in terms.items():
                    measured=set(members)&background
                    if not measured or len(measured)==len(background):continue
                    overlap=sel&measured
                    ora.append(dict(context_id=key,direction=direction,library=library,term=term,background_genes=len(background),selected_genes=len(sel),term_genes=len(measured),overlap_genes=len(overlap),overlap_symbols=';'.join(sorted(overlap)),P=float(hypergeom.sf(len(overlap)-1,len(background),len(measured),len(sel))),family=key+'|ORA|'+direction+'|'+library))
    if ora:
        od=adjust(pd.DataFrame(ora));save(od,O/'tables/all_ORA_statistics.tsv.gz');save(od[od.nominal],O/'tables/nominal_ORA.tsv.gz')
    # Every nominal signal remains; add paired and doublet evidence instead of filtering.
    for df,kind,item in [(genes,'genes','feature'),(programs,'programs','program')]:
        main=df[df.sensitivity.eq('all_cells')&df.nominal].copy();sens=df[df.sensitivity.eq('exclude_doublets')][['layer','dataset','design',item,'mean_lung_minus_bone','P','FDR']].rename(columns={'mean_lung_minus_bone':'doublet_excluded_effect','P':'doublet_excluded_P','FDR':'doublet_excluded_FDR'})
        main=main.merge(sens,on=['layer','dataset','design',item],how='left',validate='one_to_one')
        paired=df[df.sensitivity.eq('all_cells')&df.design.eq('paired')][['layer','dataset',item,'mean_lung_minus_bone','bone_N','P','FDR']].rename(columns={'mean_lung_minus_bone':'paired_effect','bone_N':'paired_N','P':'paired_P','FDR':'paired_FDR'})
        main=main.merge(paired,on=['layer','dataset',item],how='left',validate='many_to_one')
        main['paired_direction_agrees']=np.where(main.paired_effect.notna(),main.mean_lung_minus_bone*main.paired_effect>0,np.nan)
        save(main,O/'tables'/f'main_nominal_{kind}_with_sensitivity.tsv.gz')
    checks=dict(status='computation_complete_QA_pending',completed_at=now(),layers={k:len(v) for k,v in layers.items()},raw_library_check=True,original_inputs_unchanged=all(sha(paths[k])==v['sha256'] for k,v in before.items()),main_testable_contexts=sum(r['status']=='testable' and r['sensitivity']=='all_cells' for r in contexts),main_nominal_gene_records=int((genes.sensitivity.eq('all_cells')&genes.nominal).sum()),main_FDR_gene_records=int((genes.sensitivity.eq('all_cells')&genes.FDR_significant).sum()),main_nominal_program_records=int((programs.sensitivity.eq('all_cells')&programs.nominal).sum()),main_FDR_program_records=int((programs.sensitivity.eq('all_cells')&programs.FDR_significant).sum()),main_nominal_proportions=int((pr.sensitivity.eq('all_cells')&pr.nominal).sum()),new_analysis=True,FDR_family_rule='valid P only; legacy q diagnostic retained')
    dump(checks,O/'COMPUTATION_COMPLETE.json');print(json.dumps(checks,ensure_ascii=False),flush=True)
    print(pd.DataFrame(summary).query("sensitivity=='all_cells'")[['layer','dataset','design','bone_N','lung_N','valid_genes','nominal_up','nominal_down','FDR_genes','nominal_programs','FDR_programs']].to_string(index=False),flush=True)

if __name__=='__main__':
    try:main()
    except BaseException:
        (O/'logs/analysis.exit_code').write_text('1\n');traceback.print_exc();raise
    else:(O/'logs/analysis.exit_code').write_text('0\n')
