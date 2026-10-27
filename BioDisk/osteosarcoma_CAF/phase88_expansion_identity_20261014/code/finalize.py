from pathlib import Path
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='4'
import ast,json,datetime,hashlib,gzip,math
import numpy as np,pandas as pd,anndata as ad,scipy.sparse as sp
from scipy.io import mmread
from threadpoolctl import threadpool_limits
threadpool_limits(4)
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase88_expansion_identity_20261014';Q=O/'preprocess'
def stamp():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def clean(x):
    if isinstance(x,dict):return {k:clean(v) for k,v in x.items()}
    if isinstance(x,list):return [clean(v) for v in x]
    if isinstance(x,float) and not math.isfinite(x):return None
    return x
def dump(p,x):p.write_text(json.dumps(clean(x),ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()
def main():
    manifest=pd.read_csv(Q/'manifest.tsv',sep='\t').fillna('')
    input_summary=pd.read_csv(O/'review/input_summary.tsv',sep='\t').set_index('sample')
    jobs=json.loads((O/'TOOL_RESULTS.json').read_text(encoding='utf-8'));jobmap={(j['method'],j['sample']):j for j in jobs}
    panels=json.loads((Q/'review/lineage_rules.json').read_text(encoding='utf-8'))['panels'];panels['OC_specific']=['OCSTAMP','DCSTAMP','ATP6V0D2','SIGLEC15','SUCNR1','CALCR']
    # Reuse exact inherited evidence and competing/mixed-lineage definitions, without executing the old pipeline.
    src=P/'geo_identity_20260919/code/integrate.py';tree=ast.parse(src.read_text(encoding='utf-8'))
    defs=[node for node in tree.body if isinstance(node,ast.FunctionDef) and node.name in ['evidence','competing','mixed']]
    assert len(defs)==3
    ns={'np':np,'pd':pd,'panels':panels};exec(compile(ast.Module(body=defs,type_ignores=[]),str(src),'exec'),ns)
    evidence,competing,mixed=[ns[k] for k in ['evidence','competing','mixed']]
    dump(O/'resources/IDENTITY_RULE_FUNCTIONS.json',{'source':str(src),'sha256':sha(src),'functions':['evidence','competing','mixed'],'additional_gate':'exclude retained any-method doublet positives','nuclei':'modality validation open; no formal CAF admission'})
    allcells=[];summaries=[];refrows=[];clusterrows=[];checks=[];issues=[]
    for row in manifest.to_dict('records'):
        sm=row['sample'];d=O/'inputs'/sm
        a=ad.read_h5ad(Q/'accepted'/sm/'analysis_ready.h5ad');x=a.layers['counts'].tocsr();syms=a.var.symbol.astype(str).to_numpy()
        o=pd.read_csv(d/'cells.tsv',sep='\t',dtype={'cluster':str}).set_index('cell_id',drop=False);o.index.name=None
        assert o.index.tolist()==a.obs_names.tolist() and o.index.is_unique
        genes=pd.read_csv(d/'genes.tsv',sep='\t');assert genes.gene_id.tolist()==a.var_names.tolist()
        raw,tot=evidence(x,syms);raw.index=o.index;assert np.array_equal(tot,o.n_counts.to_numpy())
        for k in ['Fib','Osteo','Mural','CTA']:assert np.array_equal(raw['det_'+k],o['det_'+k+'_review']), (sm,k)
        for col in raw:o['raw_'+col]=raw[col]
        for method,prefix in [('infercnv','infer'),('copykat','copy')]:
            o[prefix+'_call']='not_assessed';o[prefix+'_score']=np.nan
            j=jobmap[(method,sm)]
            if j['status']=='complete':
                out=O/method/sm
                assert j['exit_code']==0 and (out/'COMPLETE').exists()
                r=pd.read_csv(out/'cell_results.tsv',sep='\t',dtype={'cluster':str}).set_index('cell_id')
                assert r.index.is_unique and set(r.index)==set(o.index) and len(r)==len(o)
                r=r.loc[o.index];assert r.tool_id.tolist()==o.tool_id.tolist()
                o[prefix+'_call']=r['call'].fillna('not_evaluable').to_numpy();o[prefix+'_score']=r.signal_strength.to_numpy()
                if prefix=='infer':o['infer_cutoff']=r.signal_cutoff.to_numpy()
            elif j['status']=='failed':issues.append({'sample':sm,'issue':method+'_failed','detail':'See per-tool log; no normal calls imputed'})
        low=o.infer_call.eq('low_signal')&o.copy_call.eq('diploid')
        high=o.infer_call.eq('elevated_signal')&o.copy_call.eq('aneuploid')
        o['joint']='discordant_or_ambiguous';o.loc[low,'joint']='both_low';o.loc[high,'joint']='both_high'
        missing=o.infer_call.isin(['not_assessed','not_evaluable'])|o.copy_call.isin(['not_assessed','not_evaluable'])
        o.loc[missing,'joint']='not_jointly_evaluable'
        good={}
        for method,normal in [('infer','low_signal'),('copy','diploid')]:
            roles=sorted(o.loc[o.cnv_role.str.startswith('ref_'),'cnv_role'].unique());ok=bool(roles)
            for role in roles+[r.replace('ref_','holdout_') for r in roles]:
                cc=o.loc[o.cnv_role==role,method+'_call'];ev=cc[~cc.isin(['not_evaluable','not_assessed'])]
                frac=float(ev.eq(normal).mean()) if len(ev) else np.nan;passed=len(ev)>=10 and frac>=.9
                refrows.append({'sample':sm,'method':method,'role':role,'total':len(cc),'evaluable':len(ev),'normal_fraction':frac,'pass':passed});ok=ok and passed
            good[method]=bool(ok)
        o['reference_QC_pass']=good['infer'] and good['copy'] and bool(input_summary.at[sm,'cnv_eligible'])
        o['cluster_both_high_fraction']=np.nan;o['cluster_any_high_fraction']=np.nan
        for cl,g in o.groupby('cluster'):
            binary=g.infer_call.isin(['low_signal','elevated_signal'])&g.copy_call.isin(['diploid','aneuploid'])
            evaluable=g.infer_call.isin(['low_signal','elevated_signal','ambiguous_signal'])&g.copy_call.isin(['diploid','aneuploid'])
            both=float(g.loc[evaluable,'joint'].eq('both_high').mean()) if evaluable.any() else np.nan
            anyhi=float((g.loc[binary,'infer_call'].eq('elevated_signal')|g.loc[binary,'copy_call'].eq('aneuploid')).mean()) if binary.any() else np.nan
            o.loc[g.index,'cluster_both_high_fraction']=both;o.loc[g.index,'cluster_any_high_fraction']=anyhi
            clusterrows.append({'sample':sm,'cluster':cl,'lineage':g.lineage_phase3.iloc[0],'cells':len(g),'both_high_fraction':both,'any_high_fraction':anyhi,'binary_evaluable':int(binary.sum())})
        fib=o.lineage_phase3.eq('FibLike');osteo=(raw.det_Osteo>=3)&(raw.pct_Osteo>.2)
        strict=fib&(raw.det_Fib>=2)&low&~mixed(raw)&~o.CTA_context&(raw.det_CTA<2)&~competing(raw)&~osteo&(o.cluster_both_high_fraction.fillna(1)<.2)&(o.cluster_any_high_fraction.fillna(1)<.2)&(o.n_genes>=500)&o.reference_QC_pass&~o.doublet_positive_any.astype(bool)
        o['raw_strict_candidate']=strict;o['ambient_checked']=False;o['corrected_lineage_supported']=False;o['corrected_positive_sensitivity']=False
        corrected=None
        j=jobmap[('decontx',sm)]
        if j['status']=='complete':
            out=O/'decontx'/sm;assert (out/'COMPLETE').exists() and j['exit_code']==0
            corrected=mmread(str(out/'corrected_counts.mtx.gz')).T.tocsr()
            assert corrected.shape==x.shape and np.isfinite(corrected.data).all() and (corrected.data>=0).all()
            diff=corrected-x;assert diff.nnz==0 or diff.data.max()<1e-5
            ct=pd.read_csv(out/'cell_contamination.tsv',sep='\t');assert ct.cell_id.tolist()==o.index.tolist()
            cor,ctot=evidence(corrected,syms);cor.index=o.index
            assert np.allclose(ctot,ct.corrected_UMI) and np.allclose(tot,ct.original_UMI)
            assert ct.contamination.between(0,1).all()
            o['ambient_checked']=True;o['decontx_contamination']=ct.contamination.to_numpy()
            for col in cor:o['corrected_'+col]=cor[col]
            o['corrected_lineage_supported']=(cor.det_Fib>=2)&~competing(cor,osteo=True)&~mixed(cor)&(cor.det_CTA<2)&(ctot>0)
            o['corrected_positive_sensitivity']=(cor.pos_Fib>=2)&~competing(cor,det='pos_',osteo=True)&~mixed(cor,det='pos_')&(cor.pos_CTA<2)&(ctot>0)
        else:issues.append({'sample':sm,'issue':'ambient_not_completed','detail':j['status']})
        o['supported_CAF_candidate']=strict&o.corrected_lineage_supported
        o['supported_positive_sensitivity']=strict&o.corrected_positive_sensitivity
        o['expanded_Fib_expression_candidate']=o.lineage_phase3.isin(['FibLike','MesAmb','Unresolved'])&(raw.det_Fib>=2)
        o['CAF_independently_confirmed']=False;o['technology_validation_open']='single_nucleus' in row['technology']
        o['exclusion_reasons']=''
        masks=[(~fib,'outside_reviewed_FibLike'),(raw.det_Fib<2,'fewer_than2_Fib_genes'),(~low,'not_dual_low'),(~o.reference_QC_pass,'reference_QC_failed_or_unavailable'),(o.n_genes<500,'fewer_than500_genes'),(mixed(raw),'mixed_program'),(competing(raw),'competing_lineage'),(osteo,'osteo_program'),(o.CTA_context|(raw.det_CTA>=2),'CTA_context_or_expression'),(o.cluster_both_high_fraction.fillna(1)>=.2,'cluster_both_high_risk'),(o.cluster_any_high_fraction.fillna(1)>=.2,'cluster_any_high_risk'),(o.doublet_positive_any.astype(bool),'retained_doublet_model_positive'),(~o.ambient_checked,'ambient_not_assessed'),(o.ambient_checked&~o.corrected_lineage_supported,'corrected_lineage_not_supported')]
        for mask,reason in masks:o.loc[mask,'exclusion_reasons']+=reason+';'
        assert o.loc[o.supported_CAF_candidate,'exclusion_reasons'].eq('').all()
        assert (~o.supported_CAF_candidate|(low&o.reference_QC_pass&o.ambient_checked)).all()
        for key in ['dataset','site','patient_id','technology','treatment_raw','disease_stage_raw','clinical_pairing','patient_independence_status']:o[key]=row[key]
        o.to_csv(O/'tables'/(sm+'_identity.tsv.gz'),sep='\t',index=False,encoding='utf-8')
        n=int(o.supported_CAF_candidate.sum())
        if n:
            mask=o.supported_CAF_candidate.to_numpy();b=a[mask].copy();b.obs=o.loc[mask].drop(columns='cell_id').copy()
            for col in b.obs.select_dtypes(include='object'):b.obs[col]=b.obs[col].fillna('').astype(str)
            b.layers['decontx_counts']=corrected[mask].copy();b.uns['identity_status']='operational CAF candidate; no formal baseline admission; source/modality and patient review pending'
            path=O/'objects'/(sm+'_supported_CAF_candidates.h5ad');b.write_h5ad(path,compression='gzip')
            test=ad.read_h5ad(path);assert test.obs_names.equals(b.obs_names) and (test.layers['counts']!=x[mask]).nnz==0 and (test.layers['decontx_counts']!=corrected[mask]).nnz==0
            del b,test
        summary={'sample':sm,'dataset':row['dataset'],'patient_id':row['patient_id'],'site':row['site'],'technology':row['technology'],'cells':len(o),'reviewed_FibLike':int(fib.sum()),'expanded_Fib_expression_candidates':int(o.expanded_Fib_expression_candidate.sum()),'dual_low_FibLike':int((fib&low).sum()),'raw_strict_candidates':int(strict.sum()),'supported_CAF_candidates':n,'positive_count_sensitivity':int(o.supported_positive_sensitivity.sum()),'reference_QC_pass':bool(o.reference_QC_pass.iloc[0]),'ambient_checked':corrected is not None,'CNV_input_eligible':bool(input_summary.at[sm,'cnv_eligible'])}
        summaries.append(summary);allcells.append(o);checks.append({'sample':sm,'cell_order_pass':True,'raw_panel_reproduction_pass':True,'CNV_mapping_pass':True,'ambient_value_checks_pass':corrected is not None,'supported_object_readback_pass':True,'supported_cells':n})
        if not summary['reference_QC_pass']:issues.append({'sample':sm,'issue':'reference_calibration_not_passed','detail':'No supported CAF admission from this sample under frozen rules'})
        print('IDENTITY',json.dumps(summary),flush=True);del a,x,corrected
    cells=pd.concat(allcells,ignore_index=True);assert cells.cell_id.is_unique
    assert len(cells)==sum(input_summary.cells)
    cells.to_csv(O/'tables/all_cell_identity_evidence.tsv.gz',sep='\t',index=False,encoding='utf-8')
    cells.loc[cells.supported_CAF_candidate,['cell_id','sample','patient_id','dataset','site','technology','technology_validation_open']].to_csv(O/'tables/supported_CAF_candidate_cells.tsv',sep='\t',index=False,encoding='utf-8')
    pd.DataFrame(summaries).to_csv(O/'tables/sample_identity_summary.tsv',sep='\t',index=False,encoding='utf-8')
    pd.DataFrame(refrows).to_csv(O/'tables/reference_calibration.tsv',sep='\t',index=False,encoding='utf-8')
    pd.DataFrame(clusterrows).to_csv(O/'tables/cluster_CNV_risk.tsv',sep='\t',index=False,encoding='utf-8')
    for row in manifest.to_dict('records'):
        issues.append({'sample':row['sample'],'issue':'clinical_or_patient_independence_review_pending','detail':'NCC retains within-source Patient4 pre/post pairing; cross-study independence unresolved. CELLxGENE limb is not Lung; treatment/capture mapping pending.'})
    pd.DataFrame(issues).to_csv(O/'tables/unresolved_issues.tsv',sep='\t',index=False,encoding='utf-8')
    rostercols=['sample','patient_id','dataset','site','technology','treatment_raw','disease_stage_raw','clinical_pairing','patient_independence_status','source_evidence','source_library','source_Cellgen_ID','author_patient_status','author_histological_subtype']
    manifest[[c for c in rostercols if c in manifest]].to_csv(O/'tables/patient_sample_roster.tsv',sep='\t',index=False,encoding='utf-8')
    hashes=json.loads((O/'resources/SOURCE_HASHES.json').read_text(encoding='utf-8'));sourcechecks=[]
    for r in hashes:sourcechecks.append({'path':r['path'],'unchanged':sha(Path(r['path']))==r['sha256']})
    assert all(r['unchanged'] for r in sourcechecks),'Input changed during run'
    baselinechecks=[]
    for r in json.loads((O/'resources/BASELINE_HASHES.json').read_text(encoding='utf-8')):baselinechecks.append({'path':r['source'],'unchanged':sha(Path(r['source']))==r['sha256']})
    dump(O/'resources/FINAL_INTEGRITY.json',{'sample_checks':checks,'source_checks':sourcechecks,'baseline_state_checks':baselinechecks})
    failed=[j for j in jobs if j['status']=='failed'];skipped=[j for j in jobs if j['status']=='skipped_input_gate']
    a={'status':'COMPUTATION_COMPLETE_WITH_TOOL_FAILURES_REVIEW_REQUIRED' if failed else 'COMPUTATION_COMPLETE_MANUAL_REVIEW_REQUIRED','completed_utc':stamp(),'input_libraries':12,'accepted_cells':len(cells),'reviewed_FibLike':int(cells.lineage_phase3.eq('FibLike').sum()),'supported_CAF_candidate_cells':int(cells.supported_CAF_candidate.sum()),'samples_with_supported_candidates':int((pd.DataFrame(summaries).supported_CAF_candidates>0).sum()),'formal_new_CAF_patients':0,'successful_tool_jobs':sum(j['status']=='complete' for j in jobs),'failed_tool_jobs':len(failed),'skipped_tool_jobs':len(skipped),'cell_accounting_pass':True,'input_immutability_pass':True,'output_readback_checks_pass':True,'baseline_state_unchanged':all(r['unchanged'] for r in baselinechecks),'manual_review_pending':True,'visual_QA_claimed':False,'new_inferential_tests':False,'no_baseline_admission':True}
    dump(O/'ACCEPTANCE.json',a)
    lines=['CAF 扩充数据后台身份复核：自动计算验收摘要',f'完成时间 UTC：{a["completed_utc"]}',f'状态：{a["status"]}','',f'共 12 份文库：NCC 9 份样本（8 个患者标签，BS4/BS14 同患者治疗前后配对）及 CELLxGENE 3 份文库（2 个供者；EXT_OS_1 的两个文库属于同一供者）。',f'本轮处理后细胞/细胞核总数：{len(cells):,}。',f'描述性谱系复核后 FibLike：{a["reviewed_FibLike"]:,}；同时通过双方法 CNV、参考细胞、原始表达和 decontX 支持门槛的操作性 CAF 候选：{a["supported_CAF_candidate_cells"]:,}。',f'工具任务成功 {a["successful_tool_jobs"]}，失败 {len(failed)}，因输入门槛跳过 {len(skipped)}。','跳过、校准失败及方法不一致均保留在表中，不按正常细胞填补。','原始计数与标识符、CNV 结果映射、校正矩阵及输出对象回读检查均已执行。','', '本轮只完成身份证据计算，人工验收尚未完成；未执行患者级差异检验，未修改正式 CAF 基线，正式新增患者数为 0。','CELLxGENE 为 10x multiome 单核数据，供者取材字段为 limb；未将其归为肺组织。','NCC 与其他研究、CELLxGENE 与其他研究的患者独立性仍需结合临床元数据核实。NCI/S23 既有患者映射与取材冲突未在本轮静默改写。','', '明天重点查看：','1. tables/sample_identity_summary.tsv：各样本候选、排除及方法完成情况。','2. tables/reference_calibration.tsv 与 cluster_CNV_risk.tsv：正常参考留出校准和肿瘤相关风险。','3. tables/*_identity.tsv.gz：逐细胞原始/校正表达、双方法 CNV 及排除原因。','4. tables/patient_sample_roster.tsv 与 unresolved_issues.tsv：配对、取材和未决事项。','5. objects/*_supported_CAF_candidates.h5ad：符合保守操作性门槛的候选对象（若无通过细胞，则不生成空的候选对象）。','6. ACCEPTANCE.json、JOB_FINISHED.json、logs/*：自动验收、退出状态及日志。']
    (O/'delivery/OVERNIGHT_SUMMARY_ZH.txt').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print(json.dumps(a),flush=True)
if __name__=='__main__':main()
