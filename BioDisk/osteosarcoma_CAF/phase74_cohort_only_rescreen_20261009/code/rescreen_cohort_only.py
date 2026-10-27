from pathlib import Path
import json,hashlib,datetime,traceback,zipfile,shutil
import numpy as np,pandas as pd

P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF')
S=P/'phase73_F134_bone_lung_20261009/strict_v2'
O=P/'phase74_cohort_only_rescreen_20261009'
def read(f,**kw):return pd.read_csv(f,sep='\t',low_memory=False,**kw)
def save(d,name):d.to_csv(O/'tables'/name,sep='\t',index=False,na_rep='NA')
def dump(f,x):f.write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str)+'\n',encoding='utf-8')
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def sha(f):
 h=hashlib.sha256()
 with f.open('rb') as s:
  for b in iter(lambda:s.read(8388608),b''):h.update(b)
 return h.hexdigest()
def mean(a):
 n=np.isfinite(a).sum(0);o=np.full(a.shape[1],np.nan);np.divide(np.nansum(a,axis=0),n,out=o,where=n>0);return o

def enrich_evidence(selected,full,kind):
 if selected.empty:return selected
 feature='feature' if kind=='gene' else 'program';npkey='features' if kind=='gene' else 'programs';arraykey='expression' if kind=='gene' else 'programs'
 frames=[]
 for key,g in selected.groupby('context_id',sort=False):
  g=g.copy();z=np.load(S/'objects'/(key+'.npz'));idx={s:i for i,s in enumerate(z[npkey].astype(str))};ix=np.array([idx[s] for s in g[feature]])
  x=z['bone_'+arraykey][:,ix];y=z['lung_'+arraykey][:,ix];paired=g.design.iloc[0]=='paired';values=[]
  if paired:
   delta=y-x
   for i in range(len(delta)):
    if len(delta)>1:values.append(mean(np.delete(delta,i,axis=0)))
  else:
   for i in range(len(x)):
    if len(x)>1:values.append(mean(y)-mean(np.delete(x,i,axis=0)))
   for i in range(len(y)):
    if len(y)>1:values.append(mean(np.delete(y,i,axis=0))-mean(x))
  vals=np.array(values);g['LOPO_effect_min']=np.nanmin(vals,axis=0);g['LOPO_effect_max']=np.nanmax(vals,axis=0)
  g['LOPO_direction_consistent']=np.all(np.sign(vals)==np.sign(g.mean_lung_minus_bone.to_numpy())[None,:],axis=0)
  g['LOPO_role']='descriptive effect sensitivity, no new P values'
  skey=key.replace('all_cells__','exclude_doublets__',1)
  other=full[full.context_id.eq(skey)][[feature,'mean_lung_minus_bone','P','FDR']].rename(columns={'mean_lung_minus_bone':'doublet_excluded_effect','P':'doublet_excluded_P','FDR':'doublet_excluded_FDR'})
  g=g.merge(other,on=feature,how='left',validate='one_to_one')
  g['doublet_direction_consistent']=np.sign(g.mean_lung_minus_bone)==np.sign(g.doublet_excluded_effect)
  g['doublet_nominal_retained']=g.doublet_excluded_P.lt(.05)
  pkey=key.rsplit('__',1)[0]+'__paired'
  if not paired and (S/'objects'/(pkey+'.npz')).exists():
   pz=np.load(S/'objects'/(pkey+'.npz'));pi={s:i for i,s in enumerate(pz[npkey].astype(str))};px=pz['bone_'+arraykey][:,[pi[s] for s in g[feature]]];py=pz['lung_'+arraykey][:,[pi[s] for s in g[feature]]];dd=py-px
   g['paired_followup_N']=np.isfinite(dd).sum(0);g['paired_followup_mean_effect']=mean(dd);sgn=np.sign(g.mean_lung_minus_bone.to_numpy())
   g['paired_followup_same_direction_patients']=np.sum(np.sign(dd)==sgn[None,:],axis=0)
   g['paired_followup_opposite_direction_patients']=np.sum(np.sign(dd)==-sgn[None,:],axis=0)
  else:
   g['paired_followup_N']=0;g['paired_followup_mean_effect']=np.nan;g['paired_followup_same_direction_patients']=0;g['paired_followup_opposite_direction_patients']=0
  frames.append(g)
 return pd.concat(frames,ignore_index=True)

def main():
 assert not (O/'COMPLETE.json').exists()
 for d in ['code','resources','tables','QA','delivery']:(O/d).mkdir(parents=True,exist_ok=True)
 gene_files=sorted((S/'tables').glob('*__all_gene_statistics.tsv.gz'))
 inputs=gene_files+[S/'tables'/n for n in ['all_program_statistics.tsv.gz','all_proportion_statistics.tsv','all_ORA_statistics.tsv.gz','expression_comparison_registry.tsv','expression_patient_roster.tsv','proportion_test_patient_roster.tsv','program_coverage.tsv','frozen_target_registry.tsv']]
 assert all(not any(s in f.name.lower() for s in ['surface','marker','flow']) for f in inputs)
 dump(O/'resources/INPUT_MANIFEST.json',[dict(path=str(f),sha256=sha(f),bytes=f.stat().st_size) for f in inputs])
 protocol=dict(frozen_at=now(),scope='Cohort/site differences only; no surface-marker nomination or flow decisions',source='Phase73 accepted strict_v2 full statistics',labels_and_grouping='Unchanged original supported/expanded classes and fixed internal states; strict Bone/Lung; paired and independent separated within source',selection='Valid target raw P<0.05, both signs; keep original effect/N/P/FDR/family',proportions='Separate from expression',genes='All measured gene features, without membrane, antibody, marker or gene-location filtering',programs='All1867 frozen programs; same raw P gate; no functional-theme gate',ORA='Gene-list overrepresentation results kept separately as annotation-model statistics, not independent-patient P',FDR='Original full-family FDR retained, not recomputed on selected lists',sensitivity='LOPO directions for every nominal gene/program; repaired-doublet results; paired effects; other-source contrasts with N',interpretation='Main class and nested state findings not independent replications; layers overlap; assess evidence at population level before any flow step',no_new_inferential_test=True,no_new_clustering_or_gene_sets=True,no_FDR_gate=True,no_surface_inputs=True)
 dump(O/'resources/RESCREEN_PROTOCOL_FROZEN.json',protocol)
 # Verify upstream ordering: full gene statistics are written before the downstream surface intersection.
 code=(S/'code/analyze_bone_lung_v2.py').read_text();gene_write=code.index("'__all_gene_statistics.tsv.gz'") if "'__all_gene_statistics.tsv.gz'" in code else code.index('__all_gene_statistics.tsv.gz')
 surface_intersection=code.index('canon.symbol.isin(surface_symbols)')
 assert gene_write<surface_intersection
 audit={'upstream_full_gene_tests_precede_surface_intersection':True,'original_numerical_gene_screen_not_surface_prefiltered':True,'previous_reporting_mixed_clinical_and_detection_priorities':True,'current_selected_inputs_contain_no_marker_surface_flow_tables':True,'original_statistics_preserved':True,'code_sha256':sha(S/'code/analyze_bone_lung_v2.py')}
 dump(O/'QA/SCREEN_CHAIN_AUDIT.json',audit)
 context_summary=[];gene_selected=[];all_base_genes=[]
 for f in gene_files:
  d=read(f);base=d[d.sensitivity.eq('all_cells')].copy();all_base_genes.append(base)
  good=base[base.P.notna()&base.P.lt(.05)&base.mean_lung_minus_bone.notna()].copy()
  assert (good.bone_N>=2).all() and (good.lung_N>=2).all()
  gene_selected.append(enrich_evidence(good,d,'gene'))
  for key,g in base.groupby('context_id'):
   context_summary.append(dict(context_id=key,layer=g.layer.iloc[0],dataset=g.dataset.iloc[0],target=g.target.iloc[0],parent=g.parent.iloc[0],level=g.level.iloc[0],design=g.design.iloc[0],bone_N=int(g.bone_N.max()),lung_N=int(g.lung_N.max()),all_feature_rows=len(g),valid_P_features=int(g.P.notna().sum()),nominal_up=int((g.P.lt(.05)&g.mean_lung_minus_bone.gt(0)).sum()),nominal_down=int((g.P.lt(.05)&g.mean_lung_minus_bone.lt(0)).sum()),FDR_significant=int(g.FDR.lt(.05).sum())))
  print('RESCREENED',f.name,len(good),flush=True)
 genes=pd.concat([g for g in gene_selected if not g.empty],ignore_index=True)
 # Other-source comparisons: no pooling and no exclusion just because the other source P is unavailable.
 allbase=pd.concat(all_base_genes,ignore_index=True)
 keys=['layer','target','design','feature']
 primary=genes[keys+['context_id','dataset','symbol','mean_lung_minus_bone']].rename(columns={'context_id':'primary_context','dataset':'primary_dataset','mean_lung_minus_bone':'primary_effect'})
 other=allbase[keys+['context_id','dataset','bone_N','lung_N','mean_lung_minus_bone','P','FDR']].rename(columns={'context_id':'other_context','dataset':'other_dataset','mean_lung_minus_bone':'other_effect','P':'other_P','FDR':'other_FDR'})
 cross=primary.merge(other,on=keys,how='inner');cross=cross[cross.primary_dataset.ne(cross.other_dataset)].copy()
 cross['same_direction']=np.sign(cross.other_effect)==np.sign(cross.primary_effect)
 cross['evidence_role']='other-source direction, not pooled inference'
 save(cross,'gene_other_source_evidence.tsv.gz')
 save(genes,'cohort_nominal_genes_ALL.tsv.gz')
 for target,g in genes.groupby(['layer','target']):save(g,'__'.join(target)+'__cohort_nominal_genes.tsv')
 gp=read(S/'tables/all_program_statistics.tsv.gz');base=gp[gp.sensitivity.eq('all_cells')];prog=base[base.P.notna()&base.P.lt(.05)&base.mean_lung_minus_bone.notna()]
 prog=enrich_evidence(prog,gp,'program');save(prog,'cohort_nominal_programs_ALL.tsv')
 save(base,'cohort_all_program_statistics.tsv.gz')
 proportions=read(S/'tables/all_proportion_statistics.tsv');proportions=proportions[proportions.sensitivity.eq('all_cells')];save(proportions,'cohort_all_proportion_statistics.tsv');save(proportions[proportions.P.lt(.05)],'cohort_nominal_proportions.tsv')
 ora=read(S/'tables/all_ORA_statistics.tsv.gz');ora=ora[ora.sensitivity.eq('all_cells')];save(ora,'cohort_all_ORA_statistics.tsv.gz');save(ora[ora.P.lt(.05)],'cohort_nominal_ORA.tsv')
 registry=read(S/'tables/expression_comparison_registry.tsv');registry=registry[registry.sensitivity.eq('all_cells')];save(registry,'cohort_comparison_registry.tsv')
 roster=read(S/'tables/expression_patient_roster.tsv');roster=roster[roster.sensitivity.eq('all_cells')];save(roster,'cohort_expression_patient_roster.tsv')
 summaries=[]
 for r in context_summary:
  g=genes[genes.context_id.eq(r['context_id'])];ps=prog[prog.context_id.eq(r['context_id'])]
  ctx=registry[registry.context_id.eq(r['context_id'])].iloc[0]
  for side in ['bone','lung']:
   donors=str(ctx[side+'_patients']).split(';');rr=roster[roster.layer.eq(r['layer'])&roster.target.eq(r['target'])&roster.dataset.eq(r['dataset'])&roster.site.eq(side.title())&roster.donor_record.isin(donors)]
   r[side+'_cells_by_patient']=';'.join(f'{x.donor_record}:{x.cells}' for x in rr.itertuples())
   r[side+'_total_cells']=int(rr.cells.sum())
  r.update(nominal_programs=len(ps),nominal_gene_LOPO_same_direction=int(g.LOPO_direction_consistent.sum()) if len(g) else 0,nominal_gene_doublet_nominal_retained=int(g.doublet_nominal_retained.sum()) if len(g) else 0,nominal_genes_with_paired_same_mean_direction=int((np.sign(g.mean_lung_minus_bone)==np.sign(g.paired_followup_mean_effect)).sum()) if len(g) else 0)
  summaries.append(r)
 summary=pd.DataFrame(summaries);save(summary,'cohort_population_evidence_summary.tsv')
 # Parent/internal overlap shows shared observations, not additional independent evidence.
 overlap=[]
 for (layer,ds,design,parent),z in genes.groupby(['layer','dataset','design','parent']):
  main=z[z.level.eq('main')];mainset=set(main.feature)
  for target,g in z[z.level.eq('internal')].groupby('target'):
   shared=mainset&set(g.feature);same=main[['feature','mean_lung_minus_bone']].merge(g[['feature','mean_lung_minus_bone']],on='feature',suffixes=('_parent','_state'))
   overlap.append(dict(layer=layer,dataset=ds,design=design,parent=parent,internal_state=target,parent_nominal=len(mainset),state_nominal=g.feature.nunique(),shared_nominal=len(shared),shared_same_direction=int((np.sign(same.mean_lung_minus_bone_parent)==np.sign(same.mean_lung_minus_bone_state)).sum()),independent_replication=False))
 save(pd.DataFrame(overlap),'parent_internal_signal_overlap.tsv')
 # The outcome packet includes all cohort differences; no receptor/antibody shortlist is read or written.
 assert len(genes)==sum(r['nominal_up']+r['nominal_down'] for r in context_summary)
 assert genes.P.lt(.05).all() and prog.P.lt(.05).all()
 assert not any(s in str(c).lower() for c in genes.columns for s in ['surface','antibody','flow'])
 # Exact P/FDR and family preservation, including all selected features.
 verify=genes[['context_id','feature','P','FDR','family']].merge(allbase[['context_id','feature','P','FDR','family']],on=['context_id','feature'],suffixes=('_new','_source'),validate='one_to_one')
 assert np.allclose(verify.P_new,verify.P_source,equal_nan=True) and np.allclose(verify.FDR_new,verify.FDR_source,equal_nan=True) and (verify.family_new==verify.family_source).all()
 qa=dict(status='PASS',checked_at=now(),gene_rows=len(genes),program_rows=len(prog),proportion_nominal_rows=int(proportions.P.lt(.05).sum()),ORA_nominal_rows=int(ora.P.lt(.05).sum()),gene_P_FDR_families_unchanged=True,surface_or_flow_inputs_used=False,no_hypotheses_removed_for_FDR_failure=True,new_P_values_created=False,all_nominal_directions_retained=True)
 dump(O/'QA/ACCEPTANCE.json',qa)
 text='''骨—肺队列差异重新筛选：仅第一阶段（2026-10-09）

本轮已将“队列差异”与“流式检测标志”分开。输入仅为完整比例、全基因、固定程序及ORA统计；未读取膜蛋白、表面marker、抗体或流式候选表。全部原始P<0.05的有效条目均保留，不因FDR未过或缺少试剂而剔除。

对上一轮的核查
原始差异计算使用完整可测基因，随后才与表面注释取交集，因此上游全基因P/FDR没有被流式条件限定。问题在于上一轮优先展示了检测候选，把两个阶段的解释混在一起。本轮重新从完整统计表筛出队列差异，统计数值和检验家族保持不变。

第一层：比例
原有F1/F3/F4及固定内部状态的预设骨—肺比例比较，当前没有原始P<0.05的条目。该结论只涉及捕获细胞比例，不能用表达差异代替比例变化。

第二层：表达与功能
扩展F1：主来源骨8/肺3，153上调、9下调；apCAF/抗原呈递相关程序为名义升高（效应+0.3343，P=0.024242，FDR=0.393467）。唯一配对患者的该程序效应为-0.1196，方向相反；N=1仅作方向对照，不推翻非配对差异，也不能称跨设计一致。其原内部C4有162条名义基因，与父类共享141条且同向，属于共享细胞的嵌套信号，不能当作独立重复；C4的抗原呈递程序在唯一配对患者同样反向。当前不据此确认普遍的肺侧抗原呈递增强，更不能直接命名为免疫抑制CAF。
支持F3：骨6/肺2，758上调、0下调；37条固定程序名义差异主要来自Reactome。预设客户面板、Hallmark和Matrisome未出现名义差异，当前不强行归结为胶原沉积或免疫抑制主题。扩展F3骨2/肺1，缺有效重复；两层F3不能互换。
扩展F4：骨11/肺2，633上调、161下调；炎症样、干扰素响应等程序提供非配对比较中的名义功能线索。炎症样程序的非配对效应+0.3147，而两名配对患者均为负向，平均-0.1421。794条名义基因中，配对平均效应同向257条、反向496条、零效应41条；其他来源SCPCP000023（3/3患者）同向436条、反向350条、零效应8条。跨设计及来源的一致性不足，不将非配对结果直接概括为普遍的肺侧F4功能增强。非配对肺侧一名患者仅1个F4细胞，另一名642个；完整方向、效应和不一致结果均保留，不合并成一个全局P。
扩展F4_C9：骨5/肺2，123上调、14下调；胶原、基底膜及其他基质程序有名义升高，差异基因列表有相应ORA富集。137条名义基因仅6条与F4父类名义列表重合，因此独立列为内部状态的观察，不用父类结果替代。但肺侧仅两名患者，共4个细胞（1和3），且没有同状态的有效配对表达对照，属于细胞覆盖稀疏的探索性线索，不能据此确认基质屏障状态。

本轮先保留以上“群体/状态—队列差异”线索，不指定任何流式marker或组合。当前所有名义基因和患者级程序均未达到相应FDR阈值。ORA的P/FDR是基因列表注释富集模型的统计，不是独立患者功能效应P。重复数不足的比较、阴性结果和矛盾方向均保留。原分类、内部标签、严格骨肺分组和原始统计均未改变。

新增的证据整理
对全部名义基因及程序重新计算逐患者剔除后的效应范围（描述性，不产生新P），并关联去双细胞后的效应/P/FDR、同来源配对患者的方向，以及其他来源的对应效应。以上均作为证据强度记录，不设新的筛选门槛。

文件入口
cohort_population_evidence_summary.tsv：各类和状态的患者数、细胞覆盖、基因/程序数量及敏感性摘要。
cohort_nominal_genes_ALL.tsv.gz：全部名义差异基因特征，未经过膜蛋白或抗体过滤。
cohort_nominal_programs_ALL.tsv：全部名义患者级功能程序。
cohort_all_ORA_statistics.tsv.gz、cohort_nominal_ORA.tsv：单独保存的基因列表富集结果。
cohort_all_proportion_statistics.tsv、cohort_comparison_registry.tsv：比例及全部可评估/不可评估比较。
gene_other_source_evidence.tsv.gz、parent_internal_signal_overlap.tsv：跨来源方向与父子状态共享信号。
完整未筛选基因统计继续保留在Phase73 strict_v2，精确路径和校验值见INPUT_MANIFEST.json；原始结果未覆盖。
'''
 (O/'delivery/COHORT_ONLY_RESULTS_ZH.txt').write_text(text,encoding='utf-8')
 dump(O/'delivery/SUMMARY.json',qa)
 files=[f for sub in ['code','resources','tables','QA','delivery'] for f in (O/sub).rglob('*') if f.is_file() and f.suffix!='.zip']
 dump(O/'delivery/FILE_MANIFEST.json',[dict(path=str(f.relative_to(O)),bytes=f.stat().st_size,sha256=sha(f)) for f in files])
 archive=O/'delivery/CAF_cohort_only_rescreen_20261009.zip'
 with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=1) as z:
  for f in files+[O/'delivery/FILE_MANIFEST.json']:z.write(f,str(f.relative_to(O)))
 with zipfile.ZipFile(archive) as z:assert z.testzip() is None
 complete=dict(status='completed_verified',phase=O.name,completed_at=now(),archive=str(archive),archive_sha256=sha(archive),acceptance=str(O/'QA/ACCEPTANCE.json'),report=str(O/'delivery/COHORT_ONLY_RESULTS_ZH.txt'),scope='Cohort differences only; flow-marker determination deferred',**{k:v for k,v in qa.items() if k not in ['status','checked_at']})
 dump(O/'COMPLETE.json',complete);print(json.dumps(complete,ensure_ascii=False),flush=True)
if __name__=='__main__':
 try:main()
 except BaseException:
  if O.exists():(O/'exit_code').write_text('1\n')
  traceback.print_exc();raise
 else:(O/'exit_code').write_text('0\n')
