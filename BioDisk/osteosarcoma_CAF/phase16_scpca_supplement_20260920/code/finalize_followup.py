"""Deterministic evidence acceptance and client-directed descriptive outputs."""
from pathlib import Path
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
import json,gzip,itertools,hashlib,datetime
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
from scipy.io import mmread
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase16_scpca_supplement_20260920'
read=lambda p:pd.read_csv(p,sep='\t',low_memory=False)
def save(df,name):df.to_csv(O/'tables'/name,sep='\t',index=False)
panels=json.loads((P/'phase10_batch_gates_20260919/resources/panels.json').read_text())
customer=json.loads((O/'resources/customer_panels.json').read_text())
surfaces=json.loads((P/'phase8_original_plan_20260919/resources/surface_candidates_frozen.json').read_text())['surface_candidates']
genes=sorted(set(sum(customer.values(),[])+surfaces+['PTPRC','PECAM1','EPCAM']))
manifest=read(O/'manifest.tsv').set_index('scpca_sample_id');inputs=read(O/'review/input_summary.tsv')
jobs=json.loads((O/'logs/followup_jobs.json').read_text());assert all(j['passed'] for j in jobs)
def evidence(x,sy):
 x=x.tocsc();den=np.asarray(x.sum(1)).ravel();d={}
 for name,gs in panels.items():
  v=np.column_stack([np.asarray(x[:,np.flatnonzero(sy==g)].sum(1)).ravel() for g in gs if (sy==g).any()])
  d['det_'+name]=(v>=1).sum(1);d['pos_'+name]=(v>0).sum(1);d['pct_'+name]=v.sum(1)/np.maximum(den,1)*100
 return pd.DataFrame(d),den
def conflict(d,frac=.2,osteo=.2,mural=True):
 out=np.zeros(len(d),bool)
 for name,n,pct in [('Mural',3,frac),('OC',2,frac/2),('T_NK',3,frac),('Myeloid',4,frac),('Endo',3,frac),('Epi',3,frac),('Myogenic',2,frac)]:
  if name=='Mural' and not mural:continue
  pp='OC_specific' if name=='OC' else name;out|=(d['det_'+name]>=n)&(d['pct_'+pp]>pct)
 out|=(d.det_Osteo>=3)&(d.pct_Osteo>osteo)
 immune=sum(d['pct_'+k] for k in ['T_NK','B_Plasma','Myeloid','DC']);mes=d.pct_Fib+d.pct_Osteo+d.pct_Mural
 mixed=((immune>frac)&(mes>frac))|((d.pct_Endo>frac)&(mes>frac))|((d.pct_Epi>frac)&(mes>frac))
 return np.asarray(out|mixed)
def extract(x,sy):
 x=x.tocsc();v=np.full((x.shape[0],len(genes)),np.nan)
 for j,g in enumerate(genes):
  ix=np.flatnonzero(sy==g)
  if len(ix):v[:,j]=np.asarray(x[:,ix].sum(1)).ravel()
 return v
summary=[];calibration=[];expr=[];screen=[];coverage=[]
for row in inputs.to_dict('records'):
 sm=row['sample'];meta=manifest.loc[sm];a=ad.read_h5ad(O/'objects'/f'{sm}_qc.h5ad');ob=read(O/'review'/f'{sm}_cells.tsv.gz').set_index('cell_id').loc[a.obs_names]
 x=a.layers['counts'].tocsr();sy=a.var.symbol.astype(str).to_numpy();raw,lib=evidence(x,sy);raw.index=a.obs_names
 ambient=bool(row['ambient_eligible']);cnv=bool(row['cnv_eligible']);cor=None;corlib=np.zeros(len(a));corr=pd.DataFrame(np.nan,index=a.obs_names,columns=raw.columns)
 if ambient:
  assert (O/'decontx'/sm/'COMPLETE').exists()
  with gzip.open(O/'decontx'/sm/'corrected_counts.mtx.gz','rb') as f:cor=sparse.csr_matrix(mmread(f).T)
  assert cor.shape==x.shape and np.isfinite(cor.data).all() and (cor.data>=0).all()
  delta=cor-x;assert not delta.nnz or delta.data.max()<=1e-4
  corr,corlib=evidence(cor,sy);corr.index=a.obs_names
  cont=read(O/'decontx'/sm/'cell_contamination.tsv');assert cont.cell_id.tolist()==a.obs_names.tolist()
  assert np.allclose(cont.original_UMI,lib) and np.allclose(cont.corrected_UMI,corlib)
 ic=pd.Series('not_assessed',index=a.obs_names);cc=ic.copy();refpass=False
 if cnv:
  passed=[]
  for method,low in [('infercnv','low_signal'),('copykat','diploid')]:
   d=O/method/sm;assert (d/'COMPLETE').exists()
   v=read(d/'cell_results.tsv');assert v.cell_id.tolist()==a.obs_names.tolist() and v.tool_id.tolist()==ob.tool_id.tolist()
   call=v.set_index('cell_id')['call'];held=ob.cnv_role.str.startswith('holdout_');valid=call.isin(['low_signal','elevated_signal','ambiguous_signal'] if method=='infercnv' else ['diploid','aneuploid'])
   n=int((held&valid).sum());lowfrac=float(call[held&valid].eq(low).mean()) if n else 0.;ok=n>=20 and lowfrac>=.9
   calibration.append({'sample':sm,'method':method,'heldout_total':int(held.sum()),'heldout_evaluable':n,'heldout_low_fraction':lowfrac,'reference_pass':ok});passed.append(ok)
   if method=='infercnv':ic=call
   else:cc=call
  refpass=all(passed)
 dual=(ic.eq('low_signal')&cc.eq('diploid')).to_numpy();risk=(ic.eq('elevated_signal')|cc.eq('aneuploid')).to_numpy()
 anyfrac=pd.Series(np.nan,index=a.obs_names);bothfrac=anyfrac.copy()
 for cl,idx in ob.groupby('cluster',observed=True).groups.items():
  binary=ic.loc[idx].isin(['low_signal','elevated_signal'])&cc.loc[idx].isin(['diploid','aneuploid'])
  evaluable=ic.loc[idx].isin(['low_signal','elevated_signal','ambiguous_signal'])&cc.loc[idx].isin(['diploid','aneuploid'])
  if binary.any():anyfrac.loc[idx]=float((ic.loc[idx].eq('elevated_signal')|cc.loc[idx].eq('aneuploid'))[binary].mean())
  if evaluable.any():bothfrac.loc[idx]=float((ic.loc[idx].eq('elevated_signal')&cc.loc[idx].eq('aneuploid'))[evaluable].mean())
 fixed=dual&refpass&~ob.CTA_context.to_numpy(dtype=bool)&(raw.det_CTA<2).to_numpy()&(ob.n_genes>=500).to_numpy()&(raw.det_Fib>=2).to_numpy()&(anyfrac.fillna(1)<.2).to_numpy()&(bothfrac.fillna(1)<.2).to_numpy()
 def gate(frac=.2,osteo=.2,mural=True):
  return fixed&ambient&~conflict(raw,frac,osteo,mural)&~conflict(corr,frac,osteo,mural)&(corr.det_Fib>=2).to_numpy()&(corr.det_CTA<2).to_numpy()&(corlib>0)
 scope=ob.lineage_phase3.eq('FibLike').to_numpy();strict=scope&gate();osteo=scope&gate(osteo=.5);competitor=scope&gate(frac=.5,osteo=.5)
 mural=ob.lineage_phase3.isin(['FibLike','MesAmb','Mural','Mural_SMC']).to_numpy()&gate(mural=False)&(raw.det_Mural>=3).to_numpy()&(raw.pct_Mural>.2).to_numpy()
 rr=extract(x,sy);cr=extract(cor,sy) if ambient else np.full_like(rr,np.nan)
 ev=ob[['sample','patient_id','site','lineage_phase3','cluster','cnv_role','n_genes']].copy()
 ev['infer_call']=ic;ev['copykat_call']=cc;ev['reference_QC_pass']=refpass;ev['ambient_checked']=ambient;ev['unified_strict']=strict;ev['osteo_sensitive']=osteo;ev['competitor_sensitive']=competitor;ev['mural_reference']=mural;ev['cluster_any_CNV_high_fraction']=anyfrac;ev['cluster_both_CNV_high_fraction']=bothfrac
 ev['exclusion_reasons']=''
 reasons={'outside_frozen_FibLike_scope':~scope,'not_dual_CNV_low':~dual,'reference_QC_fail':np.repeat(not refpass,len(a)),'ambient_unavailable':np.repeat(not ambient,len(a)),'raw_competing_or_mixed':conflict(raw),'corrected_competing_or_mixed':conflict(corr),'Fib_below2':(raw.det_Fib<2).to_numpy()|(corr.det_Fib<2).fillna(True).to_numpy(),'n_genes_below500':(ob.n_genes<500).to_numpy(),'CTA':ob.CTA_context.to_numpy(dtype=bool)|(raw.det_CTA>=2).to_numpy()|(corr.det_CTA>=2).to_numpy(),'cluster_CNV_risk_or_unassessed':(anyfrac.fillna(1)>=.2).to_numpy()|(bothfrac.fillna(1)>=.2).to_numpy(),'corrected_counts_missing':corlib<=0}
 for reason,mask in reasons.items():ev.loc[mask,'exclusion_reasons']+=reason+';'
 assert ev.loc[strict,'exclusion_reasons'].eq('').all()
 for name,gs in customer.items():
  ix=[genes.index(g) for g in gs];cov=float(np.isfinite(rr[:,ix]).all(0).mean());support=((rr[:,ix]>=1).sum(1)>=2)&((cr[:,ix]>=1).sum(1)>=2)&(cov>=.8)
  ev[name+'_support']=support;coverage.append({'sample':sm,'panel':name,'gene_coverage':cov,'strict_CAF':int(strict.sum()),'strict_supported':int((strict&support).sum()),'mural_supported':int((mural&support).sum())})
 target=strict&ev.myCAF_customer_support.to_numpy();ev['customer_myCAF']=target
 ev.to_csv(O/'tables'/f'{sm}_identity_evidence.tsv.gz',sep='\t',index_label='cell_id')
 site='lung' if 'lung' in str(meta.tissue_location).lower() else 'bone'
 for label,mask in [('strict_CAF',strict),('myCAF',target),('other_CAF',strict&~target)]:
  for kind,v,den in [('raw',rr,lib),('corrected',cr,corlib)]:
   total=float(den[mask].sum())
   for j,g in enumerate(genes):
    measured=bool(np.isfinite(v[:,j]).all());count=float(v[mask,j].sum()) if measured else np.nan
    expr.append({'project':meta.scpca_project_id,'sample':sm,'patient':meta.participant_id,'site':site,'set':label,'expression':kind,'gene':g,'cells':int(mask.sum()),'library_total':total,'gene_counts':count,'detected_cells':int((v[mask,j]>=1).sum()) if measured else np.nan,'measured':measured})
 dumpmask=(rr[:,[genes.index(g) for g in ['PTPRC','PECAM1','EPCAM']]]<1).all(1)
 for threshold in [1,2]:
  for pair in [(g,) for g in surfaces]+list(itertools.combinations(surfaces,2)):
   ix=[genes.index(g) for g in pair]
   if not np.isfinite(rr[:,ix]).all():continue
   gatecells=dumpmask&(rr[:,ix]>=threshold).all(1)
   screen.append({'project':meta.scpca_project_id,'sample':sm,'patient':meta.participant_id,'site':site,'markers':'+'.join(pair),'threshold':threshold,'all_QC_cells':len(a),'target_cells':int(target.sum()),'gate_cells':int(gatecells.sum()),'target_in_gate':int((gatecells&target).sum()),'CAF_in_gate':int((gatecells&strict).sum()),'mural_in_gate':int((gatecells&mural).sum()),'CNV_risk_in_gate':int((gatecells&risk).sum()),'unresolved_in_gate':int((gatecells&~(strict|mural|risk)).sum())})
 summary.append({'sample':sm,'project':meta.scpca_project_id,'patient':meta.participant_id,'site':site,'QC_cells':len(a),'ambient_evaluable':ambient,'CNV_evaluable':cnv,'reference_QC_pass':refpass,'frozen_FibLike_scope':int(scope.sum()),'strict_CAF':int(strict.sum()),'myCAF':int(target.sum()),'osteo_sensitive':int(osteo.sum()),'competitor_sensitive':int(competitor.sum()),'mural_reference':int(mural.sum()),'interpretation':'operational candidates' if cnv and ambient and refpass else 'identity_not_evaluable; zero is not biological absence'})
 save(pd.DataFrame(summary),'identity_summary.tsv');print('IDENTITY_COMPLETE',json.dumps(summary[-1]),flush=True)
save(pd.DataFrame(calibration),'CNV_heldout_acceptance.tsv');save(pd.DataFrame(coverage),'customer_panel_summary.tsv')
ex=pd.DataFrame(expr);save(ex,'sample_gene_expression.tsv.gz')
keys=['project','patient','site','set','expression','gene'];values=['cells','library_total','gene_counts','detected_cells']
don=ex.groupby(keys,dropna=False)[values].sum(min_count=1).reset_index();don['CP10k']=don.gene_counts/don.library_total.replace(0,np.nan)*1e4;don['detection_fraction']=don.detected_cells/don.cells.replace(0,np.nan);save(don,'patient_site_gene_expression.tsv.gz')
contrasts=[]
for key,g in don.groupby(['project','patient','set','expression','gene']):
 lung=g[g.site=='lung'];bone=g[g.site=='bone']
 if len(lung)!=1 or len(bone)!=1:continue
 l=lung.iloc[0];b=bone.iloc[0];eligible=bool(l.cells>=20 and b.cells>=20 and np.isfinite(l.CP10k) and np.isfinite(b.CP10k))
 contrasts.append(dict(zip(['project','patient','set','expression','gene'],key),lung_cells=int(l.cells),bone_cells=int(b.cells),eligible=eligible,delta_log2_CP10k=float(np.log2(l.CP10k+.1)-np.log2(b.CP10k+.1)) if eligible else np.nan,interpretation='same-patient site contrast; collection timing unresolved; no causal/metastatic-risk claim'))
paired=pd.DataFrame(contrasts);save(paired,'paired_site_gene_effects.tsv')
ss=pd.DataFrame(screen);save(ss,'surface_screen_sample.tsv.gz')
cols=['all_QC_cells','target_cells','gate_cells','target_in_gate','CAF_in_gate','mural_in_gate','CNV_risk_in_gate','unresolved_in_gate']
sd=ss.groupby(['project','patient','site','markers','threshold'])[cols].sum().reset_index()
for name,num,den in [('target_coverage','target_in_gate','target_cells'),('known_target_fraction','target_in_gate','gate_cells'),('CNV_risk_fraction','CNV_risk_in_gate','gate_cells'),('mural_fraction','mural_in_gate','gate_cells'),('unresolved_fraction','unresolved_in_gate','gate_cells')]:sd[name]=sd[num]/sd[den].replace(0,np.nan)
save(sd,'surface_screen_patient.tsv.gz')
eligible=sd[(sd.target_cells>=20)&sd.site.eq('lung')];ranks=[]
for (markers,t),g in eligible.groupby(['markers','threshold']):
 ranks.append({'markers':markers,'threshold':t,'eligible_patient_records':len(g),'median_target_coverage':g.target_coverage.median(),'median_known_target_fraction':g.known_target_fraction.median(),'worst_CNV_risk_fraction':g.CNV_risk_fraction.max(),'worst_unresolved_fraction':g.unresolved_fraction.max(),'score':(g.target_coverage*g.known_target_fraction.fillna(0)).median(),'interpretation':'RNA exploratory enrichment screen; not protein sorting purity or validated metastasis specificity'})
rank=pd.DataFrame(ranks)
if len(rank):rank=rank.sort_values(['threshold','score'],ascending=[True,False])
save(rank,'surface_candidate_ranking.tsv')
summarydf=pd.DataFrame(summary);npaired=int(paired[paired.eligible&paired['set'].eq('myCAF')][['project','patient']].drop_duplicates().shape[0]) if len(paired) else 0
lines=['# ScPCA 补充分析结果','',f'共处理12个样本，质控后44,242个细胞。严格规则接受 {summarydf.strict_CAF.sum()} 个操作性 CAF 候选，其中客户面板支持的 myCAF 为 {summarydf.myCAF.sum()} 个。',f'满足骨、肺两侧各至少20个 myCAF 的患者记录：{npaired} 个。同一患者多个肺标本合并后计算，未把细胞或标本当作独立患者。','','## 样本结果','', '|样本|QC细胞|严格CAF候选|myCAF|身份检查可用|','|---|---:|---:|---:|---|']
for r in summary:lines.append(f"|{r['sample']}|{r['QC_cells']}|{r['strict_CAF']}|{r['myCAF']}|{'是' if r['ambient_evaluable'] and r['CNV_evaluable'] and r['reference_QC_pass'] else '不足，零值不代表生物学缺失'}|")
lines+=['','## 客户目标达成情况','','已完成新增样本的身份筛选、客户五类面板、患者内骨—肺表达方向比较和表面标志物单基因/双基因组合筛选。结果为操作性候选与探索性效应量，不报告细胞层面的疾病差异 P 值。','当前数据不支持将这些候选称为“已验证的转移相关亚型”，也不能把 RNA 门控比例解释为实际流式分选纯度。最终分选方案仍需表面蛋白验证。','','## 主要结果文件','','- identity_summary.tsv：样本身份验收及CAF/myCAF数量。','- CNV_heldout_acceptance.tsv：两种CNV方法的留出正常参考检查。','- customer_panel_summary.tsv：myCAF、iCAF、apCAF、成骨、矿化面板覆盖与支持。','- paired_site_gene_effects.tsv：同患者骨—肺方向和效应量，包含不满足比较条件的记录。','- surface_candidate_ranking.tsv：肺样本候选RNA组合排序。','- 各样本 identity_evidence.tsv.gz：逐细胞保留与排除原因。','','## 解释边界','','初始大类标签来自冻结的标志物规则，边界细胞保守排除；MesAmb和周细胞不能凭收缩基因自动并入CAF。使用同一样本 broad-lineage DecontX 无背景敏感性校正，尚未进行空液滴背景模型比较。CNV双低不能独立证明肿瘤相关性，尤其肺正常成纤维细胞。项目17与项目23的患者编号未与全部旧GEO数据完成去重，因此不增加跨队列独立患者数。配对仅确认相同患者编号，未确立同步采样，未将其直接标为临床A/B1比较。', '原始数据、此前3976个CAF集合和原分析结果均未覆盖。']
if len(rank):
 lines+=['','## RNA筛选靠前的组合','']
 for r in rank[rank.threshold.eq(1)].head(5).to_dict('records'):lines.append(f"- {r['markers']}：{r['eligible_patient_records']}个患者记录，目标覆盖中位数{r['median_target_coverage']:.1%}，门内已知目标比例中位数{r['median_known_target_fraction']:.1%}。")
else:lines+=['','没有满足≥20个 myCAF 的肺患者记录可用于组合排序。']
(O/'delivery/补充分析结果.md').write_text('\n'.join(lines)+'\n')
qa={'source_QC_samples':12,'source_QC_cells':44242,'reviewed_cells':int(summarydf.QC_cells.sum()),'job_outputs_passed':all(j['passed'] for j in jobs),'strict_CAF':int(summarydf.strict_CAF.sum()),'myCAF':int(summarydf.myCAF.sum()),'eligible_paired_myCAF_records':npaired,'all_eligible_computations_complete':True,'final_client_goal_achieved':False,'validated_subtypes':0}
assert qa['reviewed_cells']==44242 and len(summarydf)==12 and (summarydf.myCAF<=summarydf.strict_CAF).all()
(O/'delivery/QA.json').write_text(json.dumps(qa,indent=2))
hashes={}
for path in list((O/'tables').glob('*'))+list((O/'delivery').glob('*')):
 if path.is_file():hashes[str(path.relative_to(O))]=hashlib.sha256(path.read_bytes()).hexdigest()
(O/'delivery/SHA256.json').write_text(json.dumps(hashes,indent=2))
(O/'ANALYSIS_COMPLETE').write_text('ALL_ELIGIBLE_ANALYSES_COMPLETE_WITH_DOCUMENTED_DATA_LIMITATIONS\n')
print('ANALYSIS_COMPLETE',json.dumps(qa),flush=True)
