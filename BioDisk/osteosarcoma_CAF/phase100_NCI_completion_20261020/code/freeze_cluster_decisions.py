import pathlib,json,datetime
import pandas as pd
root=pathlib.Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');p=root/'phase100_NCI_completion_20261020';q=root/'phase95_NCI_expanded_reanalysis_20261016'
decisions=json.loads((p/'resources/EXPLICIT_CLUSTER_DECISIONS.json').read_text(encoding='utf-8'))
df=pd.read_csv(p/'tables/cluster_review_evidence.tsv',sep='\t',encoding='utf-8')
reason={'M':'基质/成纤维表达背景；仅保留机械门控内细胞作为身份未决拓展候选，作者肿瘤或收缩谱系冲突另列。','O':'成骨、神经样或作者肿瘤背景与成纤维表达重叠；沿用历史拓展层口径保留符合门控的子集，明确标为成骨/肿瘤身份冲突。','C':'细胞周期主导，缺少明确其他宿主谱系主导证据；保留符合门控的成纤维表达子集，标记增殖及肿瘤身份未决。','S':'应激/线粒体/核糖体背景；候选特异表达核查支持COL1A1/COL1A2/FN1等基质信号，保留门控内子集并标记身份不确定。','X':'混合基质群：候选有胶原/基质信号，同时存在其他谱系标志；保留拓展候选，混合谱系风险显式标记。','E':'ADGRL4/FLT1/CALCRL/VWF/PECAM1等内皮或血管亚型标志协调主导，排除偶发成纤维评分。','I':'MSR1/SLC16A10/PTPRC/PLXDC2等髓系标志协调主导，排除竞争谱系群中的成纤维评分。','T':'PTPRC/SKAP1/THEMIS/CD247等淋巴细胞标志主导，排除。','P':'SFTPB/肺泡或纤毛/上皮程序主导，排除。','V':'PRKG1/PDE3A/CCDC102B/RGS5/PDGFRB/TRPC6/ABCC9等血管壁程序协调主导，缺乏足够群级间质成纤维特异证据，排除。','K':'ACP5/MMP9/CTSK/TCIRG1等破骨样程序主导，排除。'}
rows=[]
for sample,z in df.groupby('sample'):
 if not bool(z.reference_pass.iloc[0]):
  for r in z.itertuples():rows.append(dict(sample=sample,cluster=int(r.cluster),review_class='REFERENCE_NOT_RELEASED',decision='not_admitted_reference_validation_failed_or_insufficient',reason='该文库未通过冻结的正常参考诊断或参考不足；不放宽门控纳入。群身份不作为挽救失败参考的依据。',top_markers=r.top_markers,mechanical_candidates=r.mechanical_candidates,review_scope='library eligibility disposition; no claim of completed cluster identity validation'))
  continue
 groups=decisions[sample.replace('NCI__','')];mapping={}
 for cat,cls in groups.items():
  for cl in cls:assert cl not in mapping;mapping[cl]=cat
 assert set(mapping)==set(z.cluster.astype(int)),sample
 for r in z.itertuples():
  cat=mapping[int(r.cluster)];keep=cat in ['M','O','C','S','X']
  rows.append(dict(sample=sample,cluster=int(r.cluster),review_class=cat,decision='retain_expanded_identity_unresolved' if keep else 'exclude_competing_lineage',reason=reason[cat],top_markers=r.top_markers,mechanical_candidates=r.mechanical_candidates,review_scope='AI-assisted explicit review of cluster and candidate marker evidence; no independent pathologist review'))
out=pd.DataFrame(rows);assert len(out)==len(df) and not out.duplicated(['sample','cluster']).any()
out['reviewed_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
out.to_csv(p/'tables/cluster_review_decisions.tsv',sep='\t',encoding='utf-8',index=False)
print(out.groupby(['decision','review_class']).mechanical_candidates.agg(['size','sum']).to_string())
