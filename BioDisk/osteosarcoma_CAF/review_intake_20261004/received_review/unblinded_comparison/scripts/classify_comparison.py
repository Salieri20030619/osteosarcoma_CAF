import pathlib,pandas as pd,numpy as np,json
P=pathlib.Path(__file__).parent;O=P/'review_output';C=O/'unblinded_comparison';B=O/'blind_submission';m=pd.read_csv(C/'joined_cluster_evidence.tsv',sep='\t');bd=pd.read_csv(B/'cluster_review_details.tsv',sep='\t');cn=pd.read_csv(C/'proposed_CNV_cluster_composition.tsv',sep='\t');
coarse={**dict.fromkeys(['T','Treg','NK','T_NK'],'T_NK'),**dict.fromkeys(['B','Plasma','MixBPlasma'],'B_Plasma'),**dict.fromkeys(['Mono','Macro','Myeloid','MonoDC','DC','pDC'],'Myeloid'),**dict.fromkeys(['Endo','Lymph'],'Endothelial'),**dict.fromkeys(['Mural','SMC'],'Mural_like'),'OC':'Osteoclast_like','Mast':'Mast','Fibro':'Fibroblast_like','Osteo':'Osteogenic_like','Chondro':'Chondrocyte_like','Osteochondro':'Osteochondrogenic_like','Mesench':'Mesenchymal_unresolved','Myogenic':'Myogenic_like','Epi':'Epithelial','Renal':'Epithelial','Ery':'Erythroid','U':'Unresolved'}
alias={'Epithelial_like':'Epithelial','Endothelial_like':'Endothelial','Erythroid_like':'Erythroid'}
cat=[];actions=[];priority=[]
for r in m.itertuples():
 b=coarse.get(r.lineage_code,r.lineage_code);p=alias.get(r.proposed_lineage,r.proposed_lineage)
 if b==p:
  status='未定范围一致' if b in ['Unresolved','Mesenchymal_unresolved'] else '主要谱系方向一致'
  action='保留未定范围及逐簇复核事项' if status=='未定范围一致' else '可保留主要表达谱系；盲审亚型与状态作为细化建议'
 elif p=='Immune_unresolved' and r.lineage_code.startswith('Mix'):
  status='未定范围相容，盲审明确混合成分';action='保留免疫未定，同时记录盲审指出的混合程序'
 elif p=='Unresolved' and r.lineage_code.startswith('Mix'):
  status='未定范围相容，盲审明确混合成分';action='保留未定并加入混合表达原因'
 elif r.lineage_code.startswith('Mix'):
  status='新增混合表达风险';action='建议新增mixed_flag；在完成原始共表达复核前不将全簇作为纯谱系使用'
 elif r.lineage_code=='U':
  status='盲审更保守';action='保留未定，不强制赋予原拟定的范围标签'
 elif p in ['Unresolved','Mesenchymal_unresolved','Immune_unresolved','Neural_myogenic_unresolved']:
  status='原未定簇获得更细候选方向';action='保留原未定状态作为备选，加入盲审候选方向和把握；不自动升级为已确认身份'
 elif p=='Osteochondrogenic_like' and b in ['Osteogenic_like','Chondrocyte_like']:
  status='骨软骨范围相容，亚型取舍不同';action='可共同保留骨软骨间充质上位范围，细分用具体成组marker和患者证据复核'
 else:
  status='主要谱系不相容';action='优先复核原始计数、样本内表达与簇映射，暂不统一标签'
 if r.lineage_code=='Fibro':action+='；仅为成纤维样候选，非已确认CAF'
 if 'CNV支持恶性候选' in r.malignancy_assessment:action+='；恶性候选须另行验证，不按表达谱系回填正常基质'
 if r.lineage_code.startswith('Mix') or (r.lineage_code in ['Fibro','Mesench'] and 'CNV支持' in r.malignancy_assessment):pri='高'
 elif status in ['原未定簇获得更细候选方向','骨软骨范围相容，亚型取舍不同','盲审更保守'] or r.lineage_confidence=='低':pri='中'
 else:pri='常规'
 cat.append(status);actions.append(action);priority.append(pri)
m['blind_coarse_lineage']=m.lineage_code.map(coarse).fillna(m.lineage_code);m['proposed_coarse_lineage']=m.proposed_lineage.replace(alias);m['comparison_category']=cat;m['recommended_action']=actions;m['review_priority']=priority
cv=cn.pivot(index='cluster_id',columns='CNV_evidence_current',values='fraction_in_cluster').fillna(0).add_prefix('proposed_CNV_fraction_').reset_index();m=m.merge(cv,on='cluster_id',validate='one_to_one');m=m.merge(bd[['cluster_id','copykat_coverage','reference_comparable_patients','copykat_fold_patient_median']],on='cluster_id',validate='one_to_one');m['malignancy_comparison']='拟定包未作病理恶性确证；盲审亦未作确证。分类call与连续/区域摘要口径不同，不直接计算恶性一致率。'
cols=['source_id','dataset','cluster_id','original_cluster','curation_id','cells','patients','proposed_lineage','status','reviewed_lineage','lineage_code','lineage_confidence','comparison_category','review_priority','state_proposed','state_blind','malignancy_assessment','malignancy_comparison','copykat_coverage','reference_comparable_patients','copykat_fold_patient_median','recommended_action','supporting_genes','evidence_genes','BlueprintEncode_label','HPCA_label','BlueprintEncode_pruned_label','HPCA_pruned_label','genomic_evidence','uncertainty','reviewer_notes']+[c for c in m if c.startswith('proposed_CNV_fraction_')]
m[cols].to_csv(C/'cluster_annotation_comparison.tsv',sep='\t',index=False,encoding='utf-8-sig');m[m.comparison_category!='主要谱系方向一致'][cols].to_csv(C/'differences_and_unresolved.tsv',sep='\t',index=False,encoding='utf-8-sig')
summary=m.groupby('comparison_category').agg(clusters=('cluster_id','size'),cells=('cells','sum')).reset_index();summary['cluster_fraction']=summary.clusters/307;summary['cell_fraction']=summary.cells/876648;summary.to_csv(C/'comparison_summary.tsv',sep='\t',index=False);m.groupby(['source_id','dataset','comparison_category']).agg(clusters=('cluster_id','size'),cells=('cells','sum')).reset_index().to_csv(C/'comparison_by_source.tsv',sep='\t',index=False)
m.groupby(['proposed_lineage','lineage_code']).agg(clusters=('cluster_id','size'),cells=('cells','sum')).reset_index().to_csv(C/'lineage_cross_table.tsv',sep='\t',index=False)
fb=m[(m.proposed_lineage=='Fibroblast_like')|(m.lineage_code=='Fibro')];fb[cols].to_csv(C/'fibroblast_candidate_comparison.tsv',sep='\t',index=False,encoding='utf-8-sig')
harmonize=pd.DataFrame([{'blind_lineage_code':k,'comparison_parent_label':v,'scope':'比较用途，未改变盲审标签；亚型细分不作主要谱系错误'} for k,v in coarse.items()]);harmonize.to_csv(C/'label_harmonization.tsv',sep='\t',index=False)
print(summary.to_string(index=False));print('FIBRO',fb[['cluster_id','dataset','original_cluster','cells','proposed_lineage','lineage_confidence']].to_string(index=False));print('MIXNEW',m[m.comparison_category=='新增混合表达风险'][['cluster_id','dataset','original_cluster','cells','proposed_lineage','lineage_code']].to_string(index=False));print('REFINE',m[m.comparison_category=='原未定簇获得更细候选方向'][['cluster_id','cells','proposed_lineage','lineage_code','lineage_confidence']].to_string(index=False))
