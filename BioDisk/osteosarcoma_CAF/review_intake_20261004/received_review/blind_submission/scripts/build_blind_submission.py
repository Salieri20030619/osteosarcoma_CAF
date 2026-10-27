import pathlib,json,hashlib,datetime,shutil,platform
import pandas as pd,numpy as np
P=pathlib.Path(__file__).parent;O=P/'review_output';F=O/'blind_submission';F.mkdir(exist_ok=True)
dec=json.loads((O/'reviewer_decisions.json').read_text()); st=pd.read_csv(O/'cluster_structure_review.tsv',sep='\t').set_index('cluster_id'); qc=pd.read_csv(O/'cluster_expression_qc.tsv',sep='\t').set_index('cluster_id');cnv=pd.read_csv(O/'cnv_cluster_review.tsv',sep='\t').set_index('cluster_id');cv=pd.read_csv(O/'regional_CNV_coverage.tsv',sep='\t');reg=pd.read_csv(O/'regional_CNV_summary.tsv.gz',sep='\t');ev=pd.read_csv(O/'selected_gene_evidence.tsv.gz',sep='\t'); ev['rank_support']=ev.detection_fraction*(ev.patient_mean_log2_difference.clip(lower=0)+.1); ev=ev.sort_values('rank_support',ascending=False).drop_duplicates(['cluster_id','gene']); mi=pd.read_csv(O/'mixed_program_evidence.tsv',sep='\t');alt=pd.read_csv(O/'alternative_reference_sensitivity.tsv',sep='\t');
labels={'T':'T细胞','Treg':'T细胞（调节/活化T样）','NK':'NK细胞为主','T_NK':'T/NK淋巴细胞（亚群未分）','B':'B细胞为主','Plasma':'浆细胞样','DC':'树突状细胞（亚型见备注）','pDC':'浆细胞样树突状细胞方向','Mono':'单核细胞','Macro':'巨噬细胞','Myeloid':'髓系细胞','MonoDC':'单核/树突状髓系细胞（边界未定）','OC':'破骨细胞样','Mast':'肥大细胞','Endo':'血管内皮细胞','Lymph':'淋巴管内皮样','Mural':'周细胞/血管壁细胞','SMC':'平滑肌细胞','Fibro':'成纤维细胞样','Osteo':'成骨样间充质细胞','Chondro':'软骨样间充质细胞','Osteochondro':'骨/软骨双分化间充质细胞','Mesench':'间充质细胞（成纤维/成骨样边界未定）','Myogenic':'肌源性细胞','Epi':'上皮细胞（组织方向见备注）','Renal':'肾小管/集合管上皮样','Ery':'红系细胞','U':'未定谱系','MixPlatelet':'血小板/髓系混合信号','MixTMacro':'T/髓系混合信号','MixOsteoMacro':'成骨/巨噬混合信号','MixMyeloidStroma':'髓系/基质混合信号','MixEndoEpi':'内皮/肺上皮混合信号','MixPlasmaDC':'浆细胞/pDC混合信号','MixBPlasma':'B/浆细胞混合程序','MixBPlasmaDC':'B/浆细胞/pDC混合信号'}
pools={'T':'CD3D CD3E CD3G TRAC TRBC1 TRBC2 CD2 LCK','Treg':'FOXP3 CTLA4 IL2RA TIGIT TNFRSF18 CD3D TRAC CD2','NK':'GNLY NKG7 KLRD1 PRF1 KLRF1 XCL1 XCL2','T_NK':'CD3D CD3E TRBC2 TRAC GNLY NKG7 KLRD1','B':'MS4A1 CD79A CD79B CD22 BANK1 PAX5','Plasma':'JCHAIN MZB1 DERL3 IGKC TNFRSF17 TENT5C XBP1','DC':'CD1C FCER1A CLEC10A CLEC9A XCR1 LAMP3 IRF8 WDFY4 CLEC4C IL3RA','pDC':'LILRA4 CLEC4C IL3RA GZMB PLD4 IRF8 TCF4','Mono':'FCN1 S100A8 S100A9 LYZ LST1 CD300E CTSS','Macro':'C1QA C1QB C1QC CD163 CSF1R MSR1 MRC1 LST1 TYROBP','Myeloid':'CSF1R AIF1 TYROBP LYZ LST1 CD163 FCER1G ITGAX TFEC SLC16A10','MonoDC':'FCN1 LYZ LST1 CD1C FCER1A CLEC10A CTSS CSF1R','OC':'ACP5 CTSK MMP9 ATP6V0D2 SLC9B2 TCIRG1 SIGLEC15','Mast':'TPSB2 TPSAB1 CPA3 KIT MS4A2 HDC','Endo':'PECAM1 CDH5 VWF KDR EMCN CLDN5 FLT1 ERG','Lymph':'PROX1 FLT4 CCL21 LYVE1 RELN STAB2 MMRN1','Mural':'RGS5 PDGFRB NOTCH3 MCAM ABCC9 KCNJ8','SMC':'MYH11 CNN1 ACTA2 TAGLN LMOD1 MYLK','Fibro':'DCN LUM COL3A1 DPT FBLN1 COL14A1 PDGFRA POSTN FGF7 ITGBL1 COL6A3 SFRP2','Osteo':'RUNX2 SATB2 SP7 PTH1R IBSP IFITM5 ALPL DMP1 BGLAP MEPE','Chondro':'COL2A1 ACAN SOX9 COL9A1 COL9A2 COL9A3 COL10A1 COMP','Osteochondro':'COL2A1 COL9A1 COL9A3 IFITM5 SP7 PTH1R','Mesench':'DCN LUM COL3A1 COL6A3 PDGFRA POSTN RUNX2 SATB2 PTH1R IBSP','Myogenic':'MYOD1 MYOG MYF5 PAX7 ACTA1 TNNT3 MYH3 MYH7 NEB TTN RYR1','Epi':'EPCAM KRT8 KRT18 KRT19 KRT7 CDH1 SFTPB SFTPC SLC34A2 DNAH12 DNAH5 CFAP157 RSPH1 AGER CLIC5 EMP2','Renal':'PAX8 LRP2 CUBN AQP2 FXYD4 SLC4A9 SLC26A4 SLC17A1','Ery':'HBB HBA1 HBA2 ALAS2 AHSP GYPA','MixPlatelet':'PF4 PPBP GP9 C1QA TYROBP','MixTMacro':'CD3D TRAC CD3E C1QA CSF1R CD163','MixOsteoMacro':'RUNX2 IFITM5 SATB2 CD163 C1QA CSF1R','MixMyeloidStroma':'C1QA CD163 CSF1R DCN LUM COL3A1','MixEndoEpi':'PECAM1 CDH5 VWF SFTPB SLC34A2 SFTPA1','MixPlasmaDC':'IGKC JCHAIN MZB1 LILRA4 PLD4 GZMB','MixBPlasma':'MS4A1 PAX5 CD79A JCHAIN MZB1 IGKC','MixBPlasmaDC':'MS4A1 CD79A JCHAIN MZB1 LILRA4 GZMB'}
mesench={'Fibro','Mesench','Osteo','Chondro','Osteochondro','Myogenic'}
rows=[];detail=[]; usedgenes=[];reviewlog=[]
for cl,d in dec.items():
 r=st.loc[cl];q=qc.loc[cl];g=ev[ev.cluster_id==cl];label=d['label']; conf=d['confidence'];cand=g[g.gene.isin(pools.get(label,'').split())].copy();
 if label=='U':
  source_top=pd.read_csv(P/r.source_id/'top_descriptive_markers.tsv',sep='\t');names=source_top[source_top.cluster_id==cl].head(6).symbol.tolist(); cand=g[g.gene.isin(names)].copy()
 if label.startswith('Mix') or label=='Osteochondro':
  cand['order']=cand.gene.map({n:i for i,n in enumerate(pools[label].split())});cand=cand.sort_values('order')
 else:cand=cand.sort_values('rank_support',ascending=False)
 cand=cand[cand.detection_fraction>0].head(6)
 sup=[]
 for gg in cand.itertuples():
  sup.append(f'{gg.gene}[检出{gg.detection_fraction:.1%},Δlog2={gg.patient_mean_log2_difference:.2f},正向{int(gg.positive_patients)}/{int(gg.N_patients_matched)}患者]');usedgenes.append({'cluster_id':cl,'gene':gg.gene,'feature_id':gg.feature_id})
 support=('仅诊断状态/技术程序，非已确认谱系marker：' if label=='U' else '')+'；'.join(sup)
 cc=cnv.loc[cl] if cl in cnv.index else None;regional=cv[cv.cluster_id==cl];mtexts=[]
 for method,rreg in regional.groupby('method'):
  mtexts.append(f'{method}:分组覆盖{int(rreg.assigned_cells.sum())}/{int(r.cells)}细胞,{rreg.sample_id.nunique()}样本')
 ctext='CopyKAT未提供该簇有效连续谱' if cc is None else f'CopyKAT覆盖{int(cc.copykat_cells)}/{int(r.cells)}细胞({cc.copykat_cells/r.cells:.1%}),涉及{int(cc.copykat_patients)}患者'
 if cc is not None and cc.reference_comparable_patients>0:
  ctext+=f'；满足描述性比较条件{int(cc.reference_comparable_patients)}患者，较同样本表达参照幅度中位比{cc.median_fold:.2f}(范围{cc.min_fold:.2f}–{cc.max_fold:.2f})，≥2倍{int(cc.patients_ge2)}/{int(cc.reference_comparable_patients)}患者'
 else:ctext+='；同样本参照比较不可评估'
 ctext+='；'+('；'.join(mtexts) if mtexts else '区域CNV未提供可匹配记录')
 if cc is not None and cc.direction_comparisons>0:ctext+=f'；主方法粗粒度方向一致{cc.direction_agreement:.1%}/{int(cc.direction_comparisons)}个样本-簇-染色体比较(非独立样本量)'
 else:ctext+='；跨方法区域方向不足以比较'
 strong=cc is not None and cc.reference_comparable_patients>0 and cc.median_fold>=1.5 and cc.direction_comparisons>=3 and cc.direction_agreement>=.8 and cc.copykat_cells/r.cells>=.5
 if label=='U':mal='不可判定：谱系及/或基因组证据不足'
 elif label.startswith('Mix'):mal='未定：混合表达下CNV不能可靠归属单一细胞身份'
 elif label in mesench:
  if strong:mal='CNV支持恶性候选（探索性，非病理确诊）'
  elif cc is not None and cc.patients_ge2>0:mal='部分患者存在CNV异常线索；全簇恶性身份未定'
  else:mal='恶性身份未定；当前CNV不足以确认或排除'
 elif strong or (cc is not None and cc.patients_ge2>0):mal='存在需核查的CNV异常/谱系冲突；不据此改为恶性'
 else:mal='当前数据不支持确立恶性身份；未作正常/二倍体确证'
 altc=alt[alt.cluster_id==cl];strongalt=altc[altc.both_reported & (altc.immune_endothelial_reference.abs()>=10)&(altc.lymphoid_reference.abs()>=10)];altconf=int((strongalt.direction_agrees_both_reported==0).sum());
 if len(altc):ctext+=f'；替代参考：{altc.sample_id.nunique()}样本，双参考均报告且|净区域|≥10Mb的{len(strongalt)}项中{altconf}项方向相反'
 risks=[f'谱系把握：{conf}']
 if r.N_patients_matched<10:risks.append(f'匹配患者{int(r.N_patients_matched)}，描述性/探索性证据')
 if r.max_patient_fraction>=.8:risks.append(f'患者集中：{r.largest_patient}占{r.max_patient_fraction:.1%}')
 if r.cells<100:risks.append('稀有小簇')
 if q.matrix_mt_fraction>=.3:risks.append(f'簇汇总线粒体计数占比{q.matrix_mt_fraction:.1%}(非细胞中位数)')
 if q.counts_per_cell<1000:risks.append(f'簇计数总量/细胞数={q.counts_per_cell:.0f}，需逐细胞复杂度核查')
 if cc is None or cc.copykat_cells/r.cells<.5:risks.append('CopyKAT覆盖不足；不可将缺失视为正常')
 if r.r05_major_fraction<.8:risks.append('主簇跨低分辨率群，边界不稳定')
 if r.r20_n_ge10pct>=2:risks.append('高分辨率存在多个实质子群')
 if label.startswith('Mix') or label in ['U','Mesench']:risks.append('需要逐细胞原始计数/共表达及样本内细分')
 if altconf:risks.append('CNV对替代参考敏感')
 if r.source_id=='D08':risks.append('公共feature计数与完整文库分母不同，逐样本原始计数未复核')
 state=d['state'] or '未见足以替代主要谱系命名的主导状态程序'
 notes=f'状态：{state}。患者{int(r.patients)}，样本{int(r.samples)}；最大患者{r.largest_patient}占{r.max_patient_fraction:.1%}。r0.5主要归属{r.r05_major}({r.r05_major_fraction:.1%},反向占比{r.r05_major_reciprocal:.1%})；r2.0主要子群{r.r20_major}({r.r20_major_fraction:.1%})，占主簇≥10%的子群{int(r.r20_n_ge10pct)}个。'+d['note']
 mm=mi[mi.cluster_id==cl]
 if len(mm):
  best=mm.sort_values('minimum_codetection_fraction',ascending=False).iloc[0]
  if best.minimum_codetection_fraction>0:notes+=f' {best.gene_A}/{best.gene_B}同一条形码共检出的数学下界{best.minimum_codetection_fraction:.1%}；原因未定。'
 rows.append({'source_id':r.source_id,'cluster_id':cl,'cells':int(r.cells),'patients':int(r.patients),'reviewed_lineage':labels[label],'malignancy_assessment':mal,'supporting_genes':support,'genomic_evidence':ctext,'uncertainty':'；'.join(risks),'reviewer_notes':notes})
 detail.append({'source_id':r.source_id,'cluster_id':cl,'lineage_code':label,'lineage_confidence':conf,'state':d['state'] or '未指定','cells':int(r.cells),'patients':int(r.patients),'samples':int(r.samples),'matched_patients':int(r.N_patients_matched),'max_patient_fraction':r.max_patient_fraction,'copykat_coverage':0 if cc is None else cc.copykat_cells/r.cells,'reference_comparable_patients':0 if cc is None else int(cc.reference_comparable_patients),'copykat_fold_patient_median':np.nan if cc is None else cc.median_fold,'malignancy_assessment':mal,'alternative_reference_opposite_regions':altconf,'CAF_handling':('成纤维样候选；先复核恶性身份及患者分层' if label=='Fibro' else '成纤维/成骨边界待核查，不直接纳入CAF' if label=='Mesench' else '本次未纳入成纤维样候选集合')})
 reviewlog.append({'cluster_id':cl,'independent_lineage_decision':label,'full_marker_panel_reviewed':True,'patient_representation_reviewed':True,'resolution_structure_reviewed':True,'CNV_coverage_and_direction_reviewed':True})
res=pd.DataFrame(rows);template=pd.read_csv(P/'reviewer_response_template.tsv',sep='\t');res=res.set_index('cluster_id').loc[template.cluster_id].reset_index()[template.columns];assert res.shape==(307,10);assert res.cells.sum()==876648;assert res.cluster_id.is_unique;assert not res.isna().any().any();assert (res.astype(str).apply(lambda x:x.str.len()).min()>0).all();
res.to_csv(F/'reviewer_response_template.tsv',sep='\t',index=False,encoding='utf-8-sig');pd.DataFrame(detail).to_csv(F/'cluster_review_details.tsv',sep='\t',index=False,encoding='utf-8-sig');pd.DataFrame(reviewlog).to_csv(F/'review_coverage.tsv',sep='\t',index=False);pd.DataFrame(usedgenes).to_csv(F/'supporting_gene_index.tsv',sep='\t',index=False)
# All key fields must match the blank supplied template; no clusters dropped.
assert res[['source_id','cluster_id','cells','patients']].equals(template[['source_id','cluster_id','cells','patients']])
summary=pd.DataFrame(detail).groupby('lineage_code').agg(clusters=('cluster_id','size'),cells=('cells','sum')).sort_values('cells',ascending=False);summary['lineage']=summary.index.map(labels);summary.to_csv(F/'lineage_summary.tsv',sep='\t')
pd.DataFrame(detail).query('lineage_code in ["Fibro","Mesench"]').merge(res[['cluster_id','reviewed_lineage','supporting_genes','genomic_evidence','reviewer_notes']],on='cluster_id').to_csv(F/'fibroblast_boundary_review.tsv',sep='\t',index=False,encoding='utf-8-sig')
print(summary.to_string());print('TOTAL',len(res),res.cells.sum(),'unknown',pd.DataFrame(detail).query('lineage_code=="U"').cells.sum(),'candidate malign',res.malignancy_assessment.value_counts().to_dict());
(O/'labels.json').write_text(json.dumps(labels,ensure_ascii=False,indent=2),encoding='utf-8')
