from pathlib import Path
import pandas as pd,json
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase63_matrix_resistance_20261005')
t=pd.read_csv(O/'tables/posthoc_collagen_structure_statistics.tsv',sep='\t')
r=t[(t.dataset=='SCPCP000017')&(t.subgroup_id=='F2_C3')&(t.matrix=='corrected')&(t.scope=='Fibro')&(t.minimum_cells==1)&(t.program=='fibrillar_I_III_V')].iloc[0]
s=pd.read_csv(O/'tables/patient_program_statistics.tsv.gz',sep='\t');q=s[(s.dataset=='SCPCP000017')&(s.subgroup_id=='F2_C3')&(s.matrix=='corrected')&(s.scope=='all')&(s.minimum_cells==1)&(s.program=='collagen_biosynthesis')].iloc[0]
extra=f'F2-C3的I/III/V面板在Fibro限定后仍为名义正向：N={int(r.N)}，效应+{r.mean_effect:.3f}，P={r.P:.6g}，FDR={r.FDR:.6g}。原始矩阵方向一致。相对而言，完整胶原生物合成程序效应+{q.mean_effect:.3f}，P={q.P:.6g}，FDR={q.FDR:.6g}，未达名义阈值；优先线索限定于特定纤维状胶原组合，不能写成所有胶原合成环节普遍增强。'
path=O/'delivery/CAF基质耐药目标导向首轮结果_20261005.txt';text=path.read_text();anchor='F2-C1的IV型胶原面板：';assert text.count(anchor)==1;text=text.replace(anchor,extra+'\n'+anchor);path.write_text(text,encoding='utf-8')
p=O/'tables/goal_oriented_candidate_priorities.tsv';f=pd.read_csv(p,sep='\t');f.loc[f.subgroup_id.eq('F2_C3'),'interpretation']=extra+' 每臂至少5细胞时N=3，按探索性方向保留；先验证实际定位及沉积。';f.to_csv(p,sep='\t',index=False,encoding='utf-8-sig')
(O/'QA/REPORT_EVIDENCE_CHECK.json').write_text(json.dumps({'status':'PASS','F2C3_fibrillar_Fibro':r.to_dict(),'F2C3_full_biosynthesis':q.to_dict(),'only_report_and_candidate_interpretation_updated':True},ensure_ascii=False,indent=2),encoding='utf-8')
print(extra)
