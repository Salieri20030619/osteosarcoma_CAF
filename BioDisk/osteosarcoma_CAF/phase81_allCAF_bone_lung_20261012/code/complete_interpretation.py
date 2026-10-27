from pathlib import Path
import json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase81_allCAF_bone_lung_20261012')
def read(n):return pd.read_csv(O/'tables'/n,sep='\t',low_memory=False)
g=read('all_genes_statistics.tsv.gz');g=g[g.sensitivity.eq('all_cells')]
p=read('all_programs_statistics.tsv.gz');p=p[p.sensitivity.eq('all_cells')]
summ=read('comparison_results_summary.tsv');testable=summ[(summ.sensitivity=='all_cells')&(summ.status=='testable')].sort_values(['layer','dataset','design'])
cross=[]
for kind,df,key in [('gene',g,'feature'),('program',p,'program')]:
    for layer in ['supported','expanded']:
        base=df[(df.layer==layer)&(df.dataset=='SCPCP000017')&(df.design=='independent')&df.nominal]
        for ds in sorted(set(df.dataset)-{'SCPCP000017'}):
            other=df[(df.layer==layer)&(df.dataset==ds)&(df.design=='independent')]
            if other.empty:continue
            m=base[[key,'mean_lung_minus_bone']].merge(other[[key,'mean_lung_minus_bone','P','FDR','bone_N','lung_N']],on=key,suffixes=('_primary','_other'))
            if m.empty:continue
            prod=m.mean_lung_minus_bone_primary*m.mean_lung_minus_bone_other
            cross.append(dict(kind=kind,layer=layer,primary_source='SCPCP000017',other_source=ds,matched_records=len(m),same_direction=int((prod>0).sum()),opposite_direction=int((prod<0).sum()),zero=int((prod==0).sum()),other_valid_tests=int(m.P.notna().sum()),other_nominal=int(m.P.lt(.05).sum()),other_FDR_hits=int(m.FDR.lt(.05).sum())))
cross=pd.DataFrame(cross);cross.to_csv(O/'tables/cross_source_nominal_signal_summary.tsv',sep='\t',index=False)
lib=p[p.nominal].copy();lib['library']=lib.program.str.split('|').str[0]
libsum=lib.groupby(['layer','dataset','design','library']).size().reset_index(name='nominal_programs');libsum.to_csv(O/'tables/program_library_summary.tsv',sep='\t',index=False)
genes_main=read('main_nominal_genes_with_sensitivity.tsv.gz');programs_main=read('main_nominal_programs_with_sensitivity.tsv.gz')
paired=[]
for kind,df in [('gene',genes_main),('program',programs_main)]:
    for ctx,z in df[df.design.eq('independent')].groupby('context_id'):
        v=z[z.paired_effect.notna()];prod=v.mean_lung_minus_bone*v.paired_effect
        paired.append(dict(kind=kind,context_id=ctx,nominal_records=len(z),paired_effect_available=len(v),paired_same=int((prod>0).sum()),paired_opposite=int((prod<0).sum()),paired_zero=int((prod==0).sum()),doublet_nominal_retained=int(z.doublet_excluded_P.lt(.05).sum())))
paired=pd.DataFrame(paired);paired.to_csv(O/'tables/paired_and_doublet_consistency_summary.tsv',sep='\t',index=False)
report=O/'delivery/全部CAF候选_骨肺病灶比较_20261012.txt'
text=report.read_text()
summary='''
主要结果
整体候选池合并后检出了名义表达差异。主要来源SCPCP000017的扩展层骨12/肺3名独立患者有1759条上调、88条下调；支持层骨10/肺2名有1674条上调、4条下调。另一来源SCPCP000023扩展层骨4/肺3名有156条上调、7条下调。扩展层2对骨肺患者未出现名义差异基因，支持层仅1对保留方向描述。
预设客户程序中，支持层成骨程序肺侧较低，平均差-0.381624，P=0.045455、FDR=0.812937，属于名义显著。全部患者级名义基因与程序均未通过相应FDR阈值；不能将这一结果表述为两种组织没有表达差异，也不能直接据此确认稳定的部位特异CAF功能。
'''
text=text.replace('2026年10月12日\n','2026年10月12日\n'+summary,1)
text+='\n配对及去双细胞一致性\n'+paired.to_string(index=False)+'\n'
text+='\n其他来源方向\n'+cross.to_string(index=False)+'\n'
text+='\n名义功能程序按库汇总\n'+libsum.to_string(index=False)+'\n'
text+='\n上一步审计结论\n原亚群分析的分组、原始计数聚合、198项独立基因P、130项程序P及FDR复核未发现计算错误；4050条比例分母核查通过。FDR中把无效检验按P=1纳入确实保守，但仅纳入有效P的诊断计算仍无FDR显著结果。本轮重新定义整体CAF假设，预先采用有效检验家族，并保留旧口径诊断，结果同样没有患者级FDR显著信号。小样本、细胞覆盖和方向不一致仍是实际证据限制。\n'
report.write_text(text,encoding='utf-8')
# Fix an ambiguous plot label: FDR hit count, never an FDR value of zero.
rp=O/'code/review_allcaf.py';code=rp.read_text();code=code.replace("; FDR {r.FDR_genes}","; FDR hits: {r.FDR_genes}")
rp.write_text(code,encoding='utf-8')
segment=code[code.index('plt.rcParams.update'):code.index('cp=customer')]
exec(compile(segment,str(rp)+'#figure1','exec'),dict(plt=plt,np=np,testable=testable,O=O))
print('PAIRED\n'+paired.to_string(index=False));print('CROSS_SOURCE\n'+cross.to_string(index=False));print('PROGRAM_LIBRARIES\n'+libsum.to_string(index=False))
