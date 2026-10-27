from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,itertools,zipfile
from scipy import sparse
from scipy.stats import false_discovery_control
R=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');I=R/'phase119_Bone_Lung_functional_subgroups_20261024';A=R/'phase120_patient_supported_CAF_review_20261026';O=R/'phase121_CAF_illustrated_report_20261026';V=O/'data_validation';V.mkdir(exist_ok=True)
def rd(p):return pd.read_csv(p,sep='\t',low_memory=False)
m=rd(I/'tables/cohort_registry.tsv.gz');l=rd(I/'tables/hierarchical_cell_labels.tsv.gz');coords=rd(O/'tables/umap_coordinates.tsv.gz');sp=rd(A/'tables/node_donor_support.tsv');d=rd(A/'tables/functional_audit.tsv');e=rd(A/'tables/patient_effects.tsv.gz');ident=rd(A/'tables/cell_identity_audit.tsv.gz');pr=rd(O/'tables/marker_dotplot_values.tsv');dec=rd(A/'tables/node_decisions_ZH.tsv');nodes=rd(I/'tables/hierarchical_nodes.tsv');meta=json.loads((O/'FINAL_STATUS.json').read_text());sig=meta['retained_signal_groups'];lit=json.loads((O/'resources/LITERATURE_CROSSCHECK.json').read_text())
assert len(m)==6633 and m.cell_id.equals(l.cell_id) and m.cell_id.equals(coords.cell_id);fit=~m.doublet
cohort=m.groupby(['site','source']).agg(all_cells=('cell_id','size'),doublet_positive=('doublet','sum'),source_local_donors=('donor','nunique')).reset_index();cohort['main_cells']=cohort.all_cells-cohort.doublet_positive
ends=['program__collagen_biosynthesis','program__collagen_crosslinking','program__fibrillar_I_III_V','program__mCAF_matrix','audit__FNDC1_paper','audit__MYH11_paper','audit__CD36_paper','audit__C7_paper','audit__collagen_processing_nonCOL','audit__osteogenic_identity']
modes=['main_min1','main_min5','main_min10','identity_min5'];ds=d[d.node_id.isin(sig)&d.endpoint.isin(ends)&d['mode'].isin(modes)].copy();es=e[e.node_id.isin(sig)&e.endpoint.isin(ends)&e['mode'].isin(modes)].copy()
# Preserve all P values in corresponding BH families; never recompute FDR over just retained candidates.
famcols=['site','source','level','mode','family'];fams=ds[famcols].drop_duplicates();family=d.merge(fams,on=famcols,how='inner');family=family[family.P.notna()].copy();family['family_id']=family[famcols].astype(str).agg(' / '.join,axis=1)
family['FDR_recomputed']=np.nan
for _,ix in family.groupby(famcols).groups.items():family.loc[ix,'FDR_recomputed']=false_discovery_control(family.loc[ix,'P'].to_numpy())
family['FDR_abs_error']=abs(family.FDR-family.FDR_recomputed);assert family.FDR_abs_error.max()<1e-12
checks=[]
for row in ds.itertuples():
 z=es[es.node_id.eq(row.node_id)&es.source.eq(row.source)&es['mode'].eq(row.mode)&es.endpoint.eq(row.endpoint)];v=z.effect.to_numpy();n=len(v);assert n==row.N
 mu=v.mean() if n else np.nan;pos=int((v>1e-12).sum());assert pos==row.positive_patients
 if n:assert abs(mu-row.mean_effect)<1e-12
 p=np.nan
 if n>=2 and n<=16 and np.any(abs(v)>1e-12):
  signs=1-2*((np.arange(2**n,dtype=np.uint32)[:,None]>>np.arange(n))&1).astype(np.int8);p=float((abs(signs@v)>=abs(v.sum())-1e-10).mean());assert abs(p-row.P)<1e-12
 checks.append(dict(node_id=row.node_id,source=row.source,mode=row.mode,endpoint=row.endpoint,N_saved=row.N,N_recomputed=n,mean_saved=row.mean_effect,mean_recomputed=mu,positive_saved=row.positive_patients,positive_recomputed=pos,P_saved=row.P,P_recomputed=p,check='一致' if n<=16 else '均值与人数一致；P为原蒙特卡洛结果'))
check=pd.DataFrame(checks)
# Recompute dotplots from saved single-cell expression, using exactly final donor-equal aggregation.
genes=rd(I/'tables/genes.tsv').symbol;gi={g:i for i,g in enumerate(genes)};markers=pr.gene.unique().tolist();X=sparse.load_npz(I/'objects/log1p_cp10k.npz')[:,[gi[g] for g in markers]].toarray();shown=sorted(set(nodes.loc[nodes.level.eq('parent'),'node_id'])|set(sig));pats=[]
for node in shown:
 target=l['leaf' if '_s' in node else 'parent'].eq(node);z=m[fit&target]
 for (source,donor),g in z.groupby(['source','donor']):
  a=X[g.index]
  for j,gene in enumerate(markers):pats.append(dict(node_id=node,source=source,donor=donor,gene=gene,cells=len(g),mean_log1pCP10k=float(a[:,j].mean()),detection_fraction=float((a[:,j]>0).mean())))
pm=pd.DataFrame(pats);agg=pm.groupby(['node_id','gene'])[['mean_log1pCP10k','detection_fraction']].mean().reset_index();a=agg.merge(pr,on=['node_id','gene']);err=float(max(abs(a.mean_log1pCP10k-a['mean']).max(),abs(a.detection_fraction-a.detection).max()));assert err<1e-7
plots=[]
for label,groups in [('图3 骨父群',nodes.loc[nodes.level.eq('parent')&nodes.site.eq('Bone'),'node_id'].tolist()),('图4 肺父群',nodes.loc[nodes.level.eq('parent')&nodes.site.eq('Lung'),'node_id'].tolist()),('图5 13信号群',sig)]:
 a=pr[pr.node_id.isin(groups)].copy();a['figure']=label;a['colour_z']=a.groupby('gene')['mean'].transform(lambda x:(x-x.mean())/(x.std(ddof=0) or 1));a['bubble_size']=a.detection*130;plots.append(a)
dot=pd.concat(plots);supp=dec[dec['父群或子群'].isin(sig)].copy();supp=supp.drop(columns=['结论','未达标原因']);supp['报告推荐']='保留探索线索';supp.loc[supp['父群或子群'].eq('Lung_04'),'报告推荐']='优先推荐'
supportchecks=[]
for row in supp.itertuples(index=False):
 node=row[0];src=row[4];q=ds[ds.node_id.eq(node)&ds.source.eq(src)&ds.endpoint.eq('program__mCAF_matrix')].set_index('mode');assert int(q.loc['main_min5','N'])==int(row[6]);assert int(q.loc['identity_min5','N'])==int(row[8]);supportchecks.append([node,src,int(row[6]),int(q.loc['main_min5','N']),int(row[8]),int(q.loc['identity_min5','N']),'一致'])
cell=m[['cell_id','source','donor','sample','site','doublet']].copy();cell['parent']=l.parent;cell['leaf']=l.leaf;cell['UMAP1']=coords.UMAP1;cell['UMAP2']=coords.UMAP2;cell['identity_sensitivity']=ident.identity_sensitivity
assert cell.loc[fit,'UMAP1'].notna().all() and cell.loc[~fit,'UMAP1'].isna().all()
attr=rd(A/'tables/Lung_04_identity_attrition.tsv')
sources=[I/'tables/cohort_registry.tsv.gz',I/'tables/hierarchical_cell_labels.tsv.gz',I/'objects/log1p_cp10k.npz',A/'tables/functional_audit.tsv',A/'tables/patient_effects.tsv.gz',O/'tables/marker_dotplot_values.tsv',O/'tables/umap_coordinates.tsv.gz',O/'delivery/CAF_illustrated_report_legends.pdf']
prov=pd.DataFrame([dict(source_file=str(p.relative_to(R)),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sources])
guide=pd.DataFrame([
 ['适用报告','骨肺CAF亚群_图文分析与Lung04推荐报告_图例补充版.pdf','19页，2026-10-26'],
 ['校验范围','数据转录、计数、图件数值和统计计算','不把数值一致等同于生物学机制验证'],
 ['第1、2页','队列计数、信号群汇总','完整6633；主分析6423；210双细胞标记'],
 ['第3、4页 图1、2','细胞UMAP','按site、parent、leaf筛选；210例未进入嵌入'],
 ['第5至7页 图3至5','marker患者层、marker图值','患者均值/检出比例等权，再按图内每基因z标准化；颜色范围截断到[-2,2]'],
 ['第8、9页 图6、7','群患者计数','细胞数；行内除以父群总数得到贡献比例'],
 ['第10页 图8','信号群汇总、人数复核','min1、min5、身份min5对照条件不同'],
 ['第11、12页 图9、10','功能统计、患者效应、统计复核','点为患者效应，菱形为等权均值'],
 ['第13页 图11','Lung04配对损失','目标保留；两位患者合格对照分别少于5'],
 ['第14页','功能统计','筛选Lung_04、NCI、main_min5或identity_min5'],
 ['第15页','信号群汇总','全部13组保留；不把父群子群相加'],
 ['第17、18页','文献证据','原始论文DOI及功能对应边界'],
 ['第19页','人数复核、统计复核','主分析7位；身份层5位；Lung_04优先推荐'],
 ['缺失值','空白表示不可评估，不表示0','N<2或全零不产生P；支持人数0是真实计数'],
 ['patient / donor','来源内患者记录','同患者多个sample不计多个N；跨来源未直接累加'],
 ['效应定义','患者目标均值减同患者对照均值','parent对同病灶其余候选；child对同父群其余细胞'],
 ['精确P复核','2^N种符号翻转，双侧，容差1e-10','统计复核表从患者效应重新计算N、均值、同向数及N<=16精确P'],
 ['FDR复核','同病灶/来源/层级/分析层/终点家族的BH校正','FDR家族表保留完整有效P，不能只对13群重算'],
 ['数据类型','本表是固定报告的校验快照','所有数值保留完整精度；报告显示为四舍五入；不含原始全基因count矩阵'],
 ['独立复核结果',f'{len(check)}条统计记录核对一致；marker最大误差{err:.3g}',f'完整BH家族最大误差{family.FDR_abs_error.max():.3g}'],
 ],columns=['报告位置或项目','查阅工作表或定义','说明'])
frames={'说明与索引':guide,'队列计数':cohort,'信号群汇总':supp,'人数复核':pd.DataFrame(supportchecks,columns=['node_id','source','报告min5_N','重算min5_N','报告身份N','重算身份N','结果']),'细胞UMAP':cell,'群患者计数':sp,'功能统计':ds,'患者效应':es,'统计复核':check,'FDR家族':family[['family_id','node_id','endpoint','N','P','FDR','FDR_recomputed','FDR_abs_error']],'marker患者层':pm,'marker图值':dot,'Lung04配对损失':attr,'文献证据':pd.DataFrame(lit),'来源记录':prov}
payload=[]
for name,df in frames.items():payload.append(dict(name=name,columns=df.columns.tolist(),rows=json.loads(df.to_json(orient='values',force_ascii=False,double_precision=15))))
(V/'workbook_data.json').write_text(json.dumps(payload,ensure_ascii=False,separators=(',',':')),encoding='utf-8')
receipt=dict(status='PASS_NUMERICAL_RECONCILIATION',sheets={k:len(v) for k,v in frames.items()},marker_max_error=err,FDR_max_error=float(family.FDR_abs_error.max()),recomputed_statistics=len(check),report_sha256=prov.iloc[-1].sha256)
(V/'NUMERICAL_QA.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')
with zipfile.ZipFile(V/'validation_payload.zip','w',zipfile.ZIP_DEFLATED) as z:
 for n in ['workbook_data.json','NUMERICAL_QA.json']:z.write(V/n,n)
print(json.dumps(receipt))
