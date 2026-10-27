from pathlib import Path
import json,hashlib,math
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.cbook import boxplot_stats
from PIL import Image
P=Path('C:/Users/31937/CAF_paper_figures_20261021');B=P/'batch02'
R=B/'fig1m_boxplot_v2'
O=Path('C:/Users/31937/Downloads/CAF_reference_figures_600DPI_20261021/Fig1m_boxplot_v2')
for p in [R,O]:p.mkdir(parents=True,exist_ok=True)
patient=pd.read_csv(B/'tables/RNA_detection_patient_values.tsv',sep='\t');st=pd.read_csv(B/'tables/RNA_detection_exact_statistics.tsv',sep='\t')
hashes={f:hashlib.sha256((B/'tables'/f).read_bytes()).hexdigest() for f in ['RNA_detection_patient_values.tsv','RNA_detection_exact_statistics.tsv']}
genes=['PNKD','TOR1AIP1','GRID2','TAF10','BEST1','RGS2-AS1'];positions=[('SCPCP000017','Bone'),('SCPCP000017','Lung'),('SCPCP000023','Bone'),('SCPCP000023','Lung'),('NCI','Lung')]
colors={'Bone':'#477DA8','Lung':'#D27A79'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.linewidth':.8,'savefig.facecolor':'white'})
manifest=[];summaries=[]
for gene in genes:
 q=patient[patient.gene.eq(gene)];fig=plt.figure(figsize=(9.6,6.6),dpi=140);ax=fig.add_axes([.10,.27,.86,.48]);rng=np.random.default_rng(20261021)
 upper=max(10,math.ceil(float(q.RNA_detected_percent.max())*1.12/5)*5)
 for i,(src,site) in enumerate(positions):
  vals=q[(q.source==src)&(q.site==site)].RNA_detected_percent.to_numpy()
  bs=boxplot_stats(vals,whis=1.5,autorange=False)[0]
  q1,med,q3=np.quantile(vals,[.25,.5,.75]);iqr=q3-q1
  assert np.allclose([bs['q1'],bs['med'],bs['q3']],[q1,med,q3])
  lo=min(q1,vals[vals>=q1-1.5*iqr].min());hi=max(q3,vals[vals<=q3+1.5*iqr].max())
  assert np.allclose([bs['whislo'],bs['whishi']],[lo,hi])
  artists=ax.bxp([bs],positions=[i],widths=.50,patch_artist=True,showfliers=False,manage_ticks=False,boxprops=dict(facecolor=colors[site],edgecolor=colors[site],alpha=.35,linewidth=1.2),medianprops=dict(color='#25313A',linewidth=1.4),whiskerprops=dict(color='#53616B',linewidth=1),capprops=dict(color='#53616B',linewidth=1),zorder=2)
  assert len(artists['boxes'])==1 and len(artists['whiskers'])==2 and len(artists['caps'])==2
  ax.scatter(i+rng.uniform(-.13,.13,len(vals)),vals,c=colors[site],s=25,edgecolor='white',lw=.5,zorder=4)
  summaries.append(dict(gene=gene,source=src,site=site,n=len(vals),q1=q1,median=med,q3=q3,whisker_low=lo,whisker_high=hi,outliers=len(bs['fliers']),box_zero_height=bool(q1==q3)))
  ax.text(i,-.10,f'n={len(vals)}',transform=ax.get_xaxis_transform(),ha='center',fontsize=9)
 for x,src in [(.5,'SCPCP000017'),(2.5,'SCPCP000023')]:
  row=st[(st.source==src)&(st.gene==gene)].iloc[0]
  txt=f"Lung - Bone: {row.effect_lung_minus_bone_pp:+.2f} pp"
  if src=='SCPCP000017':txt+=f'\nP={row.P:.3g}; FDR={row.FDR:.3g}' if pd.notna(row.P) else '\nNo variation; effect only'
  else:txt+='\nExploratory (7 patients)'
  ax.text(x,1.08,txt,transform=ax.get_xaxis_transform(),ha='center',fontsize=9,linespacing=1.6)
 ax.text(4,1.08,'Lung only\nNo within-source contrast',transform=ax.get_xaxis_transform(),ha='center',fontsize=9,linespacing=1.6)
 ax.set_ylim(-upper*.035,upper);ax.set_xlim(-.5,4.5)
 step=5 if upper<=35 else (10 if upper<=60 else 20)
 ax.set_yticks(np.arange(0,upper+.1,step));ax.set_ylabel(f'{gene} RNA-detected cells (%)',fontsize=11,labelpad=10)
 ax.set_xticks(range(5),[s for _,s in positions]);ax.spines[['top','right']].set_visible(False);ax.grid(axis='y',color='#E8EBEE',lw=.6);ax.tick_params(axis='x',length=0,pad=7)
 for sep in [1.5,3.5]:ax.axvline(sep,color='#BBC2C8',lw=.8)
 for x,src in [(.5,'SCPCP000017'),(2.5,'SCPCP000023'),(4,'NCI')]:ax.text(x,-.22,src,transform=ax.get_xaxis_transform(),ha='center',color='#53616B')
 fig.text(.10,.94,f'{gene} RNA detection by anatomical site',fontsize=16,fontweight='bold')
 fig.text(.10,.89,'Patient-level boxplots with individual observations',fontsize=11,color='#53616B')
 fig.text(.10,.12,'Box: Q1-Q3   |   line: median   |   whiskers: 1.5 x IQR rule   |   points: all patients',fontsize=9)
 fig.text(.10,.065,'RNA detection: raw count >0 within the analyzed candidate pool. Not flow-cytometry positivity.',fontsize=9,color='#53616B')
 fig.text(.10,.027,'Source-specific exact permutations; BH FDR across six candidates per source. pp: percentage points.',fontsize=9,color='#53616B')
 name=f'Fig1m_{gene}_patient_RNA_detection_boxplot_v2';exports=[]
 for ext in ['png','tiff']:
  path=O/(name+'.'+ext);fig.savefig(path,dpi=600,**({'pil_kwargs':{'compression':'tiff_lzw'}} if ext=='tiff' else {}))
  with Image.open(path) as im:
   im.load();assert im.size==(5760,3960);dpi=list(map(float,im.info['dpi']));assert max(abs(x-600) for x in dpi)<.05
   exports.append(dict(file=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),pixels=im.size,dpi=dpi))
 fig.savefig(R/(name+'_preview.png'),dpi=130);plt.close(fig)
 old=(O.parent/f'Fig1m_{gene}_patient_RNA_detection_图注与口径.txt').read_text(encoding='utf-8')
 old=old.replace('深色横线为患者比例中位数，竖线为四分位距。','箱体表示第一至第三四分位数，箱内横线为中位数；须延伸到距箱体不超过1.5倍四分位距范围内的最远实际观测值。全部患者以散点叠加，离群患者也保留，不重复绘制额外离群符号。全零或Q1=Q3时箱体自然收缩为横线，不人为制造高度。')
 old=old.replace('纵轴统一为0–100%，图中保留所有患者与零值。','各基因纵轴按观察值设置；同一基因各来源共用刻度，跨基因比较需读刻度。图中保留所有患者与零值。')
 (O/(name+'_图注.txt')).write_text('修订v2：补齐真实箱体、须和端帽；原始数据、患者名单及P/FDR不变。旧版文件保留。分位数采用线性插值；若箱缘外无满足1.5×IQR规则的观测，须止于箱缘，不向箱内反折。\n\n'+old,encoding='utf-8')
 manifest.append(dict(gene=gene,name=name,exports=exports))
pd.DataFrame(summaries).to_csv(R/'boxplot_summary.tsv',sep='\t',index=False,encoding='utf-8')
for f,h in hashes.items():assert hashlib.sha256((B/'tables'/f).read_bytes()).hexdigest()==h
(R/'QA.json').write_text(json.dumps(dict(status='data_export_pass_visual_pending',input_hashes=hashes,groups_checked=30,box_quartiles_and_whiskers_checked=True,statistics_unchanged=True,old_versions_preserved=True,figures=manifest),indent=2),encoding='utf-8')
print('Six boxplots rendered; 30 groups verified; input/statistic hashes unchanged')
