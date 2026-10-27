from pathlib import Path
import json,hashlib,itertools,runpy,contextlib,math
import numpy as np,pandas as pd
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from statsmodels.stats.multitest import multipletests

P=Path('C:/Users/31937/CAF_paper_figures_20261021')
O=Path('C:/Users/31937/Downloads/CAF_reference_figures_600DPI_20261021')
B=P/'batch02';B.mkdir(exist_ok=True)
for s in ['code','tables','resources','previews']:(B/s).mkdir(exist_ok=True)
genes=['PNKD','TOR1AIP1','GRID2','TAF10','BEST1','RGS2-AS1']
d=pd.read_csv(P/'Fig1i_cell_values.tsv.gz',sep='\t',encoding='utf-8')
assert len(d)==6633 and d.cell_id.is_unique and d.donor.nunique()==39
assert d.groupby('donor').site.nunique().eq(1).all()
assert hashlib.sha256((P/'Fig1i_cell_values.tsv.gz').read_bytes()).hexdigest()=='29da3ee783fda0e7cd93cb2ed5298a9d6a3c4e16f90d821419afc16028dc5a10'
manifest=[]
base=(P/'plot_GRID2_violin.py').read_text(encoding='utf-8')
for gene in ['PNKD','TOR1AIP1','TAF10','RGS2-AS1']:
 s=base.replace('GRID2',gene).replace('raw_RGS2-AS1_counts=','raw_gene_counts=')
 vmax=d.groupby(['source','class','donor'])[gene+'__log'].mean().max()
 # Nice positive tick range with adequate space above the largest observation.
 power=10**math.floor(math.log10(vmax/4));step=next(z*power for z in [1,2,2.5,5,10] if z*power>=vmax/4)
 upper=math.ceil(vmax/step)*step
 s=s.replace('ax.set_ylim(-.09,2.65)',f'ax.set_ylim({-upper*.035},{upper*1.06})').replace('np.arange(0,2.6,.5)',f'np.arange(0,{upper+step*.1},{step})')
 s=s.replace('输出：','各基因图采用覆盖自身观测范围的纵轴；比较基因时需读取刻度，不直接比较图中高度。\n输出：')
 path=B/'code'/f'plot_{gene}_violin.py';path.write_text(s,encoding='utf-8')
 with (B/f'{gene}_run.log').open('w',encoding='utf-8') as log,contextlib.redirect_stdout(log):runpy.run_path(str(path),run_name='__main__')
 name=f'Fig1h_{gene}_patient_expression_violin'
 qa=json.loads((P/(name+'_QA.json')).read_text(encoding='utf-8'))
 manifest.append(dict(name=name,kind='patient_expression',gene=gene,exports=qa['exports'],preview=str(P/(name+'_preview.png'))))
 print('EXPRESSION_COMPLETE',gene,flush=True)

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.linewidth':.8,'savefig.facecolor':'white'})
def save(fig,name,kind,gene=None):
 exports=[]
 for ext in ['png','tiff']:
  path=O/(name+'.'+ext);kw={'pil_kwargs':{'compression':'tiff_lzw'}} if ext=='tiff' else {}
  fig.savefig(path,dpi=600,**kw)
  with Image.open(path) as im:
   im.load();dpi=[float(x) for x in im.info['dpi']]
   assert max(abs(x-600) for x in dpi)<.05
   assert im.size==tuple(round(x*600) for x in fig.get_size_inches())
   exports.append(dict(file=str(path),pixels=im.size,dpi=dpi,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
 preview=B/'previews'/(name+'.png');fig.savefig(preview,dpi=130);plt.close(fig)
 manifest.append(dict(name=name,kind=kind,gene=gene,exports=exports,preview=str(preview)))

rows=[]
for (src,site,donor),q in d.groupby(['source','site','donor'],sort=True):
 for gene in genes:
  rows.append(dict(source=src,site=site,donor=donor,gene=gene,cells=len(q),detected_cells=int(q[gene+'__raw'].gt(0).sum()),RNA_detected_percent=float(q[gene+'__raw'].gt(0).mean()*100)))
patient=pd.DataFrame(rows)
assert len(patient)==39*6
patient.to_csv(B/'tables/RNA_detection_patient_values.tsv',sep='\t',index=False,encoding='utf-8')
stats=[]
for src in ['SCPCP000017','SCPCP000023']:
 for gene in genes:
  q=patient[(patient.source==src)&(patient.gene==gene)]
  b=q[q.site.eq('Bone')].RNA_detected_percent.to_numpy();l=q[q.site.eq('Lung')].RNA_detected_percent.to_numpy();x=np.r_[b,l]
  obs=float(l.mean()-b.mean());count=math.comb(len(x),len(l));tol=1e-12
  if np.ptp(x)==0:p=np.nan;status='no_variation_effect_only'
  else:
   extreme=0
   for ix in itertools.combinations(range(len(x)),len(l)):
    ls=x[list(ix)].sum();delta=ls/len(l)-(x.sum()-ls)/len(b)
    extreme+=abs(delta)>=abs(obs)-tol
   p=extreme/count;status='exact_two_sided_permutation'
  stats.append(dict(source=src,gene=gene,Bone_N=len(b),Lung_N=len(l),Bone_mean_percent=float(b.mean()),Lung_mean_percent=float(l.mean()),effect_lung_minus_bone_pp=obs,P=p,permutations=count,status=status,family=src+'|six_candidates|RNA_detection_fraction'))
st=pd.DataFrame(stats);st['FDR']=np.nan
for src,ix in st.groupby('source').groups.items():
 valid=st.loc[ix].P.dropna().index
 st.loc[valid,'FDR']=multipletests(st.loc[valid,'P'],method='fdr_bh')[1]
st.to_csv(B/'tables/RNA_detection_exact_statistics.tsv',sep='\t',index=False,encoding='utf-8',na_rep='NA')
# Independent exact-permutation check via scipy, including the p=NA rule for invariant endpoints.
from scipy.stats import permutation_test
for r in st.itertuples():
 if pd.isna(r.P):continue
 q=patient[(patient.source==r.source)&(patient.gene==r.gene)]
 b=q[q.site.eq('Bone')].RNA_detected_percent.to_numpy();l=q[q.site.eq('Lung')].RNA_detected_percent.to_numpy()
 z=permutation_test((b,l),lambda x,y:float(y.mean()-x.mean()),permutation_type='independent',n_resamples=np.inf,alternative='greater',vectorized=False)
 checked=float((abs(z.null_distribution)>=abs(r.effect_lung_minus_bone_pp)-1e-12).mean())
 assert abs(checked-r.P)<1e-12

positions=[('SCPCP000017','Bone'),('SCPCP000017','Lung'),('SCPCP000023','Bone'),('SCPCP000023','Lung'),('NCI','Lung')]
colors={'Bone':'#477DA8','Lung':'#D27A79'}
for gene in genes:
 q=patient[patient.gene.eq(gene)];fig=plt.figure(figsize=(9.6,6.6),dpi=140);ax=fig.add_axes([.10,.27,.86,.48]);rng=np.random.default_rng(20261021)
 for i,(src,site) in enumerate(positions):
  vals=q[(q.source==src)&(q.site==site)].RNA_detected_percent.to_numpy()
  ax.scatter(i+rng.uniform(-.12,.12,len(vals)),vals,c=colors[site],s=30,edgecolor='white',lw=.5,zorder=4)
  low,med,high=np.quantile(vals,[.25,.5,.75]);ax.plot([i,i],[low,high],color='#25313A',lw=3,zorder=2)
  ax.plot([i-.19,i+.19],[med,med],color='#25313A',lw=1.3,zorder=5)
  ax.text(i,-.10,f'n={len(vals)}',transform=ax.get_xaxis_transform(),ha='center',fontsize=9)
 for x,src in [(.5,'SCPCP000017'),(2.5,'SCPCP000023')]:
  r=st[(st.source==src)&(st.gene==gene)].iloc[0]
  txt=f"Lung - Bone: {r.effect_lung_minus_bone_pp:+.2f} pp"
  if src=='SCPCP000017':txt+=f'\nP={r.P:.3g}; FDR={r.FDR:.3g}' if pd.notna(r.P) else '\nNo variation; effect only'
  else:txt+='\nExploratory (7 patients)'
  ax.text(x,1.08,txt,transform=ax.get_xaxis_transform(),ha='center',fontsize=9,linespacing=1.6)
 ax.text(4,1.08,'Lung only\nNo within-source contrast',transform=ax.get_xaxis_transform(),ha='center',fontsize=9,linespacing=1.6)
 ax.set_ylim(-3.5,103);ax.set_yticks([0,25,50,75,100]);ax.set_xlim(-.5,4.5)
 ax.set_ylabel(f'{gene} RNA-detected cells (%)',fontsize=11,labelpad=10);ax.set_xticks(range(5),[s for _,s in positions])
 ax.spines[['top','right']].set_visible(False);ax.grid(axis='y',color='#E8EBEE',lw=.6);ax.tick_params(axis='x',length=0,pad=7)
 for sep in [1.5,3.5]:ax.axvline(sep,color='#BBC2C8',lw=.8)
 for x,src in [(.5,'SCPCP000017'),(2.5,'SCPCP000023'),(4,'NCI')]:ax.text(x,-.22,src,transform=ax.get_xaxis_transform(),ha='center',color='#53616B')
 fig.text(.10,.94,f'{gene} RNA detection by anatomical site',fontsize=16,fontweight='bold')
 fig.text(.10,.89,'One point per patient; original source-specific cohorts retained',fontsize=11,color='#53616B')
 fig.text(.10,.12,'Points: patients   |   horizontal lines: medians   |   vertical bars: interquartile ranges',fontsize=9)
 fig.text(.10,.065,'RNA detection: raw count >0 within the analyzed candidate pool. Not flow-cytometry positivity.',fontsize=9,color='#53616B')
 fig.text(.10,.027,'Source-specific exact permutations; BH FDR across six candidates per source. pp: percentage points.',fontsize=9,color='#53616B')
 save(fig,f'Fig1m_{gene}_patient_RNA_detection','patient_detection',gene)
 print('DETECTION_COMPLETE',gene,flush=True)

classcolors={'F1':'#477DA8','F2':'#5A9B86','F3':'#9D82BC','F4':'#D27A79','F5':'#C3A353'}
comp=[]
for src in ['SCPCP000017','SCPCP000023']:
 q=d[d.source.eq(src)];meta=q.groupby(['donor','site'],sort=True).size().rename('cells').reset_index().sort_values(['site','donor'])
 counts=pd.crosstab(q.donor,q['class']).reindex(index=meta.donor,columns=list(classcolors),fill_value=0)
 fractions=counts.div(counts.sum(axis=1),axis=0)*100;assert np.allclose(fractions.sum(1),100)
 fig=plt.figure(figsize=(10.4,6.6),dpi=140);ax=fig.add_axes([.10,.31,.86,.47]);bottom=np.zeros(len(meta))
 for cl,color in classcolors.items():
  vals=fractions[cl].to_numpy();ax.bar(range(len(meta)),vals,bottom=bottom,color=color,width=.75,label=cl,edgecolor='white',lw=.4);bottom+=vals
  for row,value in zip(meta.itertuples(),vals):comp.append(dict(source=src,donor=row.donor,site=row.site,cells=row.cells,group=cl,fraction_percent=value))
 for i,row in enumerate(meta.itertuples()):ax.text(i,102,str(row.cells),ha='center',fontsize=8)
 split=int(meta.site.eq('Bone').sum());ax.axvline(split-.5,color='#53616B',lw=.9)
 ax.text((split-1)/2,1.13,'Bone',transform=ax.get_xaxis_transform(),ha='center',fontsize=11)
 ax.text((split+len(meta)-1)/2,1.13,'Lung',transform=ax.get_xaxis_transform(),ha='center',fontsize=11)
 labels=[str(v).split('|')[-1] for v in meta.donor]
 ax.set_xticks(range(len(meta)),labels,rotation=50,ha='right',fontsize=8);ax.set_ylim(0,100);ax.set_xlim(-.6,len(meta)-.4)
 ax.set_yticks([0,25,50,75,100]);ax.set_ylabel('Share of analyzed candidate cells (%)',fontsize=11)
 ax.spines[['top','right']].set_visible(False);ax.tick_params(axis='x',length=0)
 fig.text(.10,.94,f'CAF candidate composition | {src}',fontsize=16,fontweight='bold')
 fig.text(.10,.89,'One bar per patient; labels above bars indicate candidate-cell counts',fontsize=11,color='#53616B')
 fig.legend(*ax.get_legend_handles_labels(),loc='lower left',bbox_to_anchor=(.09,.072),ncol=5,frameon=False)
 fig.text(.10,.035,'Denominator: analyzed expanded candidates, not all tissue cells. Existing F classes retained.',fontsize=9,color='#53616B')
 save(fig,f'Fig1l_{src}_patient_candidate_composition','patient_composition')
 print('COMPOSITION_COMPLETE',src,flush=True)
pd.DataFrame(comp).to_csv(B/'tables/patient_class_composition.tsv',sep='\t',index=False,encoding='utf-8')
(B/'resources/BATCH02_MANIFEST.json').write_text(json.dumps(dict(status='data_export_PASS_visual_pending',new_figures=len(manifest),figures=manifest,input_sha256=hashlib.sha256((P/'Fig1i_cell_values.tsv.gz').read_bytes()).hexdigest(),exact_permutation_independent_check='PASS',FDR_families='6 prespecified gene-detection hypotheses per original source; valid P only',original_results_modified=False),indent=2),encoding='utf-8')
assert len(manifest)==12
print('BATCH02_COMPUTATION_COMPLETE',len(manifest),flush=True)
