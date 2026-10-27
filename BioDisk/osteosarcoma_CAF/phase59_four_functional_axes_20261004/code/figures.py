from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
import pandas as pd,numpy as np,json,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase59_four_functional_axes_20261004')
def read(p):return pd.read_csv(p,sep='\t',low_memory=False)
assert (O/'QA_COMPLETE').exists();s=read(O/'tables/axis_statistics.tsv');ob=read(O/'tables/frozen_cells_with_identity.tsv.gz');c=read(O/'tables/lineage_RNA_patient_equal_summary.tsv');ds='SCPCP000017';axes=['inflammation','contraction','matrix','osteogenic'];classes=['F1','F2','F3','F4','F5'];plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':180})
def star(r):return '**' if r.FDR<.05 else ('*' if r.P<.05 else '')
def heat(ax,d,rows,columns,rowcol='contrast',colcol='program',title='',limit=2):
 a=d.pivot(index=rowcol,columns=colcol,values='mean_delta').reindex(index=rows,columns=columns);im=ax.imshow(a.to_numpy(),cmap='RdBu_r',vmin=-limit,vmax=limit,aspect='auto');ax.set_xticks(range(len(columns)),[v.replace('_',' ') for v in columns],rotation=25,ha='right');ax.set_yticks(range(len(rows)),[v.replace('_vs_rest','').replace('conservative_Fibro','conservative Fibro').replace('_',' ') for v in rows]);ax.set_title(title,pad=12,weight='bold')
 for i,r in enumerate(rows):
  for j,col in enumerate(columns):
   q=d[(d[rowcol]==r)&(d[colcol]==col)]
   if not len(q):continue
   x=q.iloc[0];lab='NE' if not np.isfinite(x.mean_delta) else f'{x.mean_delta:+.2f}{star(x)}\n{x.positive_patients}/{x.N}';ax.text(j,i,lab,ha='center',va='center',fontsize=9,color='white' if np.isfinite(x.mean_delta) and abs(x.mean_delta)>limit*.6 else 'black')
 return im
fig,aa=plt.subplots(1,2,figsize=(11,5.4),layout='constrained');d=s[(s.dataset==ds)&(s.context=='main')&s.program.isin(axes)]
for ax,mode in zip(aa,['raw','corrected']):im=heat(ax,d[d.matrix==mode],[cl+'_vs_rest' for cl in classes],axes,title=mode.title()+' full-library pseudobulk')
fig.colorbar(im,ax=aa,label='Mean patient log2(CP10k + 0.1) difference',shrink=.8);fig.suptitle('Four overlapping functional axes | SCPCP000017',weight='bold');fig.supxlabel('Cell labels preserved. Annotation: effect, positive patients / N. * P < 0.05; ** FDR < 0.05. NE: not evaluable.',fontsize=9);fig.savefig(O/'figures/01_four_axes.png');plt.close(fig)
fig,aa=plt.subplots(1,2,figsize=(12,5.7),layout='constrained');contexts=['main','min10','Fibro_only','conservative_Fibro']
for ax,cl in zip(aa,['F2','F4']):
 d=s[(s.dataset==ds)&(s.matrix=='corrected')&s.context.isin(contexts)&(s.contrast==cl+'_vs_rest')&s.program.isin(axes)];im=heat(ax,d,contexts,axes,rowcol='context',title=cl+' identity / cell-count sensitivity')
fig.colorbar(im,ax=aa,label='Mean patient effect',shrink=.8);fig.supxlabel('Corrected counts; same-source patient/site pairing. Fibro is expression-review support, not a confirmed CAF identity.',fontsize=9);fig.savefig(O/'figures/02_F2_F4_sensitivity.png');plt.close(fig)
colors={'Fibro':'#277DA1','Osteo':'#D1495B','Mural':'#F4A261','Other_unresolved':'#B8BDC5'}
for cl,gn,number in [('F2',['ACTA2','TAGLN','MYL9','MYH11','COL1A1','POSTN','MMP11','LOX','LRRC15'],'03'),('F4',['RUNX2','SP7','ALPL','IBSP','BGLAP','OMD','OGN','PHOSPHO1','SOSTDC1'],'04')]:
 fig,aa=plt.subplots(1,2,figsize=(12,5.5),layout='constrained')
 for ax,mode in zip(aa,['raw','corrected']):
  d=c[(c.dataset==ds)&(c.class_id==cl)&(c.matrix==mode)&c.gene.isin(gn)];bottom=np.zeros(len(gn));
  for line,col in colors.items():
   v=d[d.lineage==line].set_index('gene').mean_patient_RNA_fraction.reindex(gn).to_numpy();ax.bar(range(len(gn)),v,bottom=bottom,color=col,label=line,width=.75);bottom+=np.nan_to_num(v)
  ax.set_xticks(range(len(gn)),gn,rotation=50,ha='right');ax.set_ylim(0,1.06);ax.set_ylabel('Mean within-patient RNA contribution');ax.set_title(mode.title(),weight='bold');ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1))
  for j,gene in enumerate(gn):
   n=d.loc[d.gene.eq(gene),'N_evaluable'];ax.text(j,1.015,'N='+str(int(n.max())) if len(n) else 'NE',ha='center',fontsize=8)
 aa[1].legend(loc='lower left',bbox_to_anchor=(0,1.12),ncol=2,frameon=False);fig.suptitle(cl+' signal attribution by reviewed expression lineage | '+ds,weight='bold');fig.supxlabel('Patients equally weighted after averaging paired sites; absent lineages contribute zero. No causal or malignant-cell identity claim.',fontsize=9);fig.savefig(O/f'figures/{number}_{cl}_RNA_sources.png');plt.close(fig)
# Cross-source consistency: retain all sources, without averaging them as independent validation.
sources=sorted(s.dataset.unique());fig,aa=plt.subplots(1,3,figsize=(13,7),layout='constrained')
for ax,cl,axis in zip(aa,['F2','F2','F4'],['contraction','matrix','osteogenic']):
 d=s[(s.matrix=='corrected')&(s.context=='main')&(s.contrast==cl+'_vs_rest')&(s.program==axis)];im=heat(ax,d,sources,[axis],rowcol='dataset',title=cl+' '+axis,limit=2.5)
fig.colorbar(im,ax=aa,label='Mean patient effect',shrink=.7);fig.supxlabel('NE retained; source-local donor counts. N < 10 exploratory. Source records are not assumed globally independent.',fontsize=9);fig.savefig(O/'figures/05_source_consistency.png');plt.close(fig)
fig,aa=plt.subplots(1,2,figsize=(12,7),layout='constrained')
for ax,cl in zip(aa,['F2','F4']):
 rows=[cl+'_'+line+'_vs_BG_'+bg for line in ['All','Fibro'] for bg in ['Fibro','Mural','Osteo']];d=s[(s.dataset==ds)&(s.context=='external_'+cl)&s.program.isin(axes)];im=heat(ax,d,rows,axes,title=cl+' versus external reviewed backgrounds',limit=3)
fig.colorbar(im,ax=aa,label='Raw-count patient effect',shrink=.7);fig.supxlabel('Background excludes original 8,174 cells; same-source patient/site. BG labels describe expression lineage. * P < 0.05; ** FDR < 0.05.',fontsize=9);fig.savefig(O/'figures/06_external_competitors.png');plt.close(fig)
summary={'main':s[(s.dataset==ds)&(s.context=='main')&s.program.isin(axes)].to_dict('records'),'focus':s[(s.dataset==ds)&(s.matrix=='corrected')&s.context.isin(['Fibro_only','conservative_Fibro','within_F2','within_F4','min10'])&s.program.isin(axes+['contractile_markers','matrix_markers','osteogenic_effectors','Liu_C2_OMD_OGN','Liu_C5_ossification','mineralization'])].to_dict('records'),'external':s[(s.dataset==ds)&s.context.str.startswith('external')&s.program.isin(axes)].to_dict('records'),'contributions':c[(c.dataset==ds)&(c.matrix=='corrected')&c.gene.isin(['ACTA2','MYL9','COL1A1','POSTN','RUNX2','ALPL','SP7','IBSP','OMD','OGN'])].to_dict('records'),'identity':ob.groupby(['class_id','lineage_stratum']).size().reset_index(name='cells').to_dict('records')}
(O/'delivery/SUMMARY.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2));print('FIGURES_COMPLETE',flush=True)
