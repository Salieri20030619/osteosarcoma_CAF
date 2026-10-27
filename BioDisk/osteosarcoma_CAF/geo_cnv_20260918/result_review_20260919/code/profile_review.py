from pathlib import Path
import pandas as pd,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/geo_cnv_20260918');O=P/'result_review_20260919';T=O/'tables'
z=pd.read_csv(T/'cell_review.tsv.gz',sep='\t');su=pd.read_csv(T/'sample_summary.tsv',sep='\t');pro=[];cor=[];cov=[]
def corr(x,y):
 if len(x)<3 or np.std(x)<1e-10 or np.std(y)<1e-10:return np.nan
 return float(np.corrcoef(x,y)[0,1])
for sm,g in z.groupby('sample'):
 a=pd.read_csv(P/'infercnv'/sm/'chromosome_log2_medians.tsv',sep='\t').set_index('chromosome');a.index=a.index.astype(str).str.replace('chr','',regex=False)
 b=pd.read_csv(T/(sm+'_copy_chromosomes.tsv.gz'),sep='\t').set_index('chromosome');b.index=b.index.astype(str)
 ch=[str(x) for x in range(1,23)];a=a.loc[ch];b=b.loc[ch]
 train=g.loc[g.role_group=='training','tool_id'];base_a=a[a.columns.intersection(train)].median(axis=1);base_b=b[b.columns.intersection(train)].median(axis=1)
 a=a.subtract(base_a,axis=0);b=b.subtract(base_b,axis=0)
 common=a.columns.intersection(b.columns);ga=g.set_index('tool_id').loc[common]
 av=a[common].to_numpy();bv=b[common].to_numpy();ax=av-av.mean(axis=0);bx=bv-bv.mean(axis=0);den=np.sqrt((ax*ax).sum(0)*(bx*bx).sum(0));rr=np.divide((ax*bx).sum(0),den,out=np.full(len(common),np.nan),where=den>1e-12)
 ga=ga[['cell_id','sample','lineage_phase3','joint','cnv_role']].copy();ga['chromosome_profile_r']=rr;cor.append(ga.reset_index())
 for (lin,j),gg in g.groupby(['lineage_phase3','joint']):
  ids=pd.Index(gg.tool_id).intersection(common)
  if not len(ids):continue
  ma=a[ids].median(axis=1);mb=b[ids].median(axis=1)
  for cc in ch:pro.append(dict(sample=sm,lineage=lin,joint=j,n_profiled=len(ids),chromosome=int(cc),infer_centered_median=ma[cc],copy_centered_median=mb[cc],group_profile_r=corr(ma.to_numpy(),mb.to_numpy())))
 go=pd.read_csv(P/'infercnv'/sm/'gene_positions_used.tsv',sep='\t');chrcol='chr' if 'chr' in go else 'chromosome'
 cc=go[chrcol].astype(str).str.replace('chr','',regex=False)
 for c in ch:
  part=go[cc==c];pos=part.start.sort_values().to_numpy();cov.append(dict(sample=sm,chromosome=int(c),infer_genes=len(part),max_intergene_gap_bp=float(np.diff(pos).max()) if len(pos)>1 else np.nan))
 print('PROFILE_REVIEWED',sm,flush=True)
pr=pd.DataFrame(pro);cr=pd.concat(cor,ignore_index=True);pr.to_csv(T/'group_chromosome_profiles.tsv',sep='\t',index=False);cr.to_csv(T/'cell_profile_concordance.tsv.gz',sep='\t',index=False);pd.DataFrame(cov).to_csv(T/'infer_postfilter_coverage.tsv',sep='\t',index=False)
cs=cr.groupby(['sample','lineage_phase3','joint']).agg(profiled_cells=('cell_id','size'),finite_correlations=('chromosome_profile_r','count'),median_profile_r=('chromosome_profile_r','median')).reset_index();cs.to_csv(T/'profile_concordance_summary.tsv',sep='\t',index=False)
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
keys=['both_low','both_high','infer_low_copy_high','infer_high_copy_low','copy_not_evaluable','infer_ambiguous']
labels=['Both low','Both high','infer low / CopyKAT high','infer high / CopyKAT low','CopyKAT not evaluable','infer ambiguous']
colors=['#3b9b86','#b64a55','#e3a647','#856ab5','#adb5bd','#436b9e']
fig,axes=plt.subplots(1,2,figsize=(15,8),gridspec_kw={'width_ratios':[1.3,1]})
s=su[su.fib_cells>0].copy().sort_values(['dataset','sample']);yy=np.arange(len(s));left=np.zeros(len(s))
for k,l,c in zip(keys,labels,colors):
 v=s['fib_'+k].to_numpy()/s.fib_cells.to_numpy()*100;axes[0].barh(yy,v,left=left,label=l,color=c);left+=v
axes[0].set_yticks(yy,[a.split('__')[1]+' (n='+str(n)+')' for a,n in zip(s['sample'],s.fib_cells)]);axes[0].invert_yaxis();axes[0].set_xlim(0,100);axes[0].set_xlabel('Percent of FibLike cells (all cells retained in denominator)');axes[0].set_title('FibLike candidates: paired CNV calls');axes[0].legend(loc='upper left',bbox_to_anchor=(0,-.12),ncol=2,fontsize=9)
s=su.sort_values(['dataset','sample']);yy=np.arange(len(s));axes[1].barh(yy,s.holdout_copy_high/s.holdout*100,label='CopyKAT high',color='#b64a55');axes[1].barh(yy,s.holdout_copy_ne/s.holdout*100,left=s.holdout_copy_high/s.holdout*100,label='CopyKAT not evaluable',color='#adb5bd');axes[1].set_yticks(yy,[a.split('__')[1]+' (n='+str(n)+')' for a,n in zip(s['sample'],s.holdout)]);axes[1].invert_yaxis();axes[1].set_xlim(0,100);axes[1].set_xlabel('Percent of heldout normal candidates');axes[1].set_title('Normal calibration: high and missing calls');axes[1].legend(loc='upper left',bbox_to_anchor=(0,-.12),fontsize=9)
fig.suptitle('GEO CNV result review | descriptive cell counts, not biological replicates',fontsize=15);fig.tight_layout(rect=(0,.1,1,.95));fig.savefig(O/'figures/cnv_review_overview.png',dpi=180);plt.close(fig)
fig,axes=plt.subplots(2,2,figsize=(15,8),sharex=True)
for row,sm in enumerate(['GSE318840__GSM9504303','GSE318840__GSM9504304']):
 for col,(m,lab) in enumerate([('infer_centered_median','inferCNV log2 expression deviation'),('copy_centered_median','CopyKAT relative CNA signal')]):
  ax=axes[row,col]
  for k,l,c in zip(keys,labels,colors):
   q=pr[(pr['sample']==sm)&(pr.lineage=='FibLike')&(pr.joint==k)].sort_values('chromosome')
   if len(q) and q.n_profiled.iloc[0]>=5:ax.plot(q.chromosome,q[m],marker='.',lw=1.5,color=c,label=l+' (n='+str(q.n_profiled.iloc[0])+')')
  ax.axhline(0,c='#555555',lw=.6);ax.set_title(sm.split('__')[1]+' | '+lab);ax.set_ylabel('Centered signal (method-specific units)');ax.set_xticks(range(1,23,2));ax.legend(fontsize=8,loc='best');ax.set_xlabel('Autosome')
fig.suptitle('Lung FibLike chromosome profiles | method scales differ',fontsize=15);fig.tight_layout(rect=(0,0,1,.95));fig.savefig(O/'figures/lung_chromosome_profiles.png',dpi=180);plt.close(fig)
print('LUNG PROFILES');print(cs[cs['sample'].str.endswith(('GSM9504303','GSM9504304'))&(cs.lineage_phase3=='FibLike')].to_string(index=False))
print('LUNG GROUP CORR');print(pr[pr['sample'].str.endswith(('GSM9504303','GSM9504304'))&(pr.lineage=='FibLike')][['sample','joint','n_profiled','group_profile_r']].drop_duplicates().to_string(index=False))
print('POSTFILTER_MIN_GENES');print(pd.DataFrame(cov).groupby('sample').infer_genes.min().to_string())
print('PROFILE_REVIEW_COMPLETE')
