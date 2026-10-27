#!/usr/bin/env python3
"""Evidence-gated review of fixed sample clusters and consistent plot colors."""
import argparse,hashlib,pathlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
P={'T_NK':'#165DB1','B_Plasma':'#8545AB','Myeloid':'#C15B00','OC_like':'#964329','Endo':'#007D78','Epi':'#C8283B','Mural':'#C6288B','FibLike':'#238329','OsteoLike':'#8D7500','MesAmb':'#495B9C','Mixed':'#555555','Unresolved':'#999999','Mast':'#C45C7E','DC':'#498500','Cycling':'#7A3567'}
IMM=['CD2','CD3D','CD3E','CD247','BCL11B','THEMIS','NKG7','KLRD1','GNLY','CD79A','MS4A1','MZB1','POU2AF1','LYZ','LST1','MS4A7','CD163','CSF1R','MERTK','LILRA4','CLEC4C','FLT3','FCER1A']
END=['PECAM1','VWF','CLDN5','EMCN','KDR','TEK','EGFL7']
p=argparse.ArgumentParser();p.add_argument('--results',type=pathlib.Path,required=True);args=p.parse_args()
code=pathlib.Path(__file__).resolve().parent
mapping=pd.read_csv(code/'review_map.tsv',sep='\t',dtype={'cluster':str})
s=pd.read_csv(args.results/'sample_summary.tsv',sep='\t');flags=[];audit=[];identity_rows=[];cnvrows=[]
for r in s.to_dict('records'):
 out=args.results/r['sample'];a=sc.read_h5ad(out/'annotated.h5ad')
 mk=pd.read_csv(out/'cluster_markers.tsv',sep='\t',dtype={'group':str})
 a.obs['lineage_reviewed']=a.obs.lineage_prelim.astype(str)
 for m in mapping[mapping['sample']==r['sample']].to_dict('records'):
  expected=set(m['evidence_genes'].split(';'))
  hits=expected & set(mk[mk.group==m['cluster']].symbol.astype(str))
  valid=len(hits)>=2
  if valid:a.obs.loc[a.obs.cluster.astype(str)==m['cluster'],'lineage_reviewed']=m['reviewed_label']
  audit.append({**m,'matched_markers':';'.join(sorted(hits)),'review_applied':valid})
 a.obs['lineage_reviewed']=a.obs.lineage_reviewed.astype('category')
 def ndet(genes):
  ix=np.flatnonzero(a.var.symbol.isin(genes));return np.asarray((a.layers['counts'][:,ix]>0).sum(axis=1)).ravel()
 mes=(a.obs.det_FibLike>=2)|(a.obs.det_OsteoLike>=2)|(a.obs.det_Mural>=2)
 a.obs['cnv_ref_immune_reviewed']=a.obs.lineage_reviewed.isin(['T_NK','B_Plasma','Myeloid','DC'])&(ndet(IMM)>=2)&~mes&(a.obs.det_OC<2)
 a.obs['cnv_ref_endo_reviewed']=(a.obs.lineage_reviewed=='Endo')&(ndet(END)>=2)&~mes
 a.obs['cnv_target_reviewed']=a.obs.lineage_reviewed.isin(['FibLike','OsteoLike','Mural','MesAmb','Mixed','Unresolved','Cycling'])|mes
 a.obs.to_csv(out/'cell_annotations_reviewed.tsv.gz',sep='\t',index_label='cell_id')
 a.uns['review_sha256']=hashlib.sha256((pathlib.Path(__file__).read_text()+(code/'review_map.tsv').read_text()).encode()).hexdigest()
 fig,axes=plt.subplots(1,2,figsize=(13,5))
 sc.pl.umap(a,color='cluster',ax=axes[0],show=False,title=r['short']+' clust',size=8,legend_loc='on data')
 sc.pl.umap(a,color='lineage_reviewed',ax=axes[1],show=False,title=r['short']+' lineage',size=8,palette=P)
 fig.tight_layout();fig.savefig(out/'UMAP.png',dpi=180,bbox_inches='tight');plt.close(fig)
 # Rebuild dot plot labels from the reviewed identities; expression evidence is unchanged.
 e=pd.read_csv(out/'marker_evidence.tsv',sep='\t',dtype={'cluster':str})
 dot=['CD3D','NKG7','CD79A','MZB1','LYZ','LST1','ACP5','CTSK','PECAM1','VWF','EPCAM','SFTPC','RGS5','NOTCH3','DCN','LUM','PDGFRA','FAP','RUNX2','ALPL','IBSP','ACTA2','TAGLN']
 extra=['MS4A7','CD163','MERTK','BCL11B','THEMIS','KLRD1','GNLY','MS4A2','TPSAB1','KIT','CLEC4C','LILRA4']
 more=[]
 for c in a.obs.cluster.cat.categories:
  mask=(a.obs.cluster==c).values
  for g in extra:
   ix=np.flatnonzero(a.var.symbol.values==g)
   counts_g=np.asarray(a.layers['counts'][:,ix].sum(axis=1)).ravel()
   log_g=np.asarray(a.X[:,ix].sum(axis=1)).ravel()
   more.append({'sample':r['sample'],'cluster':c,'gene':g,'pct_detected':float((counts_g[mask]>0).mean()*100),'mean_log1p':float(log_g[mask].mean())})
 e=pd.concat([e,pd.DataFrame(more)],ignore_index=True).drop_duplicates(['cluster','gene'],keep='last')
 e.to_csv(out/'marker_evidence_reviewed.tsv',sep='\t',index=False)
 dot+=extra
 clusters=list(a.obs.cluster.cat.categories)
 means=e.pivot(index='cluster',columns='gene',values='mean_log1p').reindex(index=clusters,columns=dot)
 pct=e.pivot(index='cluster',columns='gene',values='pct_detected').reindex(index=clusters,columns=dot)
 z=(means-means.mean())/means.std(ddof=0).replace(0,1)
 fig,ax=plt.subplots(figsize=(16,max(3,len(clusters)*.28+1.5)))
 xx,yy=np.meshgrid(np.arange(len(dot)),np.arange(len(clusters)))
 pts=ax.scatter(xx.ravel(),yy.ravel(),s=pct.values.ravel()*1.2,c=z.values.ravel(),cmap='RdBu_r',vmin=-2,vmax=2,edgecolors='none')
 labels={c:str(a.obs.loc[a.obs.cluster==c,'lineage_reviewed'].iloc[0]) for c in clusters}
 ax.set_xticks(range(len(dot)),dot,rotation=90);ax.set_yticks(range(len(clusters)),['C'+str(c)+' '+labels[c] for c in clusters]);ax.invert_yaxis();ax.set_title(r['short']+' markers');ax.set_xlabel('Gene');ax.set_ylabel('Cluster')
 fig.colorbar(pts,ax=ax,label='Mean z');fig.tight_layout();fig.savefig(out/'Markers.png',dpi=180,bbox_inches='tight');plt.close(fig)
 for c in clusters:
  idx=a.obs.cluster==c
  identity_rows.append({'sample':r['sample'],'cluster':c,'cells':int(idx.sum()),'lineage_reviewed':labels[c]})
 a.write_h5ad(out/'annotated.review_tmp.h5ad',compression='gzip')
 (out/'annotated.review_tmp.h5ad').replace(out/'annotated.h5ad')
 fraction=r['doublets']/r['qc_pass']
 ni=int(a.obs.cnv_ref_immune_reviewed.sum());ne=int(a.obs.cnv_ref_endo_reviewed.sum());nt=int(a.obs.cnv_target_reviewed.sum())
 cnvrows.append({'sample':r['sample'],'immune_ref_candidates':ni,'endo_ref_candidates':ne,'target_candidates':nt})
 flags.append({'sample':r['sample'],'doublet_fraction':fraction,'doublet_review':'low_call_rate_review' if fraction<.001 else ('high_call_rate_review' if fraction>.15 else 'routine_review'),'cnv_reference_review':'few_reference_candidates' if ni+ne<50 else 'review_reference_identity'})
 assert a.n_obs==r['retained'] and a.layers['counts'].shape==a.shape
 assert np.isfinite(a.layers['counts'].data).all() and (a.layers['counts'].data>=0).all()
 print('REVIEW_DONE',r['sample'],flush=True)
pd.DataFrame(audit).to_csv(args.results/'annotation_review_audit.tsv',sep='\t',index=False)
pd.DataFrame(flags).to_csv(args.results/'next_stage_review.tsv',sep='\t',index=False)
pd.DataFrame(cnvrows).to_csv(args.results/'cnv_candidates_reviewed.tsv',sep='\t',index=False)
ident=pd.DataFrame(identity_rows);ident.to_csv(args.results/'cluster_identity_reviewed.tsv',sep='\t',index=False)
counts=ident.pivot_table(index='sample',columns='lineage_reviewed',values='cells',aggfunc='sum',fill_value=0)
counts.to_csv(args.results/'lineage_counts_reviewed.tsv',sep='\t')
lines=['GSE270231 本轮实跑结果','',f"输入 {s.input_cells.sum():,} 个细胞/核；质量过滤后 {s.qc_pass.sum():,}；Scrublet 预测双细胞 {s.doublets.sum():,}；聚类保留 {s.retained.sum():,}。",'','逐样本结果']
for r in s.to_dict('records'):
 c=counts.loc[r['sample']];top=c[c>0].sort_values(ascending=False).head(3)
 lines.append(f"{r['short']}：输入 {r['input_cells']}，质量通过 {r['qc_pass']}，预测双细胞 {r['doublets']}，保留 {r['retained']}，Leiden 0.8 得到 {r['clusters']} 簇。复核后数量最多的三类："+'；'.join(f'{k}={int(v)}' for k,v in top.items())+'。')
lines+=['','分类判读','FibLike 为成纤维样标签，MesAmb 为间质身份重叠标签；均需结合 CNV 判定。当前尚未产生最终 CAF 数量、myCAF 亚群或流式分选 panel。','初筛对部分单核细胞的免疫身份检出不足，经前 50 位簇标志基因复核补充了髓系、T/NK、肥大和树突细胞等标签。自动标签和复核标签均保留，逐项修改证据见 annotation_review_audit.tsv。','下一阶段同时保留免疫与内皮参考候选，见 cnv_candidates_reviewed.tsv；参考身份及样本匹配还需结合表达与 CNV 结果复核。','低双细胞检出率样本保留自动门限及分数供复核，不把零检出解读为不存在双细胞。','各样本 UMAP 右图使用固定类别颜色，左图簇编号仅在本样本内有意义。']
(args.results/'本轮结果摘要.txt').write_text('\n'.join(lines)+'\n')
print('FINAL_REVIEW_COMPLETE',flush=True)
