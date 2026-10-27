from pathlib import Path
import numpy as np,pandas as pd,json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
B=Path('C:/codex_test');O=B/'CAF_full_report_20260921';F=O/'figures'
s=pd.read_csv(B/'CAF_expanded_20260921/tables/sample_layer_counts.tsv',sep='\t')
e=pd.read_csv(B/'CAF_expanded_20260921/tables/all_cell_layer_evidence.tsv.gz',sep='\t')
summary=s.groupby('dataset').agg(samples=('sample','size'),QC_cells=('QC_cells','sum'),strict=('strict','sum'),supported=('supported','sum'),expanded=('expanded','sum')).reset_index()
summary.to_csv(O/'inputs/cohort_summary.tsv',sep='\t',index=False)
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
def finish(fig,name):fig.savefig(F/name,dpi=190,bbox_inches='tight');plt.close(fig)
fig,ax=plt.subplots(figsize=(10,5));ax.set_xlim(0,10);ax.set_ylim(0,5);ax.axis('off')
texts=[('1  Cohort selection','52 sample records / 8 sources'),('2  QC and lineage review','Counts, doublets, ambient RNA'),('3  CNV evidence','inferCNV + second method'),('4  Three candidate layers','Strict / supported / broad'),('5  Joint expression states','Shared measured genes / Harmony'),('6  Patient-level comparisons','C7 / C9; effect, P and FDR'),('7  Analytical flow cytometry','Surface expression / co-expression')]
coords=[(.1,3.8),(3.5,3.8),(6.9,3.8),(6.9,2.15),(3.5,2.15),(.1,2.15),(3.5,.4)]
for (title,sub),(x,y) in zip(texts,coords):
 ax.add_patch(FancyBboxPatch((x,y),3,.85,boxstyle='round,pad=.08',facecolor='#edf3f7',edgecolor='#638398'))
 ax.text(x+1.5,y+.57,title,ha='center',fontsize=11,weight='bold');ax.text(x+1.5,y+.22,sub,ha='center',fontsize=9)
for start,end in [((3.2,4.2),(3.4,4.2)),((6.6,4.2),(6.8,4.2)),((8.4,3.7),(8.4,3.1)),((6.8,2.55),(6.6,2.55)),((3.4,2.55),(3.2,2.55)),((1.6,2.05),(4.9,1.35))]:ax.annotate('',xy=end,xytext=start,arrowprops={'arrowstyle':'->','color':'#4b6678','lw':1.6})
finish(fig,'01_analysis_workflow.png')
q=pd.read_csv(B/'CAF_resume_20260921/tables/qc_summary.tsv',sep='\t')
loss=pd.read_csv(O/'inputs/ScPCA_QC_loss_decomposition.tsv',sep='\t').set_index('sample');q['official_singlets']=q['sample'].map(loss.official_singlets)
fig,ax=plt.subplots(figsize=(10,4.2));x=np.arange(len(q));w=.2
for i,(c,label,col) in enumerate([('filtered_cells','Filtered input','#d0d9e2'),('official_processed','Official processed','#9cbacf'),('official_singlets','Official singlets','#5c9baf'),('retained','Project retained','#205f82')]):ax.bar(x+(i-1.5)*w,q[c],w,label=label,color=col)
ax.set_xticks(x,q['sample'].str.replace('SCPCS000','S'),rotation=45,ha='right');ax.set_ylabel('Cells');ax.legend(ncol=3,fontsize=9);fig.tight_layout();finish(fig,'02_scpca_qc.png')
states=['both_high','high_without_low','low_high_conflict','low_no_high','no_low_evidence'];cols=['#943f3d','#c98267','#e5b75c','#4487a2','#c4c9cc'];labels=['Both high','High, no low','Low/high conflict','Low, no high','No low evidence']
ct=pd.crosstab(e.dataset,e.CNV_state).reindex(columns=states,fill_value=0);ct.to_csv(O/'inputs/CNV_state_counts.tsv',sep='\t')
fig,ax=plt.subplots(figsize=(10,4.5));left=np.zeros(len(ct))
for st,col,lab in zip(states,cols,labels):
 v=ct[st].to_numpy()/ct.sum(axis=1).to_numpy()*100;ax.barh(np.arange(len(ct)),v,left=left,color=col,label=lab);left+=v
ax.set_yticks(range(len(ct)),ct.index);ax.invert_yaxis();ax.set_xlabel('Percent of all QC cells');ax.legend(loc='upper center',bbox_to_anchor=(.4,1.21),ncol=3,fontsize=9);fig.tight_layout();finish(fig,'03_CNV_states.png')
a=pd.read_csv(O/'inputs/infercnv_cells.tsv',sep='\t').set_index('cell_id');b=pd.read_csv(O/'inputs/copykat_cells.tsv',sep='\t').set_index('cell_id').loc[a.index];m=pd.read_csv(O/'inputs/chromosome_log2_medians.tsv',sep='\t',index_col=0)
strict=set(e.loc[e.strict,'cell_id']);rng=np.random.default_rng(20260921)
groups={'Training normal reference':a.cnv_role.str.startswith('ref_'),'Strict dual-low candidates':a.index.isin(strict),'Both methods high':a.call.eq('elevated_signal')&b.call.eq('aneuploid')}
selected=[];counts=[]
for name,mask in groups.items():
 pool=a.loc[mask,'tool_id'];pool=pool[pool.isin(m.columns)];chosen=rng.choice(pool.to_numpy(),min(30,len(pool)),replace=False);selected.extend(chosen);counts.append((name,len(chosen),len(pool)))
fig,ax=plt.subplots(figsize=(10,4.7));im=ax.imshow(m.loc[:,selected],aspect='auto',cmap='RdBu_r',vmin=-.3,vmax=.3,interpolation='nearest');cum=0;xt=[];xl=[]
for name,n,total in counts:
 xt.append(cum+(n-1)/2);xl.append(name+'\n'+str(n)+' shown / '+str(total)+' available');cum+=n
 if cum<len(selected):ax.axvline(cum-.5,color='black',lw=.7)
ax.set_xticks(xt,xl,fontsize=9);ax.set_yticks(range(0,len(m),2),m.index[::2]);ax.set_ylabel('Chromosome');ax.set_title('GSE162454 / GSM4952365: observed inferCNV chromosome medians',fontsize=11);fig.colorbar(im,ax=ax,label='Median log2 denoised expression',fraction=.025);fig.tight_layout();finish(fig,'04_actual_CNV_example.png')
cal=pd.read_csv(B/'CAF_resume_20260921/tables/CNV_heldout_acceptance.tsv',sep='\t');piv=cal.pivot(index='sample',columns='method',values='heldout_low_fraction')
fig,ax=plt.subplots(figsize=(9,2.8));xx=np.arange(len(piv));ax.bar(xx-.18,piv.infercnv,width=.36,label='inferCNV');ax.bar(xx+.18,piv.copykat,width=.36,label='CopyKAT');ax.axhline(.9,color='gray',ls='--');ax.set_ylim(.85,1.025);ax.set_xticks(xx,piv.index.str.replace('SCPCS000','S'));ax.set_ylabel('Low-signal fraction');ax.legend(ncol=2,fontsize=9);fig.tight_layout();finish(fig,'05_reference_calibration.png')
fig,axs=plt.subplots(1,2,figsize=(11,4.4));xx=np.arange(len(summary));w=.24
for i,(c,col) in enumerate([('strict','#205f82'),('supported','#5da293'),('expanded','#dfab4e')]):axs[0].bar(xx+(i-1)*w,summary[c],w,label=c,color=col)
axs[0].set_yscale('symlog',linthresh=100);axs[0].set_xticks(xx,summary.dataset,rotation=55,ha='right',fontsize=9);axs[0].set_ylabel('Candidate cells (symlog)');axs[0].legend(fontsize=9)
metrics=['CNV conflict','No low CNV evidence','Mural program','Ambient not checked'];positions=np.arange(4)
for i,(layer,col) in enumerate([('strict','#205f82'),('supported','#5da293'),('expanded','#dfab4e')]):
 g=e[e[layer]];v=[g.CNV_state.eq('low_high_conflict').mean(),g.CNV_state.eq('no_low_evidence').mean(),g.mural_flag.mean(),(~g.ambient_checked).mean()];axs[1].bar(positions+(i-1)*w,np.array(v)*100,w,color=col,label=layer)
axs[1].set_xticks(positions,metrics,rotation=25,ha='right',fontsize=9);axs[1].set_ylabel('Flagged cells (%)');fig.tight_layout();finish(fig,'06_layers_sources_and_risks.png')
xy=pd.read_csv(O/'inputs/supported_umap.tsv.gz',sep='\t',index_col=0);ob=pd.read_csv(B/'CAF_expanded_20260921/tables/supported_cell_labels.tsv.gz',sep='\t',index_col=0).loc[xy.index]
fig,axs=plt.subplots(1,2,figsize=(11,4.6));v=xy.to_numpy();axs[0].scatter(v[:,0],v[:,1],s=1,color='#d3d6d8',rasterized=True)
for cl,col in [(7,'#c65b4c'),(9,'#2779a0')]:
 mask=ob.cluster==cl;axs[0].scatter(v[mask,0],v[mask,1],s=2,color=col,label='C'+str(cl),rasterized=True);cent=np.median(v[mask],axis=0);axs[0].text(*cent,'C'+str(cl),fontsize=14,weight='bold',bbox={'facecolor':'white','alpha':.8,'edgecolor':'none'})
for name,g in ob.groupby('dataset',observed=True):
 ii=ob.index.get_indexer(g.index);axs[1].scatter(v[ii,0],v[ii,1],s=1,label=name,rasterized=True)
for ax in axs:ax.set_xticks([]);ax.set_yticks([])
axs[0].set_title('Supported layer: C7 and C9');axs[0].legend(fontsize=10);axs[1].set_title('Data source');axs[1].legend(fontsize=8,markerscale=6,loc='upper left',bbox_to_anchor=(1,1));fig.tight_layout();finish(fig,'07_supported_umap.png')
prof=pd.read_csv(B/'CAF_evidence_20260921/tables/donor_cluster_programs.tsv',sep='\t');prof=prof[(prof.layer=='supported')&prof.known_patient&(prof.cells>=10)];pp=prof.groupby(['cluster','program']).mean_score.mean().unstack().reindex(columns=['myCAF_customer','iCAF_customer','apCAF_customer','Osteogenic_customer','Mineral_customer']);z=(pp-pp.mean())/pp.std()
fig,ax=plt.subplots(figsize=(6.5,5.5));im=ax.imshow(z,cmap='RdBu_r',vmin=-2,vmax=2,aspect='auto');ax.set_yticks(range(len(z)),['C'+str(i) for i in z.index]);ax.set_xticks(range(5),['myCAF','iCAF','apCAF','Osteogenic','Mineral'],rotation=25,ha='right');fig.colorbar(im,ax=ax,label='Column-standardized donor mean',fraction=.04);fig.tight_layout();finish(fig,'08_program_heatmap.png')
direct=pd.read_csv(B/'CAF_discovery_stats_20260921/tables/C7_vs_C9_patient_differences.tsv',sep='\t');g=direct[direct.population=='all_candidate_cells'];prs=['myCAF_customer','iCAF_customer','Osteogenic_customer','Mineral_customer'];fig,ax=plt.subplots(figsize=(9,3.8))
for j,pr in enumerate(prs):
 gg=g[g.program==pr];ax.scatter(gg.delta_C7_minus_C9,np.full(len(gg),j)+rng.uniform(-.08,.08,len(gg)),s=35,color='#357a99');ax.plot(np.median(gg.delta_C7_minus_C9),j,'|',color='black',markersize=17)
ax.axvline(0,color='gray',ls='--');ax.set_yticks(range(4),['myCAF','iCAF','Osteogenic','Mineral']);ax.set_xlabel('C7 minus C9, within the same patient');fig.tight_layout();finish(fig,'09_C7_C9_paired.png')
surface=pd.read_csv(B/'CAF_evidence_20260921/tables/cluster_surface_ranking.tsv',sep='\t');r=surface[(surface.layer=='supported')&surface.cluster.isin([7,9])&(surface.markers=='FAP+THY1')].sort_values('cluster');fig,ax=plt.subplots(figsize=(8,3.3));xx=np.arange(2)
for i,(key,label,col) in enumerate([('median_coverage','Target coverage','#347e9f'),('median_known_target_fraction','Known target in RNA gate','#e1a34b'),('median_CNV_high_fraction','CNV-high in RNA gate','#bb6256')]):ax.bar(xx+(i-1)*.22,r[key]*100,.22,label=label,color=col)
ax.set_xticks(xx,['C7','C9']);ax.set_ylabel('Median percentage');ax.legend(fontsize=9,ncol=1,loc='upper right');ax.set_ylim(0,65);fig.tight_layout();finish(fig,'10_surface_screen.png')
assert summary.QC_cells.sum()==299483 and len(e)==299483
print('FIGURES_COMPLETE',len(list(F.glob('*.png'))),counts)
