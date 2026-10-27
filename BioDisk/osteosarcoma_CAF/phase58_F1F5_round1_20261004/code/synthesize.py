from common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import shutil
S=O/'delivery'; S.mkdir(exist_ok=True)
F=O/'figures'; F.mkdir(exist_ok=True)
def records(x): return json.loads(x.to_json(orient='records'))
pr=read(O/'tables/all_program_statistics.tsv.gz'); pc=read(O/'tables/corrected_program_statistics.tsv.gz')
pm=pr[pr.kind=='main'].merge(pc,on=['dataset','class_id','contrast','program'],suffixes=('_raw','_corrected'),how='outer')
pm['same_direction']=np.sign(pm.mean_delta_raw)==np.sign(pm.mean_delta_corrected)
pm.to_csv(O/'tables/program_correction_comparison.tsv.gz',sep='\t',index=False)
gs=[]
for f in sorted((O/'tables/genes').glob('*__corrected_statistics.tsv.gz')):
 ds=f.name.split('__')[1]; cr=read(f); raw=read(O/'tables/genes'/f'expanded__{ds}__all_statistics.tsv.gz');raw=raw[raw.kind=='main']
 z=raw.merge(cr,on=['dataset','class_id','contrast','feature'],suffixes=('_raw','_corrected'))
 z['same_direction']=np.sign(z.mean_delta_raw)==np.sign(z.mean_delta_corrected)
 z.to_csv(O/'tables/genes'/f'expanded__{ds}__raw_corrected_comparison.tsv.gz',sep='\t',index=False)
 for k,g in z.groupby('class_id'):
  nom=g[g.P_raw.lt(.05)];up=nom[nom.mean_delta_raw>0]
  gs.append(dict(dataset=ds,class_id=k,raw_nominal=len(nom),raw_nominal_same_corrected_direction=int(nom.same_direction.sum()),raw_nominal_corrected_P_lt_05=int(nom.P_corrected.lt(.05).sum()),raw_up=len(up),raw_up_corrected_same_direction=int(up.same_direction.sum()),raw_up_corrected_P_lt_05=int(up.P_corrected.lt(.05).sum())))
pd.DataFrame(gs).to_csv(O/'tables/gene_correction_summary.tsv',sep='\t',index=False)
gsea=[]
for f in sorted((O/'tables/gsea').glob('*.tsv.gz')):
 parts=f.name.replace('__GSEA.tsv.gz','').split('__');x=read(f);x['class_id']=parts[2].split('_')[0];x['dataset']=parts[1];x['library']=parts[3];x['source_file']=f.name
 gsea.append(x)
gg=pd.concat(gsea,ignore_index=True);gg.to_csv(O/'tables/all_GSEA_results.tsv.gz',sep='\t',index=False)
counts=gg.assign(nominal_model=gg.GSEA_model_P.lt(.05),model_BH_lt_05=gg.model_BH_FDR.lt(.05),positive_NES=gg.NES.gt(0)).groupby(['dataset','class_id','library']).agg(terms=('Term','size'),nominal_model=('nominal_model','sum'),model_BH_lt_05=('model_BH_lt_05','sum'),positive_NES=('positive_NES','sum')).reset_index()
counts.to_csv(O/'tables/GSEA_result_counts.tsv',sep='\t',index=False)
sc=read(O/'tables/surface_competition_statistics.tsv');sr=read(O/'tables/surface_candidate_registry.tsv');ss=[]
for (ds,k,gene),g in sc.groupby(['dataset','class_id','gene']):
 ss.append(dict(dataset=ds,class_id=k,gene=gene,backgrounds=len(g),evaluable_backgrounds=int(g.P.notna().sum()),positive_nominal_backgrounds=int(((g.mean_effect>0)&g.P.lt(.05)).sum()),negative_nominal_backgrounds=int(((g.mean_effect<0)&g.P.lt(.05)).sum()),higher_in_all_evaluable=bool((g.mean_effect[g.P.notna()]>0).all()) if g.P.notna().any() else False))
ss=pd.DataFrame(ss);ss.to_csv(O/'tables/surface_competition_summary.tsv',sep='\t',index=False)
reg=read(O/'tables/executed_contrast_registry.tsv');states=read(O/'tables/internal_state_composition.tsv')
internal=pr[(pr.kind=='internal')&pr.program.str.startswith('customer|')].copy();internal.to_csv(O/'tables/internal_customer_programs.tsv.gz',sep='\t',index=False)
nom=read(O/'tables/all_canonical_nominal_genes.tsv.gz');main=nom[nom.kind=='main']
qc=read(O/'tables/class_patient_QC_summary.tsv');lin=read(O/'tables/class_lineage_review_composition.tsv')
out={'main_N':records(reg[reg.kind=='main']), 'gene_summary':records(read(O/'tables/gene_statistics_summary.tsv').query("kind == 'main'")), 'QC':records(qc), 'states':records(states), 'lineage':records(lin), 'program_comparison':records(pm[pm.program.str.startswith('customer|')]), 'gene_correction':gs, 'surface':records(sr), 'surface_summary':records(ss), 'GSEA_counts':records(counts), 'GSEA_execution':read(O/'tables/GSEA_execution_summary.tsv').state.value_counts().to_dict(), 'internal_main_source':records(internal[internal.dataset=='SCPCP000017']), 'internal_registry':records(reg[reg.kind=='internal']), 'class_results':{}}
for k in ['F1','F2','F3','F4','F5']:
 a=main[(main.dataset=='SCPCP000017')&(main.class_id==k)&(main.mean_delta>0)].sort_values(['P','mean_delta'],ascending=[True,False])
 gh=gg[(gg.dataset=='SCPCP000017')&(gg.class_id==k)&(gg.library=='hallmark')]
 gr=gg[(gg.dataset=='SCPCP000017')&(gg.class_id==k)&(gg.library=='reactome')]
 out['class_results'][k]={'nominal_up_examples':records(a.head(12)),'hallmark':records(gh), 'reactome_positive_examples':records(gr[gr.NES>0].sort_values(['GSEA_model_P','NES'],ascending=[True,False]).head(10))}
(S/'report_data.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8')
plt.rcParams.update({'font.size':9,'pdf.fonttype':42,'savefig.dpi':180})
classes=['F1','F2','F3','F4','F5']
def heat(data,labels,title,path,vmin=None,vmax=None,cmap='RdBu_r',fmt='.2f',figsize=(8,5)):
 fig,ax=plt.subplots(figsize=figsize);im=ax.imshow(data,aspect='auto',cmap=cmap,vmin=vmin,vmax=vmax)
 ax.set_xticks(range(len(classes)),classes);ax.set_yticks(range(len(labels)),labels);ax.set_title(title,pad=14)
 for i in range(data.shape[0]):
  for j in range(data.shape[1]):
   v=data[i,j]
   if np.isfinite(v):ax.text(j,i,format(v,fmt),ha='center',va='center',fontsize=8,color='black',bbox=dict(facecolor='white',alpha=.55,edgecolor='none',pad=.3))
 fig.colorbar(im,ax=ax,shrink=.8);fig.tight_layout();fig.savefig(F/path,bbox_inches='tight');plt.close(fig)
r=reg[reg.kind=='main'].pivot(index='dataset',columns='target',values='N').reindex(columns=classes)
heat(r.to_numpy(float),list(r.index),'Independent matched patients per comparison','01_patient_N.png',0,15,'Blues','.0f',(8,4.5))
ps=['myCAF_customer','iCAF_customer','apCAF_customer','Osteogenic_customer','Mineral_customer','mCAF_matrix','stress_descriptive','cycle_descriptive']
for tag,df in [('raw',pr[pr.kind=='main']),('corrected',pc)]:
 t=df[(df.dataset=='SCPCP000017')&df.program.str.startswith('customer|')].pivot(index='program',columns='class_id',values='mean_delta').reindex(index=['customer|'+s for s in ps],columns=classes)
 heat(t.to_numpy(float),ps,'Patient-level program effects | '+tag,'02_program_'+tag+'.png',-2.5,2.5,figsize=(8,4.7))
hall=gg[(gg.dataset=='SCPCP000017')&(gg.library=='hallmark')].pivot(index='Term',columns='class_id',values='NES').reindex(columns=classes)
heat(hall.to_numpy(float),[s.replace('HALLMARK_','') for s in hall.index],'Full-rank GSEA | Hallmark NES | SCPCP000017','03_GSEA_all_hallmark.png',-2.5,2.5,figsize=(9,14))
fig,axs=plt.subplots(1,3,figsize=(10,3.3));q=qc[qc.dataset=='SCPCP000017'].set_index('class_id').reindex(classes)
for ax,col,title in zip(axs,['median_UMI','median_detected','median_MT_fraction'],['Median library UMI','Median detected genes','Median mitochondrial fraction']):
 ax.bar(classes,q[col],color=['#4477AA','#228833','#AA3377','#CCBB44','#66CCEE']);ax.set_title(title,fontsize=9)
fig.tight_layout();fig.savefig(F/'04_QC.png',bbox_inches='tight');plt.close(fig)
fig,ax=plt.subplots(figsize=(9,8));lab=states.round1_state.str.replace('expanded_','');ax.barh(lab,states.max_donor_fraction,color=[{'F1':'#4477AA','F2':'#228833','F3':'#AA3377','F4':'#CCBB44','F5':'#66CCEE'}[x] for x in states.class_id]);ax.invert_yaxis();ax.set_xlim(0,1);ax.set_xlabel('Fraction from the largest contributing donor');ax.set_title('Internal states: donor concentration (all sources)');fig.tight_layout();fig.savefig(F/'05_internal_donor_concentration.png',bbox_inches='tight');plt.close(fig)
for rel in ['phase30_two_layer_classification_20260926/figures/expanded_UMAP.png','phase30_two_layer_classification_20260926/figures/expanded_myCAF_UMAP.png','phase42_F3F4_identity_states_20261001/figures/expanded_F3_internal_UMAP.png','phase42_F3F4_identity_states_20261001/figures/expanded_F4_internal_UMAP.png']:
 src=P/rel
 if src.exists():shutil.copy2(src,F/src.name)
# Include provenance for the clinical files imported after the original freeze.
manifest=[]
for f in sorted((O/'tables/reused').iterdir()):
 if f.is_file():manifest.append(dict(path=str(f.relative_to(O)),bytes=f.stat().st_size,SHA256=sha(f)))
pd.DataFrame(manifest).to_csv(S/'reused_clinical_files_SHA256.tsv',sep='\t',index=False)
(O/'SYNTHESIS_COMPLETE').write_text('SUCCESS\n');print('SYNTHESIS_COMPLETE',len(gg),'GSEA term rows',flush=True)
