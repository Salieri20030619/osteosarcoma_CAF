from pathlib import Path
import os,json,hashlib,hmac,random,shutil,gc
import numpy as np,pandas as pd,h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase53_de_novo_reannotation_20261002');B=R/'blind_review';U=R/'unblinding'
maps=json.loads((U/'identifier_maps.json').read_text());cfg=json.loads((U/'anonymization_secret.json').read_text());salt=bytes.fromhex(cfg['salt'])
def cid(x):return 'C_'+hmac.new(salt,str(x).encode(),hashlib.sha256).hexdigest()[:24]
def sread(a):return np.array([x.decode() if isinstance(x,bytes) else str(x) for x in np.asarray(a).ravel()])
def bin_sort(x):
 chrom,pos=x.split(':');c=chrom.replace('chr','');return (int(c) if c.isdigit() else 23 if c=='X' else 24 if c=='Y' else 99,int(pos))
clusters=[];keys=[];templates=[];figure_manifest=[];allobs=[];h5str=h5py.string_dtype('utf-8')
for ds,source_id in sorted(maps['sources'].items(),key=lambda x:x[1]):
 d=R/'tables'/ds;dest=B/source_id;dest.mkdir(exist_ok=True)
 obs=pd.read_csv(d/'cell_structure.tsv.gz',sep='\t',dtype={'r05':str,'r10':str,'r20':str})
 cmaps={}
 for res in ['r05','r10','r20']:
  raw=sorted(obs[res].unique());coded=[f'{source_id}_{res.upper()}_K{j+1:03}' for j in range(len(raw))]
  rr=random.Random(int.from_bytes(hmac.new(salt,(ds+res).encode(),hashlib.sha256).digest()[:8],'big'));rr.shuffle(coded)
  cmaps[res]=dict(zip(raw,coded))
  keys.extend([dict(dataset=ds,source_id=source_id,resolution=res,original_cluster=a,cluster_id=b) for a,b in cmaps[res].items()])
 anon=pd.DataFrame({'cell_id':obs.cell_id.map(cid),'source_id':source_id,'sample_id':obs['sample'].map(maps['samples']),'patient_id':obs.patient_id.map(maps['patients']),'tissue':obs.site,'technology':obs.technology,'UMAP1':obs.UMAP1,'UMAP2':obs.UMAP2})
 for res in cmaps:anon[res]=obs[res].map(cmaps[res])
 anon.to_csv(dest/'cell_structure.tsv.gz',sep='\t',index=False)
 allobs.append(anon[['cell_id','source_id','sample_id','patient_id','r10']])
 meta=pd.read_csv(d/'primary_clusters.tsv',sep='\t',dtype={'cluster':str},keep_default_na=False)
 meta['cluster_id']=meta.cluster.map(cmaps['r10'])
 meta['full_library']=meta.cluster.map(obs.groupby('r10').full_library.sum())
 meta[['cluster_id','cells','patients','samples','max_patient_fraction','analysis_library','full_library']].to_csv(dest/'clusters.tsv',sep='\t',index=False)
 genes=pd.read_csv(d/'genes.tsv',sep='\t',keep_default_na=False)
 genes.to_csv(dest/'genes.tsv',sep='\t',index=False)
 shutil.copy2(R/'objects'/ds/'primary_cluster_counts.mtx.gz',dest/'primary_cluster_counts.mtx.gz')
 tops=[];fields=['mean_log_expression','rest_mean_log_expression','detection_fraction','patient_mean_log2_difference','positive_patients','negative_patients']
 with h5py.File(dest/'all_cluster_marker_evidence.h5','w') as h:
  h.create_dataset('feature_id',data=genes.feature.astype(str).to_numpy(object),dtype=h5str)
  h.create_dataset('gene_symbol',data=genes.symbol.astype(str).to_numpy(object),dtype=h5str)
  h.create_dataset('cluster_id',data=meta.cluster_id.to_numpy(object),dtype=h5str)
  for f in fields:h.create_dataset(f,shape=(len(genes),len(meta)),dtype='f8',chunks=(min(2048,len(genes)),1),compression='gzip')
  ns=[]
  for j,row in enumerate(meta.itertuples(index=False)):
   full=pd.read_csv(d/f'cluster_{row.cluster}_full_marker_evidence.tsv.gz',sep='\t',float_precision='round_trip',keep_default_na=False,na_values=[''])
   assert full.feature.tolist()==genes.feature.tolist()
   for field in fields:h[field][:,j]=full[field].to_numpy(float)
   ns.append(int(full.N_patients_matched.iloc[0]))
   t=pd.read_csv(d/f'cluster_{row.cluster}_top_markers.tsv',sep='\t')
   t=t.drop(columns=['dataset','cluster']);t.insert(0,'cluster_id',row.cluster_id);tops.append(t)
  h.create_dataset('N_patients_matched',data=np.array(ns,int))
 top=pd.concat(tops,ignore_index=True);top.to_csv(dest/'top_descriptive_markers.tsv',sep='\t',index=False)
 for row in meta.itertuples(index=False):
  templates.append({'source_id':source_id,'cluster_id':row.cluster_id,'cells':row.cells,'patients':row.patients,'reviewed_lineage':'','malignancy_assessment':'','supporting_genes':'','genomic_evidence':'','uncertainty':'','reviewer_notes':''})
 composition=anon.groupby(['r10','patient_id']).size().rename('cells').reset_index().rename(columns={'r10':'cluster_id'})
 composition.to_csv(dest/'cluster_patient_composition.tsv',sep='\t',index=False)
 # Expression-space maps show only anonymous identifiers.
 fig,axs=plt.subplots(1,2,figsize=(12,5.3))
 for ax,col,title in zip(axs,['r10','patient_id'],['Primary clusters','Patients']):
  levels=sorted(anon[col].unique());codes=pd.Categorical(anon[col],categories=levels).codes
  ax.scatter(anon.UMAP1,anon.UMAP2,c=codes,s=.5,cmap='turbo',rasterized=True,alpha=.65)
  ax.set_title(title);ax.set_xlabel('UMAP 1');ax.set_ylabel('UMAP 2')
  if col=='r10':
   for label,t in anon.groupby(col):
    ax.text(t.UMAP1.median(),t.UMAP2.median(),label.split('_')[-1],fontsize=6,ha='center',bbox=dict(facecolor='white',alpha=.6,edgecolor='none',pad=.2))
 fig.suptitle(source_id+' | expression-space view');fig.tight_layout();fig.savefig(dest/'expression_space.png',dpi=150);fig.savefig(dest/'expression_space.pdf');plt.close(fig)
 figure_manifest.append(str((dest/'expression_space.png').relative_to(B)))
 # Genome evidence remains quantitative; no predicted lineage labels are exposed.
 ckrows=[]
 obsidx=obs.set_index('cell_id')
 for sample,sub in obs.groupby('sample'):
  f=R/'objects/cnv'/f'{sample}_copykat_profiles.h5'
  if not f.exists():continue
  with h5py.File(f,'r') as h:
   cells=sread(h['cell_id'][:]);bins=sread(h['bin_id'][:]);v=h['profiles'][:]
  if v.shape!=(len(cells),len(bins)):v=v.T
  assert v.shape==(len(cells),len(bins)) and set(cells).issubset(obsidx.index)
  ct=obsidx.loc[cells]
  assert ct['sample'].eq(sample).all()
  for cl,indices in ct.groupby('r10').indices.items():
   vals=v[indices].mean(0)
   entry={'sample_id':maps['samples'][sample],'patient_id':maps['patients'][ct.patient_id.iloc[0]],'tissue':ct.site.iloc[0],'cluster_id':cmaps['r10'][cl],'profile_cells':len(indices)}
   entry.update(dict(zip(bins,vals)));ckrows.append(entry)
  del v;gc.collect()
 if ckrows:
  ck=pd.DataFrame(ckrows);bincols=sorted([x for x in ck if x.startswith('chr')],key=bin_sort)
  ck[['sample_id','patient_id','tissue','cluster_id','profile_cells']+bincols].to_csv(dest/'copykat_sample_cluster_profiles.tsv.gz',sep='\t',index=False)
  site=ck.groupby(['patient_id','tissue','cluster_id'])[bincols].mean().reset_index()
  patient=site.groupby(['patient_id','cluster_id'])[bincols].mean().reset_index()
  patient.to_csv(dest/'copykat_patient_cluster_profiles.tsv.gz',sep='\t',index=False)
  mean=patient.groupby('cluster_id')[bincols].mean().reindex(meta.cluster_id)
  finite=np.abs(mean.to_numpy());finite=finite[np.isfinite(finite)]
  scale=max(float(np.quantile(finite,.95)),.05) if len(finite) else .05
  fig,ax=plt.subplots(figsize=(14,max(5,len(meta)*.14)))
  cm=plt.get_cmap('RdBu_r').copy();cm.set_bad('#dddddd')
  im=ax.imshow(mean,aspect='auto',cmap=cm,vmin=-scale,vmax=scale,interpolation='none')
  starts=[];labels=[]
  for j,binid in enumerate(bincols):
   chrom=binid.split(':')[0]
   if not labels or labels[-1]!=chrom:starts.append(j);labels.append(chrom)
  for x in starts:ax.axvline(x-.5,color='black',lw=.25)
  ax.set_xticks(starts,[x.replace('chr','') for x in labels],fontsize=8);ax.set_yticks(range(len(meta)),[x.split('_')[-1] for x in meta.cluster_id],fontsize=6)
  ax.set_xlabel('Chromosome | native coordinates, 5 Mb display bins');ax.set_ylabel('Anonymous primary cluster')
  ax.set_title(source_id+' | patient-mean native CopyKAT signal')
  fig.colorbar(im,ax=ax,label='Native signal (display scale only)',shrink=.8);fig.tight_layout()
  fig.savefig(dest/'copykat_cluster_profiles.png',dpi=150);fig.savefig(dest/'copykat_cluster_profiles.pdf');plt.close(fig)
  figure_manifest.append(str((dest/'copykat_cluster_profiles.png').relative_to(B)))
 else:(dest/'COPYKAT_PROFILE_UNAVAILABLE.txt').write_text('No valid native continuous profiles available for this source. Missing evidence is not a normal call.')
 q=json.loads((d/'STRUCTURE_QA.json').read_text())
 safe={k:q[k] for k in ['cells','samples','patients','source_features_union','source_features_common','structure_HVG','clusters','resolution_ARI_05_10','resolution_ARI_10_20','nonpositive_analysis_library']}
 (dest/'technical_summary.json').write_text(json.dumps(safe,indent=2))
 print('BLIND_VIEW_COMPLETE',source_id,len(meta),len(obs),flush=True)
pd.DataFrame(keys).to_csv(U/'cluster_identifier_key.tsv',sep='\t',index=False)
pd.DataFrame(templates).sort_values(['source_id','cluster_id']).to_csv(B/'reviewer_response_template.tsv',sep='\t',index=False)
pd.concat(allobs,ignore_index=True).to_csv(B/'anonymous_cell_cluster_roster.tsv.gz',sep='\t',index=False)
# Native HMM genomic fractions/segments are anonymized; original native group names removed.
out=B/'regional_CNV';out.mkdir(exist_ok=True)
for f in sorted((R/'objects/infercnv_regions').glob('*')):
 matches=[sm for sm in maps['samples'] if f.name.startswith(sm+'__')];assert len(matches)==1
 sample=matches[0];tail=f.name[len(sample)+2:];samplecode=maps['samples'][sample]
 context=tail.replace('_groups.tsv.gz','').replace('_regions.tsv.gz','')
 t=pd.read_csv(f,sep='\t',float_precision='round_trip')
 if 'cell_id' in t:t['cell_id']=t.cell_id.map(cid)
 if 'cell_group_name' in t:
  t['group_id']=t.cell_group_name.map(lambda v:'G_'+hmac.new(salt,(sample+'|'+context+'|'+str(v)).encode(),hashlib.sha256).hexdigest()[:12])
  t=t.drop(columns=['cell_group_name'])
 if 'cnv_name' in t:t=t.drop(columns=['cnv_name'])
 t.to_csv(out/(samplecode+'__'+tail),sep='\t',index=False)
# Add continuous per-cell metrics, but not native discrete calls or historic cell labels.
cm=pd.read_csv(R/'tables/copykat_cell_continuous_metrics.tsv.gz',sep='\t')
cm['cell_id']=cm.cell_id.map(cid);cm['sample_id']=cm['sample'].map(maps['samples']);cm=cm.drop(columns='sample')
cm.to_csv(B/'copykat_cell_continuous_metrics.tsv.gz',sep='\t',index=False)
pd.DataFrame({'figure':figure_manifest}).to_csv(B/'FIGURE_INDEX.tsv',sep='\t',index=False)
qa={'primary_clusters':len(templates),'cells':sum(1 for _ in []),'sources':len(maps['sources']),'figures':len(figure_manifest),'proposed_labels_not_read_by_exporter':True,'original_identifiers_mapped':True}
qa['cells']=sum(len(x) for x in allobs);assert qa['cells']==876648 and qa['primary_clusters']==307
(R/'delivery/BLIND_VIEW_QA.json').write_text(json.dumps(qa,indent=2));(R/'BLIND_VIEWS_COMPLETE').write_text('SUCCESS');print(json.dumps(qa),flush=True)
