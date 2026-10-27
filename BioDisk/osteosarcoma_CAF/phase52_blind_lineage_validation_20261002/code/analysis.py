from pathlib import Path
import os,sys,json,datetime,hashlib,gc
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');A=P/'phase41_full_transcriptome_setup_20261001';C=P/'phase43_full_transcriptome_discovery_20261001';O=P/'phase52_blind_lineage_validation_20261002';PRE=P/'phase51_matrix_lineage_attribution_20261002'
sys.path.insert(0,str(P/'phase46_spatial_functional_20261001/code'))
from common import test_table
import numpy as np,pandas as pd,anndata as ad
panels=json.loads((O/'resources/LINEAGE_PANELS_FROZEN.json').read_text());pr=json.loads((O/'resources/PROTOCOL_FROZEN.json').read_text())
names=list(panels);genes=pr['heldout_genes'];panelgenes=sorted(set(sum(panels.values(),[])))
assert not set(genes)&set(panelgenes)
parts={'full':panels,'A':{k:v[::2] for k,v in panels.items()},'B':{k:v[1::2] for k,v in panels.items()}}
assert all(not set(parts['A'][k])&set(parts['B'][k]) for k in panels)
m=pd.read_csv(A/'tables/clinical_sample_registry.tsv',sep='\t',keep_default_na=False);m=m[m.dataset.eq('SCPCP000017')]
a=ad.read_h5ad(A/'objects/all_candidate_full_genes.h5ad',backed='r');ob=a.obs.copy();a.file.close();frozen=set(ob.index)
cw=pd.read_csv(C/'tables/HGNC_full_gene_symbol_crosswalk.tsv.gz',sep='\t',keep_default_na=False)
symbols=dict(zip(cw.feature,cw.canonical_symbol));quality=dict(zip(cw.feature,~cw.mapping_method.eq('HGNC_Ensembl_approved')))
annotation=(O/'resources/receiver_annotation_frozen.py').read_text();flag=lambda x:x.astype(str).str.lower().isin(['true','1'])
def winner(v):
 safe=np.where(np.isfinite(v),v,-np.inf);mx=safe.max(1);ix=safe.argmax(1);lab=np.asarray(names,dtype=object)[ix]
 unique=(np.sum(np.isclose(safe,mx[:,None],atol=1e-12,rtol=0),axis=1)==1)&np.isfinite(mx)&(mx>0)
 lab[~unique]='Unresolved';order=np.sort(safe,axis=1);margin=np.where(np.isfinite(order[:,-2]),mx-order[:,-2],np.nan)
 return lab,margin
cells_out=[];input_qa=[];coverage=[];samplegroups=[];features=[];qc_raw=[]
for row in m.itertuples(index=False):
 f=P/'phase26_client_reanalysis_20260922/inputs'/row.sample/'full_called_counts.h5ad'
 a=ad.read_h5ad(f);x=a.X.tocsr();meta=a.obs.copy();lib=np.asarray(x.sum(1)).ravel().astype(float)
 cells=pd.read_csv(f.parent/'cells.tsv',sep='\t',keep_default_na=False).set_index('cell_id');assert np.array_equal(lib,pd.to_numeric(cells.loc[a.obs_names,'n_counts']).to_numpy())
 nonfrozen=~meta.index.isin(frozen);exec(annotation)
 pref={};heldout_idx=[]
 for j,en in enumerate(a.var_names.astype(str)):
  g=symbols.get(en,'');key=(quality.get(en,1),en)
  if g in genes:heldout_idx.append(j)
  if g in panelgenes+genes and (g not in pref or key<pref[g][0]):pref[g]=(key,j)
 assert set(genes).issubset(pref)
 rawtarget=x[:,[pref[g][1] for g in genes]].toarray().astype(np.float64)
 excluded_counts=np.asarray(x[:,heldout_idx].astype(np.float64).sum(1)).ravel()
 cleanlib=lib-excluded_counts;assert (cleanlib>=0).all()
 # Counterfactual change in all held-out counts cannot alter the scoring library.
 added=np.arange(a.n_obs,dtype=np.float64)%101
 assert np.array_equal((lib+added)-(excluded_counts+added),cleanlib)
 raw=np.full((a.n_obs,len(panelgenes)),np.nan)
 available=[g for g in panelgenes if g in pref]
 raw[:,[panelgenes.index(g) for g in available]]=x[:,[pref[g][1] for g in available]].toarray().astype(np.float64)
 log=np.log1p(raw/np.where(cleanlib>0,cleanlib,np.nan)[:,None]*1e4)
 scores={};labels={};margins={}
 for split,panel in parts.items():
  v=[]
  for lineage,gs in panel.items():
   have=[g for g in gs if g in pref];idx=[panelgenes.index(g) for g in have]
   v.append(log[:,idx].mean(1) if idx else np.full(a.n_obs,np.nan))
   coverage.append(dict(sample=row.sample,split=split,lineage=lineage,requested=len(gs),measured=len(have),missing=';'.join(sorted(set(gs)-set(have)))))
  scores[split]=np.column_stack(v);labels[split],margins[split]=winner(scores[split])
 agree=(labels['full']==labels['A'])&(labels['A']==labels['B'])&(labels['full']!='Unresolved')
 consensus=np.where(agree,labels['full'],np.where((labels['full']=='Unresolved')|(labels['A']=='Unresolved')|(labels['B']=='Unresolved'),'Unresolved','Panel_discordant'))
 unk=provider.str.strip().str.lower().isin(['unknown','unannotated','','nan','na']).to_numpy()
 hunk=hist.str.strip().str.lower().isin(['unknown','unannotated','unresolved','','nan','na']).to_numpy()
 query=nonfrozen&unk&hunk;primary=query&mask['Other_nonfrozen_unresolved']
 # Existing annotations, never predictions, define benchmark controls.
 ctrl={'Fibroblast':known['Fibroblast_other'],'Mural':known['Mural'],'T_NK':known['T_NK'],'Myeloid':known['Myeloid'],'Endothelial':known['Endothelial']}
 ctrl.update({
 'Osteogenic':provider.str.contains(r'osteoblast|osteogenic',regex=True).to_numpy()|hist.isin(['osteolike','osteogenic','osteoblast']).to_numpy(),
 'B_Plasma':provider.str.contains(r'\bb cell\b|plasma cell|b lineage',regex=True).to_numpy()|hist.isin(['b_plasma','b','plasma']).to_numpy(),
 'Epithelial':provider.str.contains('epithelial').to_numpy()|hist.isin(['epi','epithelial']).to_numpy(),
 'Mast':provider.str.contains('mast cell').to_numpy()|hist.eq('mast').to_numpy(),
 'Myogenic':(provider.str.contains(r'muscle|myoblast|myocyte',regex=True)&~provider.str.contains(r'smooth|pericyte|perivascular',regex=True)).to_numpy()|hist.isin(['myogenic','myocyte','skeletal_muscle']).to_numpy(),
 'Osteoclast':provider.str.contains('osteoclast').to_numpy()|hist.isin(['oc_like','osteoclast']).to_numpy(),
 'Chondrocyte':provider.str.contains('chondrocyte').to_numpy()|hist.isin(['chondrocyte','chondrogenic']).to_numpy(),
 'Erythroid':provider.str.contains(r'erythrocyte|erythroid',regex=True).to_numpy()|hist.isin(['erythroid']).to_numpy()})
 bits=np.column_stack([ctrl[n] for n in names]);cn=bits.sum(1);control=np.full(a.n_obs,'Unmapped_annotation',dtype=object)
 for j,n in enumerate(names):control[bits[:,j]&(cn==1)&nonfrozen]=n
 control[(cn>1)&nonfrozen]='Annotation_conflict';control[~nonfrozen]='Frozen_candidate'
 assert not np.any(query&np.isin(control,names))
 here=ob[ob['sample'].astype(str).eq(row.sample)];vv=here[here.in_supported.astype(bool)]
 isF3=meta.index.isin(vv.index[vv.supported_class_id.astype(str).eq('F3')])
 d=pd.DataFrame({'cell_id':a.obs_names.astype(str),'sample':row.sample,'patient_id':row.donor_record,'site':row.site,'library':lib,'library_without_targets':cleanlib,'query_all_unassigned':query,'query_primary_other':primary,'control_annotation':control,'CNV_evidence':evidence.CNV_evidence.reindex(meta.index).astype(str).to_numpy(),'reference':reference,'supported_F3':isF3,'full_label':labels['full'],'half_A_label':labels['A'],'half_B_label':labels['B'],'consensus_label':consensus,'full_margin':margins['full'],'half_A_margin':margins['A'],'half_B_margin':margins['B']})
 for split in parts:
  for j,n in enumerate(names):d[split+'__'+n]=scores[split][:,j]
 for j,g in enumerate(genes):d[g+'_raw']=rawtarget[:,j]
 cells_out.append(d)
 for scheme,lab in [('full_panel',labels['full']),('split_panel_consensus',consensus)]:
  groups={'query_Fibroblast':query&(lab=='Fibroblast'),'query_Osteogenic':query&(lab=='Osteogenic'),'supported_F3':isF3}
  for group,v in groups.items():
   for j,g in enumerate(genes):samplegroups.append(dict(scheme=scheme,sample=row.sample,patient_id=row.donor_record,site=row.site,group=group,gene=g,cells=int(v.sum()),counts=float(rawtarget[v,j].sum()),library=float(lib[v].sum())))
 # Deterministic independent score-check sample: every cell of PT cohort is not required to validate arithmetic.
 # Save all cells' raw marker counts; R independently reconstructs every score and winner.
 qr=pd.DataFrame(raw,columns=panelgenes);qr.insert(0,'cell_id',a.obs_names.astype(str));qr.insert(1,'sample',row.sample);qr.insert(2,'library_without_targets',cleanlib);qc_raw.append(qr)
 for g in panelgenes+genes:features.append(dict(sample=row.sample,gene=g,present=g in pref,feature_id=pref[g][0][1] if g in pref else '',role='heldout_target' if g in genes else 'classification_marker'))
 input_qa.append(dict(sample=row.sample,cells=a.n_obs,query_all=int(query.sum()),query_primary=int(primary.sum()),controls=int(np.isin(control,names).sum()),heldout_feature_columns=len(heldout_idx),target_genes_in_panels=0,nonpositive_scoring_library=int((cleanlib<=0).sum()),input_SHA256=hashlib.sha256(f.read_bytes()).hexdigest()))
 print('SAMPLE_COMPLETE',row.sample,a.n_obs,int(query.sum()),flush=True);del a,x,raw,log;gc.collect()
allcells=pd.concat(cells_out,ignore_index=True);assert len(allcells)==167180 and allcells.patient_id.nunique()==20
allcells.to_csv(O/'tables/cell_lineage_scores_and_predictions.tsv.gz',sep='\t',index=False)
pd.concat(qc_raw,ignore_index=True).to_csv(O/'tables/raw_marker_counts_for_QA.tsv.gz',sep='\t',index=False)
pd.DataFrame(input_qa).to_csv(O/'tables/input_QA.tsv',sep='\t',index=False)
pd.DataFrame(coverage).to_csv(O/'tables/panel_gene_coverage.tsv',sep='\t',index=False)
pd.DataFrame(features).to_csv(O/'tables/feature_mapping.tsv',sep='\t',index=False)
sg=pd.DataFrame(samplegroups);sg.to_csv(O/'tables/sample_target_gene_pseudobulk.tsv',sep='\t',index=False)
ps=sg.groupby(['scheme','patient_id','site','group','gene'])[['cells','counts','library']].sum().reset_index()
ps['log_expression']=np.where((ps.cells>0)&(ps.library>0),np.log2(ps.counts/ps.library*1e4+.1),np.nan)
ps.to_csv(O/'tables/patient_site_target_gene_expression.tsv',sep='\t',index=False)
contrasts={'query_Fibroblast_minus_query_Osteogenic':('query_Fibroblast','query_Osteogenic'),'query_Fibroblast_minus_supported_F3':('query_Fibroblast','supported_F3'),'query_Osteogenic_minus_supported_F3':('query_Osteogenic','supported_F3')}
rows=[]
for (scheme,patient,gene),t in ps.groupby(['scheme','patient_id','gene']):
 for contrast,(left,right) in contrasts.items():
  a=t[t.group.eq(left)][['site','log_expression']];b=t[t.group.eq(right)][['site','log_expression']]
  j=a.merge(b,on='site',suffixes=('_L','_R'));v=(j.log_expression_L-j.log_expression_R).to_numpy();v=v[np.isfinite(v)]
  rows.append(dict(scheme=scheme,patient_id=patient,gene=gene,contrast=contrast,effect=float(v.mean()) if len(v) else np.nan))
pat=pd.DataFrame(rows);pat.to_csv(O/'tables/patient_target_gene_effects.tsv',sep='\t',index=False)
st=test_table(pat,['scheme','contrast','gene'],'effect',lambda r:'phase52_heldout_gene_'+r.scheme);assert st.groupby('family').size().eq(9).all()
st.to_csv(O/'tables/target_gene_statistics.tsv',sep='\t',index=False)
# Benchmark against source annotations: descriptive per-patient agreement, not true external accuracy.
br=[]
controls=allcells[allcells.control_annotation.isin(names)]
for (patient,site,label),t in controls.groupby(['patient_id','site','control_annotation']):
 br.append(dict(patient_id=patient,site=site,annotation=label,cells=len(t),full_agreement=float(t.full_label.eq(label).mean()),consensus_coverage=float(t.consensus_label.isin(names).mean()),correct_consensus_fraction=float(t.consensus_label.eq(label).mean())))
bench=pd.DataFrame(br);bench.to_csv(O/'tables/control_site_agreement.tsv',sep='\t',index=False)
bp=bench.groupby(['patient_id','annotation'])[['full_agreement','consensus_coverage','correct_consensus_fraction']].mean().reset_index()
bp.to_csv(O/'tables/control_patient_agreement.tsv',sep='\t',index=False)
bs=bp.groupby('annotation').agg(N_patients=('patient_id','nunique'),mean_full_agreement=('full_agreement','mean'),mean_consensus_coverage=('consensus_coverage','mean'),mean_correct_consensus_fraction=('correct_consensus_fraction','mean')).reset_index()
bs.to_csv(O/'tables/control_agreement_summary.tsv',sep='\t',index=False)
conf=controls.groupby(['control_annotation','full_label']).size().rename('cells').reset_index();conf.to_csv(O/'tables/control_cell_confusion.tsv',sep='\t',index=False)
# Attribute target transcripts to prediction categories, retaining every query cell and disagreement class.
records=[]
for scope in ['query_primary_other','query_all_unassigned']:
 for scheme,col in [('full_panel','full_label'),('split_panel_consensus','consensus_label')]:
  q=allcells[allcells[scope]]
  for (sample,patient,site,label,cnv),t in q.groupby(['sample','patient_id','site',col,'CNV_evidence'],dropna=False):
   r=dict(scope=scope,scheme=scheme,sample=sample,patient_id=patient,site=site,predicted_label=label,CNV_evidence=cnv,cells=len(t))
   for g in genes:r[g]=float(t[g+'_raw'].sum())
   records.append(r)
rec=pd.DataFrame(records);rec.to_csv(O/'tables/query_sample_prediction_counts.tsv',sep='\t',index=False)
alltot=allcells.groupby(['patient_id','site'])[[g+'_raw' for g in genes]].sum().rename(columns={g+'_raw':g+'_all' for g in genes}).reset_index()
out=[]
for (scope,scheme),sub in rec.groupby(['scope','scheme']):
 site=sub.groupby(['patient_id','site','predicted_label','CNV_evidence'])[['cells']+genes].sum().reset_index()
 cats=site[['predicted_label','CNV_evidence']].drop_duplicates();loc=alltot[['patient_id','site']].drop_duplicates()
 grid=loc.merge(cats,how='cross').merge(site,on=['patient_id','site','predicted_label','CNV_evidence'],how='left')
 grid[['cells']+genes]=grid[['cells']+genes].fillna(0);grid=grid.merge(alltot,on=['patient_id','site'])
 qt=site.groupby(['patient_id','site'])[['cells']+genes].sum().rename(columns={g:g+'_query' for g in ['cells']+genes}).reset_index()
 grid=grid.merge(qt,on=['patient_id','site'],how='left')
 for g in genes:
  t=grid[['patient_id','site','predicted_label','CNV_evidence']].copy()
  t['gene']=g;t['all_gene_share']=grid[g]/grid[g+'_all'];t['within_query_gene_share']=grid[g]/grid[g+'_query'];t['query_cell_share']=grid.cells/grid.cells_query
  t['scope']=scope;t['scheme']=scheme;out.append(t)
share=pd.concat(out,ignore_index=True)
pp=share.groupby(['scope','scheme','patient_id','predicted_label','CNV_evidence','gene'])[['all_gene_share','within_query_gene_share','query_cell_share']].mean().reset_index()
pp.to_csv(O/'tables/query_patient_prediction_shares.tsv',sep='\t',index=False)
ss=pp.groupby(['scope','scheme','predicted_label','CNV_evidence','gene']).agg(N_patients=('patient_id','nunique'),mean_all_gene_share=('all_gene_share','mean'),mean_within_query_gene_share=('within_query_gene_share','mean'),mean_query_cell_share=('query_cell_share','mean')).reset_index()
ss.to_csv(O/'tables/query_prediction_share_summary.tsv',sep='\t',index=False)
for (scope,scheme),part in pp.groupby(['scope','scheme']):
 z=part.groupby(['patient_id','gene']).within_query_gene_share.sum(min_count=1).dropna();assert np.allclose(z,1,atol=1e-12)
old=pd.read_csv(PRE/'tables/main_annotation_CNV_patient_shares.tsv',sep='\t',float_precision='round_trip')
old=old[old.group.eq('Other_nonfrozen_unresolved')&old.annotation.eq('unassigned_annotation')].groupby(['patient_id','gene']).share.sum().reset_index()
errs=[]
for scheme,part in pp[pp.scope.eq('query_primary_other')].groupby('scheme'):
 v=part.groupby(['patient_id','gene']).all_gene_share.sum().reset_index().merge(old,on=['patient_id','gene'])
 err=float(np.max(np.abs(v.all_gene_share-v.share)));assert err<1e-12;errs.append(err)
qa={'at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'samples':len(m),'cells':len(allcells),'patients':20,'query_all_cells':int(allcells.query_all_unassigned.sum()),'query_primary_cells':int(allcells.query_primary_other.sum()),'annotation_control_cells':len(controls),'control_lineages':sorted(controls.control_annotation.unique()),'Osteogenic_annotation_controls':int(controls.control_annotation.eq('Osteogenic').sum()),'panel_genes':len(panelgenes),'target_marker_overlap':0,'classification_library_excludes_all_target_features':True,'counterfactual_heldout_count_change_invariant':True,'old_primary_query_RNA_share_max_error':max(errs),'target_gene_statistics':len(st),'FDR_families':2,'hypotheses_per_family':9,'annotations_not_overwritten':True,'identity_not_confirmed_by_prediction':True}
(O/'delivery/ANALYSIS_QA.json').write_text(json.dumps(qa,indent=2));(O/'ANALYSIS_COMPLETE').write_text('SUCCESS');print(json.dumps(qa),flush=True)
