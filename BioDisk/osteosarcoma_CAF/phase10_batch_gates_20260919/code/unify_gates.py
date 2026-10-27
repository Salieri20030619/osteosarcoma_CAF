from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
import json,gzip,hashlib,datetime,warnings
import numpy as np,pandas as pd,anndata as ad
from scipy.io import mmread
from scipy import sparse as sp
from threadpoolctl import threadpool_limits
threadpool_limits(4)
warnings.simplefilter('ignore',pd.errors.PerformanceWarning)
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase10_batch_gates_20260919';assert (O/'EARLY_AMBIENT_COMPLETE').exists()
panels=json.loads((P/'geo_preprocess_20260918/review/lineage_rules.json').read_text())['panels'];panels['OC_specific']=['OCSTAMP','DCSTAMP','ATP6V0D2','SIGLEC15','SUCNR1','CALCR']
programs=json.loads((P/'phase9_joint_20260919/resources/programs.json').read_text());panels_export=dict(panels)
protocol={'frozen_before_evaluation':datetime.datetime.now(datetime.timezone.utc).isoformat(),'version':'unified_gates_v1','scope':'All 40 previously QC-retained sample records; common molecular and CNV gate code; frozen broad-lineage scope retained with explicit provenance. No rescue of unassessed CNV.','core_scope':'FibLike reviewed scope; competing Mural/SMC and MesAmb retained separately','primary':{'n_genes_min':500,'Fib_detected_genes_min':2,'raw_and_corrected_detection_counts':1,'competing_program_percent':.2,'osteogenic_program_percent':.2,'osteogenic_detected_genes_min':3,'cluster_any_CNV_high_fraction_max_exclusive':.2,'cluster_both_CNV_high_fraction_max_exclusive':.2,'both_CNV_low_and_reference_checks':True,'ambient_checked_required':True},'sensitivities':{'osteo_sensitive':'only osteogenic percent .5; other main gates fixed','competitor_sensitive':'competing/mixed/osteo percent .5, OC-specific .25; not an accepted main-set expansion','mural_reference':'FibLike/MesAmb/Mural scope, Fib>=2, dual CNV low and all other main gates; Mural gate flagged but not used to exclude this reference set'},'excluded_records':'retained with all reasons; no physical deletion','interpretation':'operational candidates; thresholds are not validated identity truth','target_name':'SR1004; do not rename S100A4'}
(O/'PROTOCOL.json').write_text(json.dumps(protocol,indent=2));(O/'resources/panels.json').write_text(json.dumps(panels,indent=2));(O/'resources/programs.json').write_text(json.dumps(programs,indent=2))
old=ad.read_h5ad(P/'phase9_joint_20260919/objects/joint_CAF_integrated.h5ad',backed='r');oldids=set(old.obs_names);old.file.close()
allrows=[];man=[];outobjects=[];checks=[];cover=[];calibration=[];summary=[];gatecounts=[]
def boolcol(s):
 if s.dtype==bool:return s.fillna(False).to_numpy()
 return s.fillna('').astype(str).str.lower().isin(['true','1','1.0']).to_numpy()
def evidence(x,sy):
 x=x.tocsc();den=np.asarray(x.sum(1)).ravel();d={}
 for name,gs in panels.items():
  vals=[np.asarray(x[:,np.flatnonzero(sy==g)].sum(1)).ravel() for g in gs if (sy==g).any()];v=np.column_stack(vals) if vals else np.zeros((x.shape[0],0))
  d['det_'+name]=(v>=1).sum(1);d['pos_'+name]=(v>0).sum(1);d['pct_'+name]=v.sum(1)/np.maximum(den,1)*100
 return pd.DataFrame(d),den

def conflict(d,frac=.2,osteo=.2,mural=True,det='det_'):
 out=np.zeros(len(d),bool)
 for name,n,pct in [('Mural',3,frac),('OC',2,frac/2),('T_NK',3,frac),('Myeloid',4,frac),('Endo',3,frac),('Epi',3,frac),('Myogenic',2,frac)]:
  if name=='Mural' and not mural:continue
  pp='OC_specific' if name=='OC' else name;out|=(d[det+name]>=n)&(d['pct_'+pp]>pct)
 out|=(d[det+'Osteo']>=3)&(d.pct_Osteo>osteo)
 immune=sum(d['pct_'+k] for k in ['T_NK','B_Plasma','Myeloid','DC']);mes=d.pct_Fib+d.pct_Osteo+d.pct_Mural
 mixed=((immune>frac)&(mes>frac))|((d.pct_Endo>frac)&(mes>frac))|((d.pct_Epi>frac)&(mes>frac))
 return np.asarray(out),np.asarray(mixed)

sources=[('early',P/'phase3_identity_20260914/identity/all_cell_identity_evidence.tsv.gz'),('discovery',P/'phase7_review_states_20260917/identity/all_cell_identity_evidence.tsv.gz'),('later',P/'geo_identity_20260919/tables/all_cell_identity_evidence.tsv.gz')]
for stage,ef in sources:
 e=pd.read_csv(ef,sep='\t',low_memory=False);e['cluster']=e.cluster.astype(str)
 for sm,g in e.groupby('sample',sort=True):
  if stage=='early':src=P/'OS_CAF_phase2/results'/sm/'annotated.h5ad';inp=O/'inputs'/sm;dc=O/'decontx'/sm;lincol='lineage_phase3'
  elif stage=='discovery':src=P/'phase5_nonpdx_qc_20260916/results'/sm/'annotated.h5ad';inp=P/'phase7_review_states_20260917/remediation_20260917/decontx_inputs'/sm;dc=P/'phase7_review_states_20260917/remediation_20260917/decontx'/sm;lincol='lineage_reviewed'
  else:src=P/'geo_preprocess_20260918/accepted'/sm/'analysis_ready.h5ad';inp=P/'geo_identity_20260919/inputs'/sm;dc=P/'geo_identity_20260919/decontx'/sm;lincol='lineage_phase3'
  a=ad.read_h5ad(src);assert set(a.obs_names)==set(g.cell_id);g=g.set_index('cell_id').loc[a.obs_names];x=a.layers['counts'].tocsr();sy=a.var.symbol.astype(str).to_numpy();raw,lib=evidence(x,sy);raw.index=a.obs_names
  assert np.allclose(lib,g.n_counts.to_numpy())
  ds=str(g.dataset.iloc[0]);technology=str(g['technology'].iloc[0] if 'technology' in g else g.technology_reported.iloc[0]);lineage=g[lincol].astype(str)
  o=pd.DataFrame(index=a.obs_names);o['dataset']=ds;o['sample']=sm;o['sample_key']=ds+'|'+sm;o['patient_id']=g.patient_id.astype(str);o['site']=g.site.astype(str);o['technology']=technology;o['study_technology']=ds+'|'+technology;o['original_lineage']=lineage;o['source_stage']=stage;o['source_file']=str(src);o['old_joint_member']=[i in oldids for i in a.obs_names];o['original_cluster']=g.cluster
  o['donor_known']=~o.patient_id.isin(['UNCONFIRMED','unknown','','nan','NA']);o['donor_record']=np.where(o.donor_known,o.dataset+'|'+o.patient_id,o.sample_key+'|UNRESOLVED');o['inferential_eligible']=o.donor_known&(o.dataset!='GSE250015');o['full_library_counts']=lib;o['n_genes']=np.asarray((x>0).sum(1)).ravel()
  if stage=='later':ic=g.infer_call;cc=g.copy_call;ref=boolcol(g.reference_QC_pass);cta=boolcol(g.CTA_context);o['cnv2_method']='copykat'
  else:
   ic=g.infercnv_call;cc=g.decision_cnv2_call.replace({'normal':'diploid','tumor':'aneuploid'});o['cnv2_method']=g.decision_cnv2_source
   ref=boolcol(g.reference_QC_pass) if stage=='discovery' else boolcol(g.infercnv_reference_check)&boolcol(g.decision_cnv2_reference_check)
   cta=boolcol(g.CTA_context) if 'CTA_context' in g else boolcol(g.cluster_CTA_program)
  o['infer_call']=ic;o['cnv2_call']=cc;o['dual_low']=(ic=='low_signal')&(cc=='diploid');o['reference_QC_pass']=ref;o['CTA_context']=cta
  both=(ic=='elevated_signal')&(cc=='aneuploid');ev=ic.isin(['low_signal','elevated_signal','ambiguous_signal'])&cc.isin(['diploid','aneuploid']);binary=ic.isin(['low_signal','elevated_signal'])&cc.isin(['diploid','aneuploid'])
  o['cluster_both_high_fraction']=np.nan;o['cluster_any_high_fraction']=np.nan
  for cl,idx in g.groupby('cluster',observed=True).groups.items():
   z=ev.loc[idx];z2=binary.loc[idx];o.loc[idx,'cluster_both_high_fraction']=float(both.loc[idx][z].mean()) if z.any() else np.nan;o.loc[idx,'cluster_any_high_fraction']=float(((ic.loc[idx]=='elevated_signal')|(cc.loc[idx]=='aneuploid'))[z2].mean()) if z2.any() else np.nan
  present=(dc/'COMPLETE').exists();cor=None;corlib=np.zeros(len(o));corr=pd.DataFrame(np.nan,index=o.index,columns=raw.columns)
  if present:
   ids=pd.read_csv(inp/'cells.tsv',sep='\t').cell_id;vg=pd.read_csv(inp/'genes.tsv',sep='\t');gene_col='gene_id' if 'gene_id' in vg else 'ensembl_id';assert set(ids)==set(a.obs_names)
   with gzip.open(dc/'corrected_counts.mtx.gz','rb') as f:cm=sp.csr_matrix(mmread(f).T)
   ix=pd.Index(ids).get_indexer(a.obs_names);vi=pd.Index(vg[gene_col].astype(str)).get_indexer(a.var_names);assert (ix>=0).all() and (vi>=0).all();cor=cm[ix][:,vi];diff=cor-x;assert not diff.nnz or diff.data.max()<1e-5
   corr,corlib=evidence(cor,sy);corr.index=a.obs_names;ct=pd.read_csv(dc/'cell_contamination.tsv',sep='\t').set_index('cell_id').loc[a.obs_names];assert np.allclose(corlib,ct.corrected_UMI);o['contamination']=ct.contamination.to_numpy()
  else:o['contamination']=np.nan
  o['ambient_checked']=present;o['corrected_library_counts']=corlib;o['cnv_role']=g.cnv_role.astype(str)
  for col in raw:o['raw_'+col]=raw[col];o['corrected_'+col]=corr[col]
  scope=lineage.eq('FibLike').to_numpy();o['core_scope']=scope
  fixed=o.dual_low.to_numpy()&ref&~cta&(raw.det_CTA<2).to_numpy()&(o.n_genes>=500).to_numpy()&(raw.det_Fib>=2).to_numpy()&(o.cluster_both_high_fraction.fillna(1)<.2).to_numpy()&(o.cluster_any_high_fraction.fillna(1)<.2).to_numpy()
  def gate(frac=.2,osteo=.2,mural=True,det='det_'):
   c,m=conflict(raw,frac=frac,osteo=osteo,mural=mural);d,mm=conflict(corr,frac=frac,osteo=osteo,mural=mural,det=det)
   return fixed&~c&~m&present&(corr[det+'Fib']>=2).to_numpy()&(corr[det+'CTA']<2).to_numpy()&(corlib>0)&~d&~mm
  molecular=gate();o['unified_strict']=scope&molecular;o['osteo_sensitive']=scope&gate(osteo=.5);o['competitor_sensitive']=scope&gate(frac=.5,osteo=.5);o['trace_count_sensitive']=scope&gate(det='pos_');o['mural_reference']=lineage.isin(['FibLike','MesAmb','Mural','Mural_SMC']).to_numpy()&gate(mural=False)&((raw.det_Mural>=3)&(raw.pct_Mural>.2)).to_numpy()
  o['osteo_conflict_raw']=(raw.det_Osteo>=3)&(raw.pct_Osteo>.2);o['osteo_conflict_corrected']=(corr.det_Osteo>=3)&(corr.pct_Osteo>.2)
  co,mix=conflict(raw);cco,cmix=conflict(corr);reasons={'outside_frozen_FibLike_scope':~scope,'not_dual_CNV_low':~o.dual_low.to_numpy(),'reference_QC_fail':~ref,'CTA_context_or_genes':cta|(raw.det_CTA>=2).to_numpy(),'n_genes_below500':(o.n_genes<500).to_numpy(),'Fib_genes_below2':(raw.det_Fib<2).to_numpy(),'cluster_CNV_risk':(o.cluster_both_high_fraction.fillna(1)>=.2).to_numpy()|(o.cluster_any_high_fraction.fillna(1)>=.2).to_numpy(),'raw_competing_or_osteo':co,'raw_mixed':mix,'ambient_unavailable':np.repeat(not present,len(o)),'corrected_Fib_CTA_or_competing':(corr.det_Fib<2).fillna(True).to_numpy()|(corr.det_CTA>=2).fillna(True).to_numpy()|cco|cmix|(corlib<=0)}
  o['exclusion_reasons']=''
  for reason,mask in reasons.items():o.loc[mask,'exclusion_reasons']+=reason+';'
  assert (o.loc[o.unified_strict,'exclusion_reasons']=='').all()
  for key in ['unified_strict','osteo_sensitive','competitor_sensitive','mural_reference','trace_count_sensitive']:
   for role in ['holdout_immune','holdout_endo']:
    rmask=o.cnv_role.eq(role).to_numpy();mg=molecular if key=='unified_strict' else gate(osteo=.5) if key=='osteo_sensitive' else gate(frac=.5,osteo=.5) if key=='competitor_sensitive' else gate(mural=False) if key=='mural_reference' else gate(det='pos_')
    calibration.append({'sample':sm,'dataset':ds,'rule':key,'role':role,'cells':int(rmask.sum()),'molecular_pass_without_scope_gate':int((mg&rmask).sum()),'note':'descriptive negative-control leakage; not independent stromal truth'})
  # Freeze program scores from full source genes and raw counts for interpretation.
  for name,gs in programs.items():
   vals=np.zeros(len(o));n=0
   for gene in gs:
    ix=np.flatnonzero(sy==gene)
    if not len(ix):continue
    n+=1;v=np.asarray(x[:,ix].sum(1)).ravel();vals+=np.log1p(v/np.maximum(lib,1)*1e4)
   o['program_'+name]=vals/max(n,1);cover.append({'sample':sm,'dataset':ds,'program':name,'measured':n,'requested':len(gs)})
  keep=o[['unified_strict','osteo_sensitive','competitor_sensitive','mural_reference','trace_count_sensitive']].any(axis=1)
  if keep.any():
   km=keep.to_numpy();b=a[km].copy();bo=o.loc[keep].copy()
   if ds=='GSE250015':gm=pd.read_csv(P/'geo_cnv_20260918/inputs'/sm/'genes.tsv',sep='\t').set_index('source_gene_id');ids=gm.ensembl_id.reindex(b.var_names).astype(str).to_numpy()
   else:ids=np.array([v.split('.')[0] for v in b.var_names])
   assert len(set(ids))==len(ids);b.var['source_gene_id']=b.var_names;b.var_names=ids;b.obs=bo
   for col in b.obs.select_dtypes(include='object'):b.obs[col]=b.obs[col].fillna('').astype(str)
   b.obs.index=pd.Index(b.obs.index.to_numpy(dtype=object),dtype=object);b.var.index=pd.Index(b.var.index.to_numpy(dtype=object),dtype=object)
   for frame in [b.obs,b.var]:
    for col in frame:
     if isinstance(frame[col].dtype,pd.StringDtype):frame[col]=frame[col].astype(object)
   b.X=x[km];b.layers['counts']=x[km];b.layers['decontx_counts']=cor[km];b.obsm.clear();b.obsp.clear();b.uns.clear();f=O/'objects'/(sm+'_candidate_union.h5ad');b.write_h5ad(f,compression='gzip');outobjects.append({'path':str(f),'sample':sm,'dataset':ds,'cells':b.n_obs,'genes':b.n_vars})
  row={'sample':sm,'dataset':ds,'source_stage':stage,'technology':technology,'retained_cells':len(o),'ambient_checked':present,'old_joint_cells':int(o.old_joint_member.sum()),'unified_strict':int(o.unified_strict.sum()),'osteo_sensitive':int(o.osteo_sensitive.sum()),'competitor_sensitive':int(o.competitor_sensitive.sum()),'mural_reference':int(o.mural_reference.sum()),'old_retained_in_new':int((o.old_joint_member&o.unified_strict).sum()),'old_removed':int((o.old_joint_member&~o.unified_strict).sum()),'new_vs_old':int((~o.old_joint_member&o.unified_strict).sum())};summary.append(row)
  o.to_csv(O/'tables'/(sm+'_unified_evidence.tsv.gz'),sep='\t',index_label='cell_id');allrows.append(o)
  for reason in reasons:gatecounts.append({'sample':sm,'dataset':ds,'reason':reason,'old_joint_excluded_cells':int((reasons[reason]&o.old_joint_member.to_numpy()).sum())})
  man.append({'sample':sm,'source_file':str(src),'source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'evidence_file':str(ef),'ambient_dir':str(dc) if present else 'unavailable_no_previous_FibLike_scope','full_counts_equal_metadata':True,'corrected_not_above_raw':True if present else None});print('UNIFIED',json.dumps(row),flush=True)
allobs=pd.concat(allrows);assert allobs.index.is_unique;assert int(allobs.old_joint_member.sum())==4300
allobs.to_csv(O/'tables/all_cell_unified_evidence.tsv.gz',sep='\t',index_label='cell_id');pd.DataFrame(summary).to_csv(O/'tables/sample_gate_changes.tsv',sep='\t',index=False);pd.DataFrame(man).to_csv(O/'resources/source_manifest.tsv',sep='\t',index=False);pd.DataFrame(outobjects).to_csv(O/'resources/object_manifest.tsv',sep='\t',index=False);pd.DataFrame(cover).to_csv(O/'tables/program_coverage.tsv',sep='\t',index=False);pd.DataFrame(calibration).to_csv(O/'tables/negative_control_gate_checks.tsv',sep='\t',index=False);pd.DataFrame(gatecounts).to_csv(O/'tables/old_joint_exclusion_reason_counts.tsv',sep='\t',index=False)
status={'status':'UNIFIED_GATE_EVIDENCE_COMPLETE','all_retained_cells':len(allobs),'sample_records':len(summary),'old_joint_cells':4300,**{c:int(allobs[c].sum()) for c in ['unified_strict','osteo_sensitive','competitor_sensitive','mural_reference','trace_count_sensitive']},'old_retained_in_new':int((allobs.old_joint_member&allobs.unified_strict).sum()),'original_labels_unchanged':True}
(O/'STATUS.json').write_text(json.dumps(status,indent=2));(O/'UNIFIED_GATES_COMPLETE').write_text('UNIFIED_GATES_COMPLETE\n');print(json.dumps(status),flush=True)