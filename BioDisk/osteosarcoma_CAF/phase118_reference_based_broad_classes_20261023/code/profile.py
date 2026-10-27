from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
import numpy as np,pandas as pd,anndata as ad,json,hashlib,datetime,shutil
from scipy import sparse
from scipy.stats import rankdata
R=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=R/'phase118_reference_based_broad_classes_20261023';B=R/'phase117_joint_restratification_20261023'
for d in ['code','resources','tables','objects','logs','references','delivery']:(O/d).mkdir(parents=True,exist_ok=True)
def dump(n,x):(O/n).write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str)+'\n',encoding='utf-8')
def save(d,n):d.to_csv(O/'tables'/n,sep='\t',index=False,encoding='utf-8')
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
panels={
 'doc_iCAF':['APOD','IGF1','ITM2A','CXCL12','IL6'],
 'doc_myCAF':['ACTA2','MYL9','TAGLN','LRRC15'],
 'doc_mCAF':['MMP11','POSTN','COL14A1','COL11A1','LOX'],
 'doc_osteogenic':['OMD','OGN','SOSTDC1','PHOSPHO1','ALPL'],
 'paper_CD36':['CD36','FABP4','APOE','STEAP4','APOLD1'],
 'paper_MYH11':['MYH11','RGS5','MCAM','MYH9','TGFBI'],
 'paper_C7':['C7','PTGDS','TNXB','FGL2','IGF2','CXCL12'],
 'paper_FNDC1':['FNDC1','COMP','MMP2','COL11A1','COL10A1','P4HA3'],
 'iCAF_extended':['APOD','IGF1','ITM2A','PDGFRA','CXCL12','IL6','DPT','CFD','C3','C7','ANGPTL7'],
 'contractile_extended':['ACTA2','MYL9','TAGLN','MYH11','TPM2','CNN1','DES'],
 'matrix_extended':['MMP11','POSTN','COL1A1','COL1A2','COL3A1','COL5A2','COL8A1','COL10A1','COL11A1','COL12A1','COMP','CTHRC1','INHBA','SULF1','VCAN','LRRC15','ASPN','COL14A1','LOX'],
 'osteogenic_extended':['OMD','OGN','SOSTDC1','PHOSPHO1','ALPL','SP7','RUNX2','IBSP','BGLAP'],
 'lipid_extended':['CD36','FABP4','APOE','STEAP4','APOLD1','PPARG','ANGPTL4','LPL','PLIN1','ADIPOQ'],
 'fibroblast_context':['DCN','LUM','COL1A1','COL1A2','PDGFRA','PDPN','FAP'],
 'mural_context':['RGS5','CSPG4','NOTCH3','MCAM','MYH11','ACTA2','PDGFRB','KCNJ8','ABCC9','DES'],
 'immune_context':['PTPRC','LST1','TYROBP','FCER1G','CD3D','CD3E','CD79A','MS4A1','GNLY','NKG7'],
 'endothelial_context':['PECAM1','VWF','KDR','CDH5','ESAM','EMCN'],
 'epithelial_context':['EPCAM','KRT8','KRT18','KRT19'],
 'erythroid_context':['HBA1','HBA2','HBB','BPGM','ALAS2','GYPA'],
 'osteogenic_lineage_context':['RUNX2','SP7','SATB2','ALPL','PTH1R','IBSP','BGLAP'],
 'chondrocyte_context':['COL2A1','ACAN','SOX9','COMP','CHAD','COL9A1','COL11A2'],
 'myogenic_context':['MYOD1','MYOG','ACTA1','MYH3','DES','MYF6'],
 'cycle_state':['MKI67','TOP2A','UBE2C','CENPF','STMN1','TYMS','TK1','PCNA'],
 'stress_state':['FOS','JUN','JUNB','ATF3','HSPA1A','HSPA1B','DDIT3']}
dump('resources/PANELS_FROZEN.json',panels)
protocol=dict(created_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),scope='Literature-guided broad functional annotation of frozen22 joint J clusters; not reclustering',source_frameworks='User DOCX section5.1 exact four core signatures; user PDF Fig1 and p6 exact four subtype markers. These are two distinct frameworks, not one shared four-class taxonomy.',primary='Four DOCX axes, paper four signatures as comparison; CD36 lipid axis separately assessed; mixed/unresolved and competing lineage retained',scores='Published UCell >=2.7.6 rank formula; Python implementation using average ties, maxRank1500; sensitivity1000/2000; full available gene universe per source object, absent features NaN not zero',score_formula='1-(sum(min(rank,maxRank))-n(n+1)/2)/(n*maxRank-n(n+1)/2), clipped[0,1]',normalization='Rank scores within each cell; log1pCP10k and detection fractions as independent descriptive readouts; full-library denominator',aggregation='Means within source/donor/cluster, then donor-equal within source, source-equal for overall cluster profiles',annotation_order='Examine marker/coexpression, source-stratified and balanced scores and existing identity evidence; freeze named broad classes before new inferential P values. No argmax forced CAF identity.',identity='Osteogenic expression does not prove osteogenic CAF; normal CNV does not prove fibroblast and abnormal/unknown CNV is not silently replaced. Mural/myogenic, immune and endothelial signals retained separately.',statistics='Later source-stratified patient/site paired mean differences; P<0.05 nominal; BH per source/arm/program_or_gene; incomplete biological replication not assigned P',no_PD1_response_inference=True,no_actual_matrix_deposition_claim=True,figure_queue='paused_by_user',old_F_and_J_labels_unchanged=True,planned_outputs=['22J-to-new broad class annotation','per-source patient evidence','oldF/old33/newclass correspondence','source document provenance','internal report'])
dump('PROTOCOL_V1.json',protocol)
m=pd.read_csv(B/'tables/frozen_registry.tsv.gz',sep='\t',low_memory=False);labels=pd.read_csv(B/'tables/joint_cell_labels.tsv.gz',sep='\t').set_index('cell_id');m['joint_cluster']=m.cell_id.map(labels.joint_cluster);assert len(m)==10510 and m.cell_id.is_unique
target=sorted(set().union(*map(set,panels.values())));ti={g:i for i,g in enumerate(target)};rank=np.full((len(m),len(target)),np.nan);rawval=rank.copy();logs=rank.copy();track=[];coverage=[];checks=[]
canonical_path=R/'phase58_F1F5_round1_20261004/resources/source_canonical_feature_map.tsv';canonical=pd.read_csv(canonical_path,sep='\t');canon={s:dict(zip(z.symbol,z.feature)) for s,z in canonical.groupby('dataset')};track.append(canonical_path)
presence_path=R/'phase41_full_transcriptome_setup_20261001/objects/sample_feature_presence.npz';presence=sparse.load_npz(presence_path).tocsr();si_path=R/'phase41_full_transcriptome_setup_20261001/tables/feature_presence_sample_index.tsv';si=pd.read_csv(si_path,sep='\t').set_index('sample').presence_row;track.extend([presence_path,si_path])
objects=[];oldpath=R/'phase41_full_transcriptome_setup_20261001/objects/all_candidate_full_genes.h5ad';objects.append((oldpath,m.index[m.cohort.eq('baseline')].to_numpy(),'baseline'))
for src in ['NCC','NCI']:
 for sm,z in m[m.source.eq(src)].groupby('sample'):
  path=R/'phase93_expanded_admission_20261016/objects'/(sm+'_admitted_raw_counts.h5ad') if src=='NCC' else R/'phase95_NCI_expanded_reanalysis_20261016/inputs'/sm/'available_counts.h5ad'
  objects.append((path,z.index.to_numpy(),src))
for path,rows,kind in objects:
 track.append(path);a=ad.read_h5ad(path,backed='r');ix=a.obs_names.get_indexer(m.loc[rows,'cell_id']);assert (ix>=0).all();x=a.X[ix,:].tocsr();den=np.asarray(x.sum(1)).ravel();assert np.allclose(den,m.loc[rows,'library_counts'])
 sy=a.var.symbol.astype(str).to_numpy();vid=pd.Index(a.var_names)
 # Resolve each target against the existing source-specific canonical map; use exact unique symbol otherwise.
 mapped={}
 for src in m.loc[rows,'source'].unique():
  smap={}
  for gene in target:
   fid=canon.get(src,{}).get(gene);j=vid.get_indexer([fid])[0] if fid is not None else -1
   if j<0:
    hits=np.flatnonzero(sy==gene);j=int(hits[0]) if len(hits)==1 else -1
   smap[gene]=j
  mapped[src]=smap
 for k,idx in enumerate(rows):
  src=m.loc[idx,'source'];sm=m.loc[idx,'sample'];cols=x.indices[x.indptr[k]:x.indptr[k+1]];values=x.data[x.indptr[k]:x.indptr[k+1]];rr=rankdata(-values,method='average');rmap=dict(zip(cols,rr));vmap=dict(zip(cols,values))
  for gene,j in mapped[src].items():
   measured=j>=0 and (kind!='baseline' or bool(presence[int(si.loc[sm]),j]))
   if not measured:continue
   t=ti[gene];v=float(vmap.get(j,0));rawval[idx,t]=v;logs[idx,t]=np.log1p(v/den[k]*1e4);rank[idx,t]=float(rmap.get(j,1e9))
  if k<2:
   dense=x[k].toarray().ravel();dense_ranks=rankdata(-dense,method='average')
   for gene,j in mapped[src].items():
    t=ti[gene]
    if j>=0 and np.isfinite(rank[idx,t]):checks.append(abs(min(dense_ranks[j],2000)-min(rank[idx,t],2000)))
 for sm,z in m.loc[rows].groupby('sample'):
  for gene in target:
   t=ti[gene];meas=int(np.isfinite(rawval[z.index,t]).sum());coverage.append(dict(source=z.source.iloc[0],sample=sm,gene=gene,cells=len(z),measured_cells=meas,feature_index=mapped[z.source.iloc[0]][gene]))
 a.file.close();print('TARGET_EXTRACTION',kind,path.parent.name,len(rows),flush=True)
assert max(checks,default=0)<1e-10
save(pd.DataFrame(coverage),'target_feature_coverage.tsv');np.savez_compressed(O/'objects/target_gene_values.npz',cell_id=m.cell_id.to_numpy(str),genes=np.array(target),raw=rawval,logCP10k=logs,ranks=rank)
scoreframes=[];pcoverage=[]
for maxrank in [1500,1000,2000]:
 frame=m[['cell_id','source','donor','site','sample','cohort','joint_cluster','doublet','old_lineage','old_lineage_confidence']].copy()
 for name,gs in panels.items():
  gidx=[ti[g] for g in gs];r=rank[:,gidx];n=np.isfinite(r).sum(1);rs=np.nansum(np.minimum(r,maxrank),axis=1);u=rs-n*(n+1)/2;den=n*maxrank-n*(n+1)/2;v=np.divide(u,den,out=np.full(len(m),np.nan),where=den>0);v=np.clip(1-v,0,1);valid=n>=max(3,int(np.ceil(.8*len(gs))));v[~valid]=np.nan;frame[name]=v
  if maxrank==1500:
   for src,z in m.groupby('source'):pcoverage.append(dict(source=src,panel=name,defined=len(gs),min_measured=int(n[z.index].min()),max_measured=int(n[z.index].max()),valid_cells=int(valid[z.index].sum()),cells=len(z)))
 save(frame,'cell_rank_scores_'+str(maxrank)+'.tsv.gz');scoreframes.append(frame)
save(pd.DataFrame(pcoverage),'panel_coverage.tsv');score=scoreframes[0]
for frame,maxrank in zip(scoreframes,[1500,1000,2000]):
 patient=frame[~frame.doublet].groupby(['source','donor','joint_cluster'])[list(panels)].mean().reset_index();source=patient.groupby(['source','joint_cluster'])[list(panels)].mean().reset_index();balanced=source.groupby('joint_cluster')[list(panels)].mean().reset_index();save(patient,'patient_panel_scores_'+str(maxrank)+'.tsv');save(source,'source_cluster_scores_'+str(maxrank)+'.tsv');save(balanced,'balanced_cluster_scores_'+str(maxrank)+'.tsv')
# Patient-equal expression and RNA detection for each individual referenced marker.
long=[]
for (src,donor,cl),z in m[~m.doublet].groupby(['source','donor','joint_cluster']):
 for j,gene in enumerate(target):
  good=np.isfinite(rawval[z.index,j]);idx=z.index[good];long.append(dict(source=src,donor=donor,joint_cluster=cl,gene=gene,cells=len(z),measured_cells=len(idx),mean_logCP10k=float(logs[idx,j].mean()) if len(idx) else np.nan,detection_fraction=float((rawval[idx,j]>0).mean()) if len(idx) else np.nan))
long=pd.DataFrame(long);save(long,'patient_marker_expression.tsv');ss=long.groupby(['source','joint_cluster','gene'])[['mean_logCP10k','detection_fraction']].mean().reset_index();save(ss,'source_marker_expression.tsv');save(ss.groupby(['joint_cluster','gene'])[['mean_logCP10k','detection_fraction']].mean().reset_index(),'balanced_marker_expression.tsv')
oldmeta=pd.read_csv(R/'phase66_doublet_repair_DSC1_20261005/tables/cells_repaired_with_confidence.tsv.gz',sep='\t',low_memory=False);save(m[~m.doublet].groupby(['joint_cluster','old_lineage','old_lineage_confidence']).size().rename('cells').reset_index(),'inherited_lineage_composition.tsv')
save(m,'cell_registry.tsv.gz');hashes=[dict(path=str(p),sha256=sha(p)) for p in dict.fromkeys(track)];dump('resources/INPUT_HASHES.json',hashes)
dump('PROFILE_COMPLETE.json',dict(status='PASS_PROFILE_READY_FOR_ANNOTATION',cells=len(m),fit_cells=int((~m.doublet).sum()),joint_clusters=m.joint_cluster.nunique(),target_genes=len(target),panels=len(panels),dense_rank_checks=len(checks),max_rank_difference=max(checks,default=0),main_maxRank=1500,sensitivities=[1000,2000],old_J_labels_modified=False,no_inferential_P_before_annotation=True,software_UCell_package_called=False,formula_source='https://code.bioconductor.org/browse/UCell/'))
print('BALANCED_PROFILES');print(scoreframes[0][~scoreframes[0].doublet].groupby('joint_cluster')[list(panels)[:8]].mean().round(3).to_string())
