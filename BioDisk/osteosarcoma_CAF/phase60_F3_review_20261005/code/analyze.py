from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
import json,hashlib,datetime,gc,itertools,shutil
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
from statsmodels.stats.multitest import multipletests
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');B=P/'phase58_F1F5_round1_20261004';Q=P/'phase59_four_functional_axes_20261004';A=P/'phase41_full_transcriptome_setup_20261001';O=P/'phase60_F3_review_20261005'
def read(p,**kw):return pd.read_csv(p,sep='\t',low_memory=False,**kw)
def dump(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str),encoding='utf-8')
def mean(x,axis=0):
 n=np.isfinite(x).sum(axis);out=np.full(np.shape(n),np.nan);np.divide(np.nansum(x,axis=axis),n,out=out,where=n>0);return out
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
for d in ['code','resources','tables','objects','figures','QA','logs','delivery']:(O/d).mkdir(parents=True,exist_ok=True)
ob=read(Q/'tables/frozen_cells_with_identity.tsv.gz').set_index('cell_id');f3=ob[ob.class_id.eq('F3')];assert len(f3)==251
sources=sorted(f3.dataset.unique());libs=json.loads((B/'resources/libraries_canonical.json').read_text());sets={};origin={}
for key in ['MYC_TARGETS_V1','MYC_TARGETS_V2','OXIDATIVE_PHOSPHORYLATION','G2M_CHECKPOINT','E2F_TARGETS','UNFOLDED_PROTEIN_RESPONSE','APOPTOSIS','HYPOXIA','INTERFERON_GAMMA_RESPONSE','INTERFERON_ALPHA_RESPONSE','INFLAMMATORY_RESPONSE']:
 sets[key]=libs['hallmark']['HALLMARK_'+key];origin[key]='hallmark|HALLMARK_'+key
for key in ['TRANSLATION','EUKARYOTIC_TRANSLATION_INITIATION','EUKARYOTIC_TRANSLATION_ELONGATION','MITOCHONDRIAL_TRANSLATION']:
 sets[key]=libs['reactome']['REACTOME_'+key];origin[key]='reactome|REACTOME_'+key
prior=json.loads((Q/'resources/GENE_SETS.json').read_text())
for k in ['contraction','matrix','osteogenic','mineralization']:sets[k]=prior['definitions'][k];origin[k]=prior['origins'][k]
panels={'MHCII':['HLA-DRA','HLA-DRB1','HLA-DPA1','HLA-DPB1','HLA-DQA1','HLA-DQB1','CD74','CIITA'],'fibroblast':['DCN','LUM','COL1A1','COL1A2','DPT','PDGFRA'],'myeloid':['PTPRC','LST1','TYROBP','FCER1G','LYZ','CD68'],'Bcell':['MS4A1','CD79A','CD79B','CD37'],'T_NK':['CD3D','CD3E','TRAC','NKG7','GNLY'],'osteogenic_markers':['RUNX2','ALPL','SP7','IBSP','BGLAP'],'cycling_markers':['MKI67','TOP2A','PCNA','BIRC5','UBE2C','CENPF'],'costimulation':['CD80','CD86'],'matrix_markers':prior['definitions']['matrix_markers']}
for k,v in panels.items():sets[k]=v;origin[k]='prespecified_diagnostic_panel_not_classifier'
for k in ['MYC_TARGETS_V1','MYC_TARGETS_V2','OXIDATIVE_PHOSPHORYLATION','TRANSLATION']:
 sets[k+'_no_MT_RP']=[g for g in sets[k] if not g.startswith(('MT-','RPL','RPS'))];origin[k+'_no_MT_RP']='original_set_minus_mitochondrial_and_ribosomal_protein_genes'
roles={k:('technical_gene_removal_sensitivity' if k.endswith('_no_MT_RP') else 'diagnostic_panels' if k in panels else 'full_programs') for k in sets};programs=list(sets);genes=sorted(set().union(*map(set,sets.values())));gix={g:i for i,g in enumerate(genes)};mem={k:np.array([gix[g] for g in v]) for k,v in sets.items()};markers=sorted(set().union(*map(set,panels.values())));markerix=[gix[g] for g in markers]
protocol={'at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope':'F3 donor/state, expression lineage, QC, MHCII and translation/metabolism/cycle follow-up; F4 existing two main expression strata exported with residual category; no new cluster fitting','preexisting_findings_known':True,'statistics':'same source/donor/site pseudobulk full-library log2(CP10k+0.1); donor equal after site equal; exact two-sided signflip; N<2 or allzero no P; N<10 exploratory; no cell-level P','families':'source|matrix|context|full_programs/diagnostic_panels/technical_gene_removal_sensitivity; marker genes separate','screening':'raw target P<0.05 nominates statistical candidates; when impossible retain descriptive effects without calling nominal candidates','coverage':'80% and at least min(3,size) measured in all contributing samples; unassayed not zero','primary_comparison':'F3 versus original other four classes; Fibro and Osteo restricted comparisons on both arms; each original F3 internal state versus rest of F3','sensitivity':'minimum10 cells in both matched arms; remove MT-/RPL/RPS genes from selected pathways; same-sample raw-library10th-90thpercentile matched BG-Fibro based on F3-Fibro target','external_background':'raw same donor/site background outside8174 in Fibro/Osteo/Myeloid; no disease-normal interpretation','MHCII_codetection':'at least2 MHCII and2 fibro genes, each count>=1; PTPRC>=1 flagged separately; descriptive only; exclude missing coverage from denominators','F4_groups':'F4_Fibro_like; F4_Osteo_like includes Osteo/Osteochondro/Chondro; F4_other_unresolved retains Mural and all residual; these are analysis strata not newly validated subclusters'}
dump(O/'resources/PROTOCOL.json',protocol);dump(O/'resources/GENE_SETS.json',{'sets':sets,'origins':origin,'roles':roles});dump(O/'resources/PREVIOUS_CURRENT.json',json.loads((P/'CAF_NEXT_STAGE_CURRENT.json').read_text()));pd.DataFrame([{'program':k,'gene':g,'role':roles[k],'origin':origin[k]} for k,v in sets.items() for g in v]).to_csv(O/'resources/gene_set_membership.tsv',sep='\t',index=False)
f4=ob[ob.class_id.eq('F4')].copy();f4['analysis_group']=np.select([f4.lineage_stratum.eq('Fibro'),f4.lineage_stratum.eq('Osteo')],['F4_Fibro_like','F4_Osteo_like'],default='F4_other_unresolved');f4.to_csv(O/'tables/F4_two_main_strata_with_residual.tsv.gz',sep='\t');f4.groupby(['analysis_group','review_lineage_code']).size().rename('cells').reset_index().to_csv(O/'tables/F4_stratum_counts.tsv',sep='\t',index=False)
qc=read(B/'tables/cell_QC_and_class.tsv.gz').set_index('cell_id');state=read(P/'phase42_F3F4_identity_states_20261001/tables/expanded_F3_internal_cell_states.tsv.gz').set_index('cell_id');f3=f3.join(qc[['UMI','detected_features','mitochondrial_fraction','ribosomal_fraction']]);f3=f3.join(state[['UMAP1','UMAP2']]);mp=read(B/'resources/source_canonical_feature_map.tsv');mapping={ds:mp[mp.dataset.eq(ds)].set_index('symbol').feature.to_dict() for ds in sources};presence=sparse.load_npz(A/'objects/sample_feature_presence.npz').tocsr();sidx=read(A/'tables/feature_presence_sample_index.tsv').set_index('sample').presence_row.to_dict();counts={};denom={};assay=np.zeros((len(ob),len(genes)),bool)
for mode,path in [('raw',A/'objects/all_candidate_full_genes.h5ad'),('corrected',B/'objects/expanded_corrected_full_genes.h5ad')]:
 a=ad.read_h5ad(path);assert a.obs_names.tolist()==ob.index.tolist();v=np.zeros((len(ob),len(genes)));dd=np.asarray(a.X.sum(1)).ravel()
 for ds in sources:
  rr=np.flatnonzero(ob.dataset.eq(ds));ii=a.var_names.get_indexer([mapping[ds].get(g,'__missing') for g in genes]);ok=ii>=0;v[np.ix_(rr,np.flatnonzero(ok))]=a.X[rr][:,ii[ok]].toarray()
  if mode=='raw':
   for sample in ob.iloc[rr]['sample'].unique():
    ri=np.flatnonzero(ob['sample'].eq(sample));assay[np.ix_(ri,np.flatnonzero(ok))]=presence[sidx[sample],ii[ok]].toarray().ravel()>0
 counts[mode]=v;denom[mode]=dd;del a;gc.collect()
fi=ob.index.get_indexer(f3.index);f3['corrected_library']=denom['corrected'][fi];f3['fraction_removed_from_library']=1-denom['corrected'][fi]/denom['raw'][fi];f3.to_csv(O/'tables/F3_cells_review_QC.tsv.gz',sep='\t')
f3.groupby(['dataset','sample','donor_record','site','round1_state','lineage_stratum'],observed=True).size().rename('cells').reset_index().to_csv(O/'tables/F3_patient_state_lineage_counts.tsv',sep='\t',index=False)
f3.groupby(['dataset','donor_record','round1_state'],observed=True).agg(cells=('sample','size'),median_UMI=('UMI','median'),median_features=('detected_features','median'),median_MT=('mitochondrial_fraction','median'),median_Ribo=('ribosomal_fraction','median'),median_removed_fraction=('fraction_removed_from_library','median')).reset_index().to_csv(O/'tables/F3_patient_state_QC.tsv',sep='\t',index=False)
np.savez_compressed(O/'objects/diagnostic_counts.npz',raw=counts['raw'],corrected=counts['corrected'],raw_library=denom['raw'],corrected_library=denom['corrected'],assayed=assay,genes=np.array(genes),cells=ob.index.to_numpy(str))
stats=[];geffects=[];effects=[];gstats=[];roster=[];scorevalues=[]
def test(v):
 n=len(v)
 if n<2 or not np.any(v):return np.nan
 assert n<=5
 signs=np.array(list(itertools.product([-1.,1.],repeat=n)));return np.mean(np.abs(signs@v/n)>=abs(v.mean())-1e-12)
def compare(ds,mode,ctx,left,right,minimum=1,bg=None):
 left=np.asarray(left,int);right=np.asarray(right,int);donor={};donorgenes={}
 for (d,site),lg in ob.iloc[left].groupby(['donor_record','site'],observed=True):
  il=ob.index.get_indexer(lg.index)
  if len(il)<minimum:continue
  if bg is None:
   rg=ob.iloc[right];rg=rg[(rg.donor_record==d)&(rg.site==site)];ir=ob.index.get_indexer(rg.index)
   if len(ir)<minimum:continue
   rc=counts[mode][ir].sum(0);rd=denom[mode][ir].sum();ra=assay[ir].all(0);nr=len(ir)
  else:
   ref=bg.get((ds,str(d),str(site)))
   if ref is None or ref['cells']<minimum:continue
   rc=ref['counts'];rd=ref['library'];ra=ref['assayed'];nr=ref['cells']
  ld=denom[mode][il].sum()
  if ld<=0 or rd<=0:continue
  valid=assay[il].all(0)&ra;ll=np.log2(counts[mode][il].sum(0)/ld*1e4+.1);rr=np.log2(rc/rd*1e4+.1);diff=ll-rr;diff[~valid]=np.nan;pv=[]
  for k,ii in mem.items():
   covered=int(valid[ii].sum());ok=covered>=max(min(3,len(ii)),int(np.ceil(.8*len(ii))));value=float(np.nanmean(diff[ii])) if ok else np.nan;pv.append(value);scorevalues.append({'dataset':ds,'matrix':mode,'context':ctx,'donor_record':d,'site':site,'program':k,'target_score':float(ll[ii][valid[ii]].mean()) if ok else np.nan,'reference_score':float(rr[ii][valid[ii]].mean()) if ok else np.nan,'effect':value,'measured_genes':covered,'definition_genes':len(ii)})
  donor.setdefault(str(d),[]).append(pv);donorgenes.setdefault(str(d),[]).append(diff[markerix]);roster.append({'dataset':ds,'matrix':mode,'context':ctx,'donor_record':d,'site':site,'target_cells':len(il),'reference_cells':nr,'minimum_cells':minimum})
 for k,key in enumerate(programs):
  vals=[]
  for d,v in donor.items():
   value=float(mean(np.array(v)[:,k]));vals.append(value);effects.append({'dataset':ds,'matrix':mode,'context':ctx,'program':key,'donor_record':d,'effect':value})
  v=np.array(vals);v=v[np.isfinite(v)];n=len(v);stats.append({'dataset':ds,'matrix':mode,'context':ctx,'program':key,'N':n,'mean_delta':v.mean() if n else np.nan,'positive_patients':int((v>0).sum()),'negative_patients':int((v<0).sum()),'P':test(v),'family':'|'.join([ds,mode,ctx,roles[key]]),'min_possible_two_sided_P':2/2**n if n>=2 else np.nan,'evidence':'exploratory' if n>=2 else 'direction_only_or_unassessable'})
 for k,gene in enumerate(markers):
  vals=[]
  for d,v in donorgenes.items():
   value=float(mean(np.array(v)[:,k]));vals.append(value);geffects.append({'dataset':ds,'matrix':mode,'context':ctx,'gene':gene,'donor_record':d,'effect':value})
  v=np.array(vals);v=v[np.isfinite(v)];gstats.append({'dataset':ds,'matrix':mode,'context':ctx,'gene':gene,'N':len(v),'mean_delta':v.mean() if len(v) else np.nan,'positive_patients':int((v>0).sum()),'negative_patients':int((v<0).sum()),'P':test(v),'family':'|'.join([ds,mode,ctx,'marker_genes'])})
def ix(m):return np.flatnonzero(m)
for ds in sources:
 same=ob.dataset.eq(ds);cl=ob.class_id.eq('F3')
 for mode in counts:
  for context,allowed,minimum in [('F3_vs_rest',np.ones(len(ob),bool),1),('F3_vs_rest_min10',np.ones(len(ob),bool),10),('F3_Fibro_vs_rest_Fibro',ob.lineage_stratum.eq('Fibro'),1),('F3_Osteo_vs_rest_Osteo',ob.lineage_stratum.eq('Osteo'),1)]:compare(ds,mode,context,ix(same&cl&allowed),ix(same&~cl&allowed),minimum)
  for st in ['expanded_F3_C0','expanded_F3_C1','expanded_F3_C2']:compare(ds,mode,st+'_vs_F3_rest',ix(same&cl&ob.round1_state.eq(st)),ix(same&cl&ob.round1_state.ne(st)))
 print('CANDIDATES',ds,flush=True)
# Independent cell-level descriptions; no cell-level inferential tests.
cellrows=[]
for mode in counts:
 for cellid,i in zip(f3.index,fi):
  row={'cell_id':cellid,'matrix':mode}
  for k in panels:
   ii=mem[k];row[k+'_detected_genes']=int((counts[mode][i,ii]>=1).sum());row[k+'_coverage']=float(assay[i,ii].mean())
  row['PTPRC_detected']=bool(counts[mode][i,gix['PTPRC']]>=1);row['MHCII_fibro_codetected']=row['MHCII_detected_genes']>=2 and row['fibroblast_detected_genes']>=2 and row['MHCII_coverage']>=.8 and row['fibroblast_coverage']>=.8;cellrows.append(row)
pd.DataFrame(cellrows).to_csv(O/'tables/F3_cell_diagnostic_codetection.tsv.gz',sep='\t',index=False)
# Focused raw references from same donors outside original layer.
anno=read(P/'phase55_annotation_freeze_20261004/tables/all_cell_annotation_v2.tsv.gz');subset=anno[anno.original_donor_record.isin(f3.donor_record)&~anno.original_cell_id.isin(ob.index)].copy();subset['background']=np.select([subset.review_lineage_code.eq('Fibro'),subset.review_lineage_code.isin(['Osteo','Chondro','Osteochondro']),subset.review_lineage_code.isin(['Macro','Myeloid','Mono','MonoDC','DC','pDC'])],['Fibro','Osteo','Myeloid'],default='Other');subset=subset[subset.background.ne('Other')];bg={k:{} for k in ['Fibro','Osteo','Myeloid','Fibro_UMI_matched']};bgq=[];bgr=[];matchlog=[]
for sampleid,sub in subset.groupby('sample_id'):
 x=ad.read_h5ad(P/'phase53_de_novo_reannotation_20261002/blind_review/raw_counts'/f'{sampleid}.h5ad');ds=sub.dataset.iloc[0];jj=x.var_names.get_indexer([mapping[ds].get(g,'__absent') for g in genes]);ok=jj>=0;target=f3[(f3['sample']==sub['sample'].iloc[0])&f3.lineage_stratum.eq('Fibro')];lo,hi=target.UMI.quantile([.1,.9]).to_numpy() if len(target) else (np.nan,np.nan)
 for (donor,site,group),part in sub.groupby(['original_donor_record','site','background']):
  rows=x.obs_names.get_indexer(part.cell_id);assert (rows>=0).all();groups=[(group,rows)]
  if group=='Fibro':
   libs=x.obs.iloc[rows].library_counts.to_numpy(float);chosen=rows[(libs>=lo)&(libs<=hi)];groups.append(('Fibro_UMI_matched',chosen));matchlog.append({'sample_id':sampleid,'donor_record':donor,'target_F3_Fibro_cells':len(target),'library_lower':lo,'library_upper':hi,'background_original_cells':len(rows),'matched_cells':len(chosen)})
  for name,ii in groups:
   if not len(ii):continue
   xx=x.X[ii];den=float(x.obs.iloc[ii].library_counts.sum());vv=np.zeros(len(genes));vv[ok]=np.asarray(xx[:,jj[ok]].sum(0)).ravel();key=(ds,str(donor),str(site));prev=bg[name].get(key)
   if prev is None:bg[name][key]={'counts':vv,'library':den,'cells':len(ii),'assayed':ok.copy()}
   else:prev['counts']+=vv;prev['library']+=den;prev['cells']+=len(ii);prev['assayed']&=ok
   symbols=x.var['symbol'].astype(str) if 'symbol' in x.var else pd.Series(x.var_names.astype(str));mt=symbols.str.startswith('MT-').to_numpy();lib=x.obs.iloc[ii].library_counts.to_numpy(float);nfeatures=xx.getnnz(axis=1);mtfrac=np.asarray(xx[:,mt].sum(1)).ravel()/np.maximum(lib,1);bgq.append({'sample_id':sampleid,'dataset':ds,'donor_record':donor,'site':site,'background':name,'cells':len(ii),'median_UMI':float(np.median(lib)),'median_features':float(np.median(nfeatures)),'median_MT':float(np.median(mtfrac))});bgr.append({'sample_id':sampleid,'dataset':ds,'donor_record':donor,'site':site,'background':name,'cells':len(ii),'library':den})
 del x;gc.collect();print('BACKGROUND',sampleid,flush=True)
for ds in sources:
 for stratum in ['All','Fibro','Osteo']:
  allowed=ob.dataset.eq(ds)&ob.class_id.eq('F3')
  if stratum!='All':allowed&=ob.lineage_stratum.eq(stratum)
  for name in bg:compare(ds,'raw','F3_'+stratum+'_vs_BG_'+name,ix(allowed),[],1,bg[name])
pd.DataFrame(bgq).to_csv(O/'tables/background_QC.tsv',sep='\t',index=False);pd.DataFrame(bgr).to_csv(O/'tables/background_roster.tsv',sep='\t',index=False);pd.DataFrame(matchlog).to_csv(O/'tables/UMI_matching_ledger.tsv',sep='\t',index=False)
np.savez_compressed(O/'objects/background_pseudobulks.npz',**{k:np.array([{'key':key,**v} for key,v in d.items()],dtype=object) for k,d in bg.items()})
def adjust(d):
 d=d.copy();d['FDR']=np.nan
 for family,g in d.groupby('family'):
  q=multipletests(g.P.fillna(1),method='fdr_bh')[1];d.loc[g.index,'FDR']=np.where(g.P.notna(),q,np.nan)
 d['nominal']=d.P.lt(.05);d['FDR_significant']=d.FDR.lt(.05);return d
adjust(pd.DataFrame(stats)).to_csv(O/'tables/program_statistics.tsv',sep='\t',index=False);adjust(pd.DataFrame(gstats)).to_csv(O/'tables/marker_statistics.tsv',sep='\t',index=False);pd.DataFrame(effects).to_csv(O/'tables/patient_program_effects.tsv',sep='\t',index=False);pd.DataFrame(geffects).to_csv(O/'tables/patient_marker_effects.tsv',sep='\t',index=False);pd.DataFrame(roster).to_csv(O/'tables/patient_site_roster.tsv',sep='\t',index=False);pd.DataFrame(scorevalues).to_csv(O/'tables/patient_site_program_scores.tsv',sep='\t',index=False)
# Reuse full-transcriptome F3 results with original inferential family and exact provenance.
allgene=[]
for ds in sources:
 for mode,suffix in [('raw','all_statistics'),('corrected','corrected_statistics')]:
  old=read(B/'tables/genes'/f'expanded__{ds}__{suffix}.tsv.gz');old=old[old.contrast.eq('F3_vs_rest')].copy();old['matrix']=mode;old['evidence_use']='descriptive_full_gene_effects_not_nominal_candidates';allgene.append(old)
pd.concat(allgene,ignore_index=True).to_csv(O/'tables/reused_F3_full_transcriptome_statistics.tsv.gz',sep='\t',index=False)
reused=[]
for path in (B/'tables/gsea').glob('*__F3_vs_rest__*.tsv.gz'):
 target=O/'tables'/('reused_'+path.name);shutil.copy2(path,target);reused.append({'source':str(path),'copy':str(target.relative_to(O)),'SHA256':sha(path),'interpretation':'existing gene-set permutation model, not patient biological treatment P'})
dump(O/'resources/REUSE_MANIFEST.json',reused);dump(O/'resources/INPUT_MANIFEST.json',[{'path':str(p),'bytes':p.stat().st_size,'SHA256':sha(p)} for p in [A/'objects/all_candidate_full_genes.h5ad',B/'objects/expanded_corrected_full_genes.h5ad',Q/'tables/frozen_cells_with_identity.tsv.gz',P/'phase42_F3F4_identity_states_20261001/tables/expanded_F3_internal_cell_states.tsv.gz']]);dump(O/'delivery/ANALYSIS_ACCEPTANCE.json',{'cells_F3':len(f3),'F3_donor_records':f3.donor_record.nunique(),'F3_source_counts':f3.dataset.value_counts().to_dict(),'F4_groups':f4.analysis_group.value_counts().to_dict(),'programs':len(sets),'genes':len(genes),'program_rows':len(stats),'gene_rows':len(gstats),'new_clustering':False,'new_CNV':False,'QA':'pending'});(O/'ANALYSIS_COMPLETE').write_text('SUCCESS\n');print('ANALYSIS_COMPLETE',flush=True)
