from caf119_common import *
import anndata as ad,shutil,importlib.metadata
event('01_inputs','started');assert (O/'PLAN_ZH.txt').exists()
cpath=R/'phase111_reference_figure_adaptation_20261021/tables/Fig1i_cell_values.tsv.gz';fpath=R/'phase117_joint_restratification_20261023/tables/frozen_registry.tsv.gz';c=pd.read_csv(cpath,sep='\t',low_memory=False);full=pd.read_csv(fpath,sep='\t',low_memory=False).set_index('cell_id');assert len(c)==6633 and c.cell_id.is_unique
m=c[['cell_id','source','donor','site','sample']].copy();x=full.loc[m.cell_id];assert np.array_equal(m.source.to_numpy(),x.source.to_numpy()) and np.array_equal(m.donor.to_numpy(),x.donor.to_numpy())
for col in ['doublet','library_counts','cohort','original_class','original_state','old_lineage','old_lineage_confidence','truth_lineage','conservative_fibro']:m[col]=x[col].to_numpy()
m['donor_key']=m.source+'::'+m.donor;m['fibro_supported']=m.conservative_fibro.astype(bool)&~m.doublet;m.loc[m.source.eq('NCI'),'fibro_supported']=~m.loc[m.source.eq('NCI'),'doublet']
assert set(m.source)=={'SCPCP000017','SCPCP000023','NCI'} and m.groupby('donor_key').site.nunique().max()==1
assert m.site.value_counts().to_dict()=={'Lung':3845,'Bone':2788};assert int((~m.doublet).sum())==6423
old_members=R/'phase103_pooled_NCI_sensitivities_20261020/tables/old_cell_membership.tsv.gz';ncifile=R/'phase102_NCI_author_CAF_20261020/tables/author_CAF_cell_registry.tsv.gz';om=pd.read_csv(old_members,sep='\t');nn=pd.read_csv(ncifile,sep='\t');assert set(m.cell_id)==set(om.loc[om.cellset.eq('expanded'),'cell_id'])|set(nn.cell_id)
paths=[cpath,fpath,old_members,ncifile];canonical_path=R/'phase58_F1F5_round1_20261004/resources/source_canonical_feature_map.tsv';can=pd.read_csv(canonical_path,sep='\t');paths.append(canonical_path)
presence_path=R/'phase41_full_transcriptome_setup_20261001/objects/sample_feature_presence.npz';presence=sparse.load_npz(presence_path).tocsr();si_path=R/'phase41_full_transcriptome_setup_20261001/tables/feature_presence_sample_index.tsv';si=pd.read_csv(si_path,sep='\t').set_index('sample').presence_row;paths.extend([presence_path,si_path])
objects=[]
for src in ['SCPCP000017','SCPCP000023']:
 objects.append(dict(path=R/'phase41_full_transcriptome_setup_20261001/objects/all_candidate_full_genes.h5ad',source=src,rows=m.index[m.source.eq(src)].to_numpy()))
for sm,g in m[m.source.eq('NCI')].groupby('sample'):objects.append(dict(path=R/'phase95_NCI_expanded_reanalysis_20261016/inputs'/sm/'available_counts.h5ad',source='NCI',rows=g.index.to_numpy()))
shared=None;mapping_records=[]
for o in objects:
 p=o['path'];paths.append(p);a=ad.read_h5ad(p,backed='r');vid=pd.Index(a.var_names)
 if o['source']!='NCI':
  cs=can[can.dataset.eq(o['source'])];assert not cs.symbol.duplicated().any();mp={row.symbol:int(vid.get_indexer([row.feature])[0]) for row in cs.itertuples()};samples=m.loc[o['rows'],'sample'].unique();rr=si.loc[samples].to_numpy();measured=np.asarray(presence[rr].sum(0)).ravel()==len(rr);mp={g:j for g,j in mp.items() if j>=0 and measured[j]}
 else:
  sy=a.var.symbol.astype(str);ok=~sy.duplicated(keep=False)&sy.ne('')&~sy.str.startswith('ENSG');mp={s:i for i,s in enumerate(sy) if ok.iloc[i]}
 o['map']=mp;shared=set(mp) if shared is None else shared&set(mp);a.file.close();mapping_records.append(dict(path=str(p),source=o['source'],cells=len(o['rows']),usable_symbols=len(mp)))
genes=sorted(shared);assert len(genes)>=10000;parts=[];order=[]
for o in objects:
 a=ad.read_h5ad(o['path'],backed='r');ix=a.obs_names.get_indexer(m.loc[o['rows'],'cell_id']);assert (ix>=0).all();raw=a.X[ix,:].tocsr();den=np.asarray(raw.sum(1)).ravel();assert np.allclose(den,m.loc[o['rows'],'library_counts']);assert np.max(abs(raw.data-np.round(raw.data)),initial=0)<1e-6
 z=raw[:,[o['map'][g] for g in genes]].astype(np.float32).multiply((10000/den)[:,None]).tocsr();z.data=np.log1p(z.data);parts.append(z);order.extend(o['rows']);a.file.close();print('READ',o['source'],len(o['rows']),flush=True)
xx=sparse.vstack(parts,format='csr');ix=pd.Index(order).get_indexer(np.arange(len(m)));xx=xx[ix];assert xx.shape==(6633,len(genes));sparse.save_npz(O/'objects/log1p_cp10k.npz',xx);save(m,'cohort_registry.tsv.gz');save(pd.DataFrame({'symbol':genes}),'genes.tsv');save(m.groupby(['site','source']).agg(cells=('cell_id','size'),patients=('donor','nunique'),software_positive=('doublet','sum')).reset_index(),'cohort_counts.tsv')
program_path=R/'phase115_incremental_function_20261023/resources/PROGRAMS.json';programs=json.loads(program_path.read_text());paths.append(program_path);geneidx={g:i for i,g in enumerate(genes)};score=m[['cell_id','source','donor','site']].copy();coverage=[]
for name,gs in programs.items():
 idx=[geneidx[g] for g in gs if g in geneidx];valid=len(idx)>=min(3,len(gs)) and len(idx)/len(gs)>=.8;score['program__'+name]=np.asarray(xx[:,idx].mean(1)).ravel() if valid else np.nan;coverage.append(dict(endpoint='program__'+name,defined=len(gs),measured=len(idx),coverage=len(idx)/len(gs),valid=valid,missing=';'.join(set(gs)-set(genes))))
rankpath=R/'phase118_reference_based_broad_classes_20261023/tables/cell_rank_scores_1500.tsv.gz';rank=pd.read_csv(rankpath,sep='\t').set_index('cell_id').loc[m.cell_id];panels=json.loads((R/'phase118_reference_based_broad_classes_20261023/resources/PANELS_FROZEN.json').read_text());paths.append(rankpath)
for name in panels:score['rank__'+name]=rank[name].to_numpy()
dictionary_path=R/'phase62_F1F5_internal_20261005/tables/F1F5_subgroup_summary.tsv';dictionary=pd.read_csv(dictionary_path,sep='\t').fillna('');paths.append(dictionary_path);save(dictionary,'original33_dictionary.tsv')
for row in dictionary.itertuples():
 gs=(row.stable_markers or row.display_markers).split(';');idx=[geneidx[g] for g in gs if g in geneidx];valid=len(idx)>=3 and len(idx)/len(gs)>=.8;name='signature__'+row.subgroup_id;score[name]=np.asarray(xx[:,idx].mean(1)).ravel() if valid else np.nan;coverage.append(dict(endpoint=name,defined=len(gs),measured=len(idx),coverage=len(idx)/len(gs),valid=valid,missing=';'.join(set(gs)-set(genes))))
targets=['COL1A1','COL1A2','COL3A1','COL5A1','COL5A2','COL6A3','POSTN','CD248','FN1','MMP11','CTHRC1','LOX','LOXL2','PCOLCE','SERPINH1','P4HA1','PLOD2','ACTA2','TAGLN','MYH11','RGS5','OMD','OGN','RUNX2','SP7','ALPL','IBSP','PHOSPHO1','PTPRC','PECAM1','CD36','FABP4']
for g in targets:score['gene__'+g]=xx[:,geneidx[g]].toarray().ravel() if g in geneidx else np.nan
save(score,'cell_function_scores.tsv.gz');save(pd.DataFrame(coverage),'endpoint_coverage.tsv');dump('resources/PROGRAMS.json',programs);dump('resources/OBJECTS.json',mapping_records);dump('resources/INPUT_HASHES.json',[dict(path=str(p),sha256=sha(p)) for p in dict.fromkeys(paths)])
dump('PROTOCOL_V1.json',dict(created_at=now(),plan_sha256=sha(O/'PLAN_ZH.txt'),all_cells=6633,fit_cells=6423,sites=['Bone','Lung'],source_gene_intersection=len(genes),canonical_per_source=True,NCC_included=False,normal_reference_in_clustering=False,no_repeat_overall_Bone_Lung_DE=True,NCI_premerge_feature_masks='unavailable; available author object coverage only',parent_resolution=.6,child_resolution=.5,child_eligibility=dict(cells=80,patient_records=3,patients_with5cells=2),candidate_statistic='Patient-level within-site/source contrast; child versus siblings primary, site context separately',figure_queue='paused_by_user'))
dump('resources/SOFTWARE.json',{k:importlib.metadata.version(k) for k in ['numpy','pandas','scipy','scanpy','anndata','harmonypy','leidenalg']})
dump('INPUTS_COMPLETE.json',dict(status='PASS',cells=len(m),fit_cells=int((~m.doublet).sum()),genes=len(genes),cohort_ID_union_exact=True,donor_site_labels_unchanged=True,normal_cells_added=False,NCC_added=False));event('01_inputs','completed',genes=len(genes))
