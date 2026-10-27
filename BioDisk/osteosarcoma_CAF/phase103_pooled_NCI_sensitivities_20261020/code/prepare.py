import pathlib,json,hashlib,datetime
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
R=pathlib.Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');P=R/'phase103_pooled_NCI_sensitivities_20261020'
for d in ['code','resources','tables','report','figures','logs']: (P/d).mkdir(parents=True,exist_ok=True)
B=R/'phase81_allCAF_bone_lung_20261012';Q=R/'phase73_F134_bone_lung_20261009/strict_v2'
genes=json.loads((R/'phase96_all14_literature_20261016/project_registry.json').read_text())
meta_file=Q/'tables/expanded_fixed_cells_with_sites.tsv.gz';ob=pd.read_csv(meta_file,sep='\t')
roster=[];expected={}
for src in ['SCPCP000017','SCPCP000023']:
 f=B/'objects'/('all_cells__expanded__'+src+'__All_CAF__independent.npz');a=np.load(f,allow_pickle=False)
 for site in ['Bone','Lung']:
  for i,donor in enumerate(a[site.lower()+'_patients']):
   roster.append(dict(source=src,donor=donor,site=site))
   for g in genes:expected[(donor,site,g['gene'])]=float(a[site.lower()+'_expression'][i,np.flatnonzero(a['features']==g['feature_ids'])[0]])
roster=pd.DataFrame(roster);assert len(roster)==22 and roster.donor.is_unique
ob=ob.merge(roster,left_on=['dataset','donor_record','strict_site'],right_on=['source','donor','site'],validate='many_to_one',suffixes=('_legacy',''))
assert ob.cell_id.is_unique
full=ad.read_h5ad(R/'phase41_full_transcriptome_setup_20261001/objects/all_candidate_full_genes.h5ad',backed='r')
ix=full.obs_names.get_indexer(ob.cell_id);assert (ix>=0).all()
gi=full.var_names.get_indexer([g['feature_ids'] for g in genes]);assert (gi>=0).all()
# Read each target column to avoid backed sparse fancy-index ambiguity.
target=np.column_stack([np.asarray(full.X[:,int(j)].toarray()).ravel()[ix] for j in gi]);full.file.close()
presence=sparse.load_npz(R/'phase41_full_transcriptome_setup_20261001/objects/sample_feature_presence.npz').tocsr()
sample_index=pd.read_csv(R/'phase41_full_transcriptome_setup_20261001/tables/feature_presence_sample_index.tsv',sep='\t').set_index('sample').presence_row
truth=lambda s:s.astype(str).str.lower().isin(['true','1','1.0'])
masks={'expanded':np.ones(len(ob),bool),'singlets':~truth(ob.doublet_positive_any).to_numpy(),'conservative':truth(ob.conservative_fibro).to_numpy() & ~truth(ob.doublet_positive_any).to_numpy()}
strict_donors=set(ob.loc[masks['conservative'],'donor'])
masks['expanded_matched_strict_donors']=ob.donor.isin(strict_donors).to_numpy()
rows=[];members=[]
for label,mask in masks.items():
 z=ob.loc[mask].copy()
 for (src,donor,site),c in z.groupby(['source','donor','site']):
  ids=c.index.to_numpy();raw=target[ids,:].sum(0);den=float(c.full_library_counts.sum())
  valid=np.asarray(presence[sample_index.loc[c['sample'].unique()].to_numpy()][:,gi].sum(0)).ravel()==c['sample'].nunique()
  assert valid.all()
  expr=np.log2(raw/den*1e4+.1)
  for j,g in enumerate(genes):
   rows.append(dict(cellset=label,source=src,donor=donor,site=site,gene=g['gene'],feature=g['feature_ids'],cells=len(c),raw_count=raw[j],library_UMIs=den,expression=expr[j]))
   if label=='expanded':assert abs(expr[j]-expected[(donor,site,g['gene'])])<1e-12
  members.append(c[['cell_id','source','donor','site','sample']].assign(cellset=label))
nci=pd.read_csv(R/'phase102_NCI_author_CAF_20261020/tables/author_CAF_patient_expression.tsv',sep='\t')
for label in masks:
 st='author_singlets' if label in ['singlets','conservative'] else 'author_all'
 n=nci[nci.stratum.eq(st)].copy()
 for _,r in n.iterrows():rows.append(dict(cellset=label,source='NCI',donor=r.patient_id,site='Lung',gene=r.gene,feature=r.feature,cells=r.cells,raw_count=r.raw_count,library_UMIs=r.library_UMIs,expression=r.expression))
out=pd.DataFrame(rows)
out.to_csv(P/'tables/all_patient_expression.tsv',sep='\t',encoding='utf-8',index=False)
pd.concat(members).to_csv(P/'tables/old_cell_membership.tsv.gz',sep='\t',encoding='utf-8',index=False)
coverage=out[out.gene.eq('EBF1')].groupby(['cellset','source','site']).agg(patients=('donor','nunique'),cells=('cells','sum')).reset_index()
coverage.to_csv(P/'tables/cohort_coverage.tsv',sep='\t',encoding='utf-8',index=False)
ids=out[out.cellset.eq('expanded') & out.gene.eq('EBF1')][['source','donor','site','cells']].copy()
ids['identifier_without_source']=ids.donor.str.split('|').str[-1].str.replace('NCI_PATIENT_','',regex=False)
assert not ids.identifier_without_source.duplicated().any()
ids.to_csv(P/'tables/patient_identity_audit.tsv',sep='\t',encoding='utf-8',index=False)
d=dict(status='PASS',original_expanded_reproduced_values=len(expected),original_full_frozen_cells=8174,main_old_cells=int(coverage[coverage.cellset.eq('expanded')&coverage.source.ne('NCI')].cells.sum()),main_NCI_cells=1530,main_bone_patients=16,main_lung_patients=23,strict_donors=len(strict_donors),no_shared_public_patient_identifiers=True,genotype_or_global_author_crosswalk_available=False,old_paired_donors_kept_excluded=True)
(P/'resources/EXTRACTION_QA.json').write_text(json.dumps(d,indent=2),encoding='utf-8')
files=[meta_file,R/'phase102_NCI_author_CAF_20261020/tables/author_CAF_patient_expression.tsv']
(P/'resources/INPUT_HASHES.json').write_text(json.dumps([dict(path=str(f),sha256=hashlib.sha256(f.read_bytes()).hexdigest()) for f in files],indent=2),encoding='utf-8')
print(coverage.to_string(index=False));print(json.dumps(d))
