from common import *
import os,select,subprocess,anndata as ad
Q=O/'independent_QA';Q.mkdir(exist_ok=True)
for stage in ['gsea_pipeline','correction','characterize']:
 if not (O/(stage+'.exit_code')).exists():
  try:subprocess.run(['/usr/bin/python3','-c','import os,select,sys; f=os.pidfd_open(int(sys.argv[1])); select.select([f],[],[]); os.close(f)',(O/(stage+'.pid')).read_text().strip()],check=True)
  except ProcessLookupError:pass
 assert (O/(stage+'.exit_code')).read_text().strip()=='0'
a=ad.read_h5ad(A/'objects/all_candidate_full_genes.h5ad',backed='r');features=a.var_names.astype(str).to_numpy();a.file.close();expected=[];effects=[];rng=np.random.default_rng(20261004)
# Prespecified target genes and a reproducible random sample of feature-level P values per contrast; all retained in output.
for f in sorted((O/'tables/genes').glob('*.tsv.gz')):
 x=read(f);corrected='corrected' in f.name
 for contrast,g in x.groupby('contrast',sort=False):
  ds=g.dataset.iloc[0];ctx=f'expanded__{ds}__{contrast}';ef=O/'objects'/(ctx+('__corrected_patient_effects.npz' if corrected else '__patient_effects.npz'))
  if not ef.exists():continue
  z=np.load(ef);d=z['effects'];donors=z['donors'].astype(str);fe=features[z['feature_indices']];idx={v:i for i,v in enumerate(fe)};chosen=set(g.loc[g.symbol.isin(['COL1A1','LOX','MMP14','RUNX2','PLXDC2','EXT1','ITGB1','MYH11','MT-CO2','NCAM1']),'feature']);chosen.update(rng.choice(g.feature,size=min(16,len(g)),replace=False));chosen.update(g.sort_values('P').head(4).feature)
  tag=('corrected' if corrected else 'raw')+'|'+ctx
  for row in g[g.feature.isin(chosen)].itertuples(index=False):
   item='gene|'+row.feature;expected.append({'context_id':tag,'item':item,'P':row.P,'N':int(row.N),'type':'gene_P_sample'})
   for j,donor in enumerate(donors):effects.append({'context_id':tag,'item':item,'donor_record':donor,'effect':d[j,idx[row.feature]]})
# All Hallmark/Matrisome/customer programs plus reproducibly sampled Reactome programs, for every evaluable context.
for file,corrected in [(O/'tables/all_program_statistics.tsv.gz',False),(O/'tables/corrected_program_statistics.tsv.gz',True)]:
 x=read(file)
 for (ds,contrast),g in x.groupby(['dataset','contrast'],sort=False):
  ctx=f'expanded__{ds}__{contrast}';ef=O/'objects'/(ctx+('__corrected_program_effects.npz' if corrected else '__program_effects.npz'))
  if not ef.exists():continue
  z=np.load(ef);names=z['programs'].astype(str);donors=z['donors'].astype(str);d=z['effects'];idx={v:i for i,v in enumerate(names)};chosen=set(n for n in names if not n.startswith('reactome|'));other=[n for n in names if n.startswith('reactome|')];chosen.update(rng.choice(other,size=min(20,len(other)),replace=False));tag=('corrected' if corrected else 'raw')+'|'+ctx
  for row in g[g.program.isin(chosen)].itertuples(index=False):
   item='program|'+row.program;expected.append({'context_id':tag,'item':item,'P':row.P,'N':int(row.N),'type':'program_P_sample'})
   for j,donor in enumerate(donors):effects.append({'context_id':tag,'item':item,'donor_record':donor,'effect':d[j,idx[row.program]]})
s=read(O/'tables/surface_competition_statistics.tsv');e=read(O/'tables/surface_competition_patient_effects.tsv.gz')
for row in s.itertuples(index=False):
 tag='surface|'+row.dataset+'|'+row.class_id+'|'+row.background;expected.append({'context_id':tag,'item':row.gene,'P':row.P,'N':int(row.N),'type':'all_surface_P'});z=e[(e.dataset==row.dataset)&(e.class_id==row.class_id)&(e.background==row.background)&(e.gene==row.gene)]
 for r in z.itertuples(index=False):effects.append({'context_id':tag,'item':row.gene,'donor_record':r.donor_record,'effect':r.effect})
ex=pd.DataFrame(expected);ex.to_csv(Q/'P_expected.tsv',sep='\t',index=False);pd.DataFrame(effects).to_csv(Q/'P_patient_effects.tsv.gz',sep='\t',index=False)
for n in sorted(ex.N.unique()):
 if n>15:np.savetxt(Q/f'MC_signs_N{n}.tsv.gz',np.random.default_rng(20261001).choice([-1.,1.],size=(19999,n)),fmt='%d',delimiter='\t')
# GSEA rank and leading-edge provenance checks, separate from patient tests.
libs=json.loads((O/'resources/libraries_canonical.json').read_text());gle=[]
for f in (O/'tables/gsea').glob('*.tsv.gz'):
 x=read(f);parts=f.name.replace('__GSEA.tsv.gz','').split('__');ctx='__'.join(parts[:-1]);library=parts[-1];rank=read(O/'resources'/f'{ctx}__rank.tsv');bg=set(rank.symbol);assert rank.symbol.is_unique
 assert pd.to_numeric(x['GSEA_model_P']).between(0,1).all();assert np.isfinite(pd.to_numeric(x.NES)).all();assert x.Term.is_unique
 for row in x.itertuples(index=False,name=None):pass
 for term,lead in zip(x.Term,x.Lead_genes.fillna('')):
  genes=set(str(lead).split(';'))-{''};assert genes<=bg and genes<=set(libs[library][term]),(f,term)
 gle.append({'file':str(f),'rows':len(x),'ranking_genes':len(rank),'leading_edge_gene_membership_valid':True,'P_range_valid':True,'finite_NES':True})
pd.DataFrame(gle).to_csv(Q/'GSEA_input_output_QA.tsv',sep='\t',index=False)
(O/'resources/INDEPENDENT_QA_SCOPE.json').write_text(json.dumps({'P_checks':ex.type.value_counts().to_dict(),'selection_seed':20261004,'full_FDR_reconstruction':True,'GSEA_all_outputs_rank_and_leading_edge_checked':True},indent=2))
env=os.environ.copy();env.update(OMP_NUM_THREADS='4',OPENBLAS_NUM_THREADS='4',MKL_NUM_THREADS='4')
with (Q/'independent_R.log').open('w') as f:r=subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript','--vanilla',str(O/'code/independent_QA.R')],stdout=f,stderr=subprocess.STDOUT,env=env)
(Q/'independent_R.exit_code').write_text(str(r.returncode));assert r.returncode==0;(O/'INDEPENDENT_QA_COMPLETE').write_text('SUCCESS\n');print('INDEPENDENT_QA_COMPLETE',flush=True)