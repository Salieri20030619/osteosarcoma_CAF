from pathlib import Path
import pandas as pd,numpy as np,json,subprocess,hashlib
from scipy import sparse
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase61_F3_internal_portraits_20261005';A=P/'phase41_full_transcriptome_setup_20261001'
def read(p):return pd.read_csv(p,sep='\t',low_memory=False)
main=read(O/'tables/primary_245_cells.tsv.gz');allcells=read(O/'tables/F3_cells.tsv.gz');assert len(main)==245 and len(allcells)==251 and main.cell_id.is_unique
zf=np.load(O/'objects/primary_full_transcriptome.npz');z={k:zf[k] for k in zf.files};zf.close();assert main.cell_id.tolist()==z['cells'].tolist();assert np.allclose(main.UMI,z['raw_library']) and np.allclose(main.corrected_library,z['corrected_library'])
# Confirm that source-wide canonical selection is measured in each of the three samples.
import anndata as ad
a=ad.read_h5ad(A/'objects/all_candidate_full_genes.h5ad',backed='r');jj=a.var_names.get_indexer(z['features']);a.file.close();presence=sparse.load_npz(A/'objects/sample_feature_presence.npz').tocsr();si=read(A/'tables/feature_presence_sample_index.tsv').set_index('sample').presence_row.to_dict();assert all(np.asarray(presence[si[s],jj].toarray()).all() for s in main['sample'].unique()),'Non-common measured features require comparison-specific masking'
registry=read(O/'tables/contrast_registry.tsv');checks=[]
for row in registry.itertuples(index=False):
 name=row.contrast;mode=row.matrix;allowed=np.ones(len(main),bool)
 if 'Fibro' in name:allowed&=main.lineage_stratum.eq('Fibro').to_numpy()
 if 'same_patient' in name:allowed&=main.donor_record.eq('SCPCP000017|17-439_00067').to_numpy();left=allowed&main.state.eq('C1').to_numpy();right=allowed&main.state.eq('C2').to_numpy()
 else:st=name.split('_')[0];left=allowed&main.state.eq(st).to_numpy();right=allowed&main.state.ne(st).to_numpy()
 assert left.sum()==row.target_cells and right.sum()==row.reference_cells
 frame=read(O/'tables/genes'/f'{row.context}.tsv.gz');lc=z[mode][left].sum(0)/z[mode+'_library'][left].sum()*1e4;rc=z[mode][right].sum(0)/z[mode+'_library'][right].sum()*1e4;effect=np.log2((lc+.1)/(rc+.1));assert np.allclose(effect,frame.mean_delta,atol=1e-10);assert frame.symbol.tolist()==z['symbols'].tolist()
 # Descriptive AUC checked by direct pairwise comparisons for deterministic probe genes.
 ix=np.random.default_rng(20261005).choice(len(frame),20,replace=False)
 for j in ix:
  l=z[mode][left,j]/z[mode+'_library'][left];r=z[mode][right,j]/z[mode+'_library'][right];auc=np.mean((l[:,None]>r[None,:])+.5*(l[:,None]==r[None,:]));assert abs(auc-frame.descriptive_AUC.iloc[j])<1e-9
 rank=read(O/'resources'/f'{row.context}__rank.tsv');assert rank.symbol.is_unique and len(rank)==int(frame.any_expression.sum());checks.append({'context':row.context,'features':len(frame),'rank_genes':len(rank),'max_effect_error':float(np.max(np.abs(effect-frame.mean_delta))),'AUC_direct_checks':20})
pd.DataFrame(checks).to_csv(O/'QA/expression_reconstruction.tsv',sep='\t',index=False);(O/'PREPARE_QA_COMPLETE').write_text('SUCCESS\n')
if not (O/'GSEA_COMPLETE').exists():print('PREPARE_QA_COMPLETE_GSEA_PENDING',flush=True);raise SystemExit(0)
assert (O/'gsea.exit_code').read_text().strip()=='0';libs=json.loads((O/'resources/libraries_canonical.json').read_text());jobs=json.loads((O/'resources/GSEA_JOBS.json').read_text());rows=[];allq=[]
for i,j in enumerate(jobs):
 assert (O/'logs'/f'gsea_{i:02d}.exit_code').read_text().strip()=='0';d=read(Path(j['output']));r=read(Path(j['rank']));bg=set(r.symbol);assert d.Term.is_unique;assert d.GSEA_model_P.between(0,1).all() and np.isfinite(d.NES).all()
 for term,lead in zip(d.Term,d.Lead_genes.fillna('')):
  genes=set(str(lead).split(';'))-{''};assert genes<=bg and genes<=set(libs[j['library']][term])
 rows.append({'context':j['context'],'library':j['library'],'terms':len(d),'rank_genes':len(r),'leading_edge_valid':True});q=d[['Term','GSEA_model_P','model_BH_FDR']].copy();q['family']=j['context']+'|'+j['library'];allq.append(q)
pd.concat(allq).to_csv(O/'QA/GSEA_FDR_input.tsv',sep='\t',index=False)
rscript="x<-read.delim('"+str(O/ 'QA/GSEA_FDR_input.tsv')+"',check.names=FALSE); err<-0;for(f in unique(x$family)){g<-x[x$family==f,];err<-max(err,max(abs(p.adjust(g$GSEA_model_P,'BH')-g$model_BH_FDR)))};stopifnot(err<1e-10);write.table(data.frame(rows=nrow(x),families=length(unique(x$family)),max_FDR_error=err),'"+str(O/'QA/independent_R_FDR.tsv')+"',sep='\\t',row.names=FALSE,quote=FALSE)"
(O/'QA/check_FDR.R').write_text(rscript)
with (O/'QA/independent_R.log').open('w') as f:r=subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript','--vanilla',str(O/'QA/check_FDR.R')],stdout=f,stderr=subprocess.STDOUT)
(O/'QA/independent_R.exit_code').write_text(str(r.returncode));assert r.returncode==0
qa={'full_transcriptome_contrasts':len(checks),'common_measurement_all3samples':True,'all_gene_effects_reconstructed':True,'AUC_pairwise_checks':20*len(checks),'GSEA_jobs':len(rows),'GSEA_rows':sum(x['terms'] for x in rows),'leading_edges_verified':True,'R_FDR':read(O/'QA/independent_R_FDR.tsv').to_dict('records'),'biological_patient_P_generated':False};(O/'delivery/QA_ACCEPTANCE.json').write_text(json.dumps(qa,indent=2));pd.DataFrame(rows).to_csv(O/'QA/GSEA_output_QA.tsv',sep='\t',index=False);(O/'QA_COMPLETE').write_text('SUCCESS\n');print(json.dumps(qa),flush=True)
