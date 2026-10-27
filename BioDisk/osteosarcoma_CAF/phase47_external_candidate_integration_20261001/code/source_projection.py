import sys
sys.path.insert(0,'/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase46_spatial_functional_20261001/code')
from common import *
E=P/'phase47_external_candidate_integration_20261001';pr=json.loads((O/'resources/FUNCTIONAL_SIGNATURES_FROZEN.json').read_text());cw=pd.read_csv(C/'tables/HGNC_full_gene_symbol_crosswalk.tsv.gz',sep='\t',keep_default_na=False);cw['priority']=~cw.mapping_method.eq('HGNC_Ensembl_approved');cw=cw[cw.canonical_symbol.ne('')].sort_values(['priority','feature']).drop_duplicates('canonical_symbol');sf=dict(zip(cw.canonical_symbol,cw.feature));var=pd.read_csv(A/'tables/gene_id_crosswalk.tsv.gz',sep='\t',nrows=1) if False else None
contexts=json.loads((C/'resources/EFFECT_CONTEXTS.json').read_text());global_features=cw.feature.tolist();import anndata as ad
a=ad.read_h5ad(A/'objects/all_candidate_full_genes.h5ad',backed='r');feature_names=np.array(a.var_names);a.file.close();rows=[];coverage=[]
for ctx in contexts:
 if ctx['contrast'] not in ['F3_vs_rest','F4_vs_rest']:continue
 npz=np.load(C/ctx['effect_file'],allow_pickle=True);effects=npz['effects'];features=feature_names[npz['feature_indices']];gm={f:i for i,f in enumerate(features)};donors=npz['donors'].astype(str)
 ends={('program',k):{'up':v,'down':[]} for k,v in pr['programs'].items()};ends.update({('signature',k):v for k,v in pr['signatures'].items()})
 for (kind,name),sg in ends.items():
  ui=[gm[sf[g]] for g in sg['up'] if g in sf and sf[g] in gm];di=[gm[sf[g]] for g in sg['down'] if g in sf and sf[g] in gm]
  def average(ids):
   if not ids:return np.full(len(donors),np.nan)
   v=effects[:,ids];num=np.isfinite(v).sum(1);return np.divide(np.nansum(v,axis=1),num,out=np.full(len(donors),np.nan),where=num>0)
  v=average(ui)-(average(di) if sg['down'] else 0)
  if sg['down'] and not di:v[:]=np.nan
  for d,value in zip(donors,v):rows.append(dict(patient_id=d,layer=ctx['layer'],dataset=ctx['dataset'],class_id=ctx['target'],kind=kind,endpoint=name,effect=value,discovery_source=ctx['dataset']=='SCPCP000017',independent_external_confirmation=False))
  coverage.append(dict(context_id=ctx['context_id'],kind=kind,endpoint=name,up_requested=len(sg['up']),up_present=len(ui),down_requested=len(sg['down']),down_present=len(di)))
f=pd.DataFrame(rows);f.to_csv(E/'tables/source_projection_patient_effects.tsv',sep='\t',index=False);pd.DataFrame(coverage).to_csv(E/'tables/source_projection_coverage.tsv',sep='\t',index=False);s=test_table(f,['layer','dataset','class_id','kind','endpoint'],'effect',lambda r:'source_projection_'+r.layer+'_'+r.dataset+'_'+r.class_id+'_'+r.kind);s.to_csv(E/'tables/source_projection_statistics.tsv',sep='\t',index=False);(E/'SOURCE_PROJECTION_COMPLETE').write_text('SUCCESS');print('SOURCE_PROJECTION_DONE',len(s),int(s.P.lt(.05).sum()))
