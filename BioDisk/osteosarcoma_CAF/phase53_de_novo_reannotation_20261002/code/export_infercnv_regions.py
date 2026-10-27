from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,re,gzip
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');R=P/'phase53_de_novo_reannotation_20261002';B=P/'phase26_client_reanalysis_20260922'
m=pd.read_csv(P/'phase41_full_transcriptome_setup_20261001/tables/clinical_sample_registry.tsv',sep='\t',keep_default_na=False)
out=R/'objects'/'infercnv_regions';out.mkdir(parents=True,exist_ok=True)
meta=[];metrics=[]
for row in m.itertuples(index=False):
 sm=row.sample;cells=pd.read_csv(B/'inputs'/sm/'cells.tsv',sep='\t')
 provider=B/'inputs'/sm/'provider_infercnv_table.tsv.gz'
 if provider.exists():
  t=pd.read_csv(provider,sep='\t',float_precision='round_trip')
  t['cell_id']=sm+'_'+t.barcodes.astype(str).str.strip()
  assert t.cell_id.is_unique
  matched=t.cell_id.isin(cells.cell_id);unmatched=int((~matched).sum());t=t[matched].copy()
  cols=[x for x in t if x.startswith('proportion_loss_chr') or x.startswith('proportion_dupli_chr') or x.startswith('proportion_cnv_chr')]
  assert cols
  t[['cell_id']+cols].to_csv(out/f'{sm}__provider_arm_fractions.tsv.gz',sep='\t',index=False)
  numeric=t[cols].to_numpy(float);assert np.isfinite(numeric).all() and (numeric>=0).all() and (numeric<=1+1e-8).all()
  x=t[['cell_id']].copy();x['sample']=sm;x['method_source']='provider_HMM';x['mean_region_fraction']=t[[c for c in cols if c.startswith('proportion_cnv_chr')]].mean(1);x['max_region_fraction']=t[[c for c in cols if c.startswith('proportion_cnv_chr')]].max(1)
  metrics.append(x);meta.append(dict(sample=sm,source='provider_HMM',input_cells=len(cells),covered_cells=len(t),unmatched_provider_rows=unmatched,regions=len(cols),file=str(provider),SHA256=hashlib.sha256(provider.read_bytes()).hexdigest()))
 roots=[('native_HMM',B/'cnv'/'infercnv'/sm)]
 if sm in ['SCPCS000610','SCPCS000614']:
  for variant in ['immune_endothelial_reference','lymphoid_reference']:
   roots.append((variant,P/'phase34_reference_rerun_20260927/cnv_reference'/sm/variant))
 for name,d in roots:
  gf=list(d.glob('17_HMM*.cell_groupings'));rf=list(d.glob('HMM_CNV_predictions*Pnorm_*.pred_cnv_regions.dat'))
  if len(gf)!=1 or len(rf)!=1:
   if name=='native_HMM' and not provider.exists():meta.append(dict(sample=sm,source=name,input_cells=len(cells),covered_cells=0,unmatched_provider_rows=0,regions=0,file='',SHA256=''))
   continue
  groups=pd.read_csv(gf[0],sep='\t');regions=pd.read_csv(rf[0],sep='\t')
  assert {'cell_group_name','cell'}.issubset(groups.columns)
  groups=groups.merge(cells[['tool_id','cell_id']],left_on='cell',right_on='tool_id',how='left')
  assert groups.cell_id.notna().all() and groups.cell_id.is_unique
  assert {'cell_group_name','state','chr','start','end'}.issubset(regions.columns)
  assert regions.cell_group_name.isin(groups.cell_group_name).all()
  groups[['cell_id','cell_group_name']].to_csv(out/f'{sm}__{name}_groups.tsv.gz',sep='\t',index=False)
  regions.to_csv(out/f'{sm}__{name}_regions.tsv.gz',sep='\t',index=False)
  # Region counts and total called segment span are descriptive, not malignancy thresholds.
  alt=regions[regions.state.ne(3)].copy();alt['span']=alt['end']-alt.start+1
  ag=alt.groupby('cell_group_name').agg(altered_segments=('state','size'),altered_span_bp=('span','sum'))
  x=groups[['cell_id','cell_group_name']].merge(ag,on='cell_group_name',how='left').fillna({'altered_segments':0,'altered_span_bp':0})
  x['sample']=sm;x['method_source']=name;metrics.append(x.drop(columns='cell_group_name'))
  meta.append(dict(sample=sm,source=name,input_cells=len(cells),covered_cells=len(groups),unmatched_provider_rows=0,regions=len(regions),file=str(rf[0]),SHA256=hashlib.sha256(rf[0].read_bytes()).hexdigest()))
 print('REGIONS_COMPLETE',sm,flush=True)
pd.DataFrame(meta).to_csv(R/'tables/infercnv_region_coverage.tsv',sep='\t',index=False)
pd.concat(metrics,ignore_index=True).to_csv(R/'tables/infercnv_cell_region_metrics.tsv.gz',sep='\t',index=False)
(R/'INFER_CNV_REGIONS_COMPLETE').write_text('SUCCESS')
print(json.dumps({'samples':len(m),'source_records':len(meta),'covered_records':sum(x['covered_cells'] for x in meta),'native_states_neutral':3,'no_new_malignancy_threshold':True}),flush=True)
