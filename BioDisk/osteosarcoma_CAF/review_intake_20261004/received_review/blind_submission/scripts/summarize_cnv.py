import pathlib,json
import numpy as np,pandas as pd
P=pathlib.Path(__file__).parent; O=P/'review_output'; s=pd.read_pickle(O/'cell_structure.pkl'); lookup=s.set_index('cell_id'); meta=pd.read_csv(P/'sample_technical_metadata.tsv',sep='\t').set_index('sample_id')
cm=pd.read_csv(P/'copykat_cell_continuous_metrics.tsv.gz',sep='\t'); assert not cm.cell_id.duplicated().any(); cm=cm.merge(s[['cell_id','source_id','patient_id','r10']],on='cell_id',validate='one_to_one',how='left'); assert cm.r10.notna().all()
cc=cm.groupby(['source_id','r10','sample_id','patient_id']).agg(copykat_cells=('cell_id','size'),mean_abs_signal_median=('mean_abs_native_signal','median'),mean_abs_signal_mean=('mean_abs_native_signal','mean'),binned_regions_min=('binned_regions','min')).reset_index().rename(columns={'r10':'cluster_id'})
cc.to_csv(O/'copykat_sample_cluster_metrics.tsv',sep='\t',index=False)
profiles=[]
for d in sorted(P.glob('D0?')):
 f=pd.read_csv(d/'copykat_sample_cluster_profiles.tsv.gz',sep='\t'); idc=[c for c in f if not c.startswith('chr')]; bins=[c for c in f if c.startswith('chr')];
 print(d.name,idc,flush=True)
 for chrom in sorted(set(c.split(':')[0] for c in bins)):
  cols=[c for c in bins if c.startswith(chrom+':')]; z=f[idc].copy(); z['chromosome']=chrom; z['native_chromosome_mean']=f[cols].mean(axis=1); z['native_chromosome_mean_abs']=f[cols].abs().mean(axis=1); z['valid_bins']=f[cols].notna().sum(axis=1); z['source_id']=d.name; profiles.append(z)
pd.concat(profiles).to_csv(O/'copykat_chromosome_profiles.tsv.gz',sep='\t',index=False)
regional=[]; regcov=[]; filestatus=[]
for f in sorted((P/'regional_CNV').glob('*.tsv.gz')):
 if f.name.endswith('_regions.tsv.gz'): continue
 sample,kind=f.name.split('__',1); kind=kind.replace('.tsv.gz',''); x=pd.read_csv(f,sep='\t'); sample_s=s[s.sample_id==sample]; merged=x.merge(sample_s[['cell_id','source_id','patient_id','r10']],on='cell_id',how='left',validate='many_to_one'); assert merged.r10.notna().all(),f
 filestatus.append({'sample_id':sample,'kind':kind,'rows':len(x),'unique_cells':x.cell_id.nunique(),'duplicate_cells':int(x.cell_id.duplicated().sum())})
 if kind.endswith('_groups'):
  method=kind[:-7]; r=pd.read_csv(str(f).replace('_groups.tsv.gz','_regions.tsv.gz'),sep='\t'); r['chr']=r.chr.astype(str).str.replace('chr','',regex=False); assert (r.end>=r.start).all();
  r['chromosome']='chr'+r.chr; r['direction']=np.select([r.state<3,r.state>3],['loss','gain'],default='neutral')
  def union_len(g):
   intervals=sorted(zip(g.start,g.end)); total=0; a=b=None
   for u,v in intervals:
    if b is None: a,b=u,v
    elif u<=b: b=max(b,v)
    else: total+=b-a+1;a,b=u,v
   return total+(b-a+1 if b is not None else 0)
  event=[]
  for (group,ch,di),rg in r.groupby(['group_id','chromosome','direction']): event.append({'group_id':group,'chromosome':ch,'direction':di,'reported_region_Mb':union_len(rg)/1e6})
  e=pd.DataFrame(event)
  cnt=merged.groupby(['source_id','patient_id','r10','group_id']).size().rename('group_cells').reset_index()
  total=merged.groupby(['source_id','patient_id','r10']).size().rename('assigned_cells').reset_index(); total['sample_id']=sample;total['method']=method; regcov.append(total)
  if len(e):
   z=cnt.merge(e,on='group_id').merge(total,on=['source_id','patient_id','r10']); z['weighted_Mb']=z.reported_region_Mb*z.group_cells/z.assigned_cells
   z=z.groupby(['source_id','patient_id','r10','sample_id','method','chromosome','direction'],as_index=False).agg(reported_region_Mb=('weighted_Mb','sum')); z['metric']='cell_weighted_reported_region_Mb'; z=z.rename(columns={'reported_region_Mb':'value'}); regional.append(z)
 else:
  method='provider_arm'; total=merged.groupby(['source_id','patient_id','r10']).size().rename('assigned_cells').reset_index(); total['sample_id']=sample;total['method']=method; regcov.append(total)
  cols=[c for c in x if c.startswith('proportion_')]; z=merged.groupby(['source_id','patient_id','r10'])[cols].mean().reset_index().melt(id_vars=['source_id','patient_id','r10'],var_name='measurement',value_name='value'); z['sample_id']=sample;z['method']=method; z['region']=z.measurement.str.split('_').str[-1]; z['chromosome']=z.region.str.replace('[pq]$','',regex=True); z['direction']=z.measurement.str.split('_').str[1].replace({'dupli':'gain','cnv':'any'}); z['metric']='cell_mean_arm_fraction'; z=z.drop(columns='measurement');regional.append(z)
pd.concat(regional).rename(columns={'r10':'cluster_id'}).to_csv(O/'regional_CNV_summary.tsv.gz',sep='\t',index=False)
pd.concat(regcov).rename(columns={'r10':'cluster_id'}).to_csv(O/'regional_CNV_coverage.tsv',sep='\t',index=False)
pd.DataFrame(filestatus).to_csv(O/'regional_file_status.tsv',sep='\t',index=False)
# Primary structure merges and subdivisions retained without using UMAP separation as proof.
trans=[]
for res in ['r05','r20']:
 z=s.groupby(['source_id','r10',res]).size().rename('cells').reset_index(); z['fraction_of_primary']=z.cells/z.groupby('r10').cells.transform('sum'); z['fraction_of_other_resolution']=z.cells/z.groupby(res).cells.transform('sum'); z['resolution']=res;z=z.rename(columns={res:'other_cluster_id','r10':'cluster_id'}); trans.append(z)
pd.concat(trans).to_csv(O/'resolution_transitions.tsv',sep='\t',index=False)
(O/'cnv_completed.json').write_text(json.dumps({'status':'complete','copykat_cells':len(cm),'regional_files':len(filestatus),'cell_roster':len(s)}),encoding='utf-8')
print('COMPLETED',len(cm),len(filestatus),flush=True)
