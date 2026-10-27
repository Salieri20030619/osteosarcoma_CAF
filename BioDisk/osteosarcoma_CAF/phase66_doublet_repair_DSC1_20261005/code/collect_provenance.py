from pathlib import Path
import pandas as pd,anndata as ad,json,datetime
from current_metadata import load_current_cells
B=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=B/'phase66_doublet_repair_DSC1_20261005'
ob=load_current_cells();out=[];inputs=[]
for sample,t in ob.groupby('sample'):
 p=B/'phase26_client_reanalysis_20260922/inputs'/sample/'full_called_counts.h5ad'
 a=ad.read_h5ad(p,backed='r');obs=a.obs.loc[t.cell_id].copy();a.file.close()
 cols=[c for c in obs if any(v in c.lower() for v in ['doublet','scrublet','scdblfinder'])]
 d=obs[cols].copy().add_prefix('input_');d.index.name='cell_id';d['source_sample']=sample;d['source_file']=str(p)
 native=B/'phase31_customer_completion_20260926/doublet'/sample
 if (native/'COMPLETE').is_file() and (native/'exit_code').read_text().strip()=='0':
  e=pd.read_csv(native/'doublet_results.tsv',sep='\t').set_index('cell_id').loc[t.cell_id]
  for c in e:d['native_'+c]=e[c]
  d['native_source_file']=str(native/'doublet_results.tsv')
 out.append(d.reset_index());inputs.append({'sample':sample,'source':str(p),'cells':len(t),'fields':';'.join(cols)})
d=pd.concat(out,ignore_index=True);assert len(d)==8174 and d.cell_id.is_unique
d.to_csv(O/'tables/doublet_method_details.tsv.gz',sep='\t',index=False);pd.DataFrame(inputs).to_csv(O/'resources/DOUBLET_METHOD_INPUTS.tsv',sep='\t',index=False)
hits=[]
for phase in ['phase58_F1F5_round1_20261004','phase59_four_functional_axes_20261004','phase62_F1F5_internal_20261005']:
 for p in (B/phase/'code').glob('*.py'):
  ss=p.read_text().splitlines()
  for i,line in enumerate(ss):
   if any(k in line for k in ['doublet_assessed_any','target_cell_doublet_evidence','cell_review_decisions']):
    hits.append({'file':str(p),'line':i+1,'context':'\n'.join(ss[max(0,i-3):i+5])})
(O/'resources/METADATA_INHERITANCE_AUDIT.json').write_text(json.dumps(hits,indent=2))
(O/'PROVENANCE_COMPLETE.json').write_text(json.dumps({'success':True,'cells':len(d),'samples':len(inputs),'at':datetime.datetime.now(datetime.timezone.utc).isoformat()},indent=2))
print('PROVENANCE_COMPLETE',len(d),len(inputs));print(json.dumps(hits,ensure_ascii=False)[:8000])
