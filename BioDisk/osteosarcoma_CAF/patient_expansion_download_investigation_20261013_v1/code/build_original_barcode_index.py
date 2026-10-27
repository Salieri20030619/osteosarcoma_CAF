
from pathlib import Path
import pandas as pd,h5py,re,json,gzip,hashlib
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'patient_expansion_download_investigation_20261013_v1';D=O/'resources/source_barcode_index';D.mkdir(exist_ok=True)
m=pd.read_csv(P/'phase26_client_reanalysis_20260922/tables/analysis_source_manifest.tsv',sep='\t');cm=pd.read_csv(P/'phase81_allCAF_bone_lung_20261012/tables/clinical_sample_manifest.tsv',sep='\t').set_index('sample');pat=re.compile(r'([ACGT]{16})');rows=[]
for r in m.itertuples():
 f=P/'phase26_client_reanalysis_20260922/inputs'/r.sample/'full_called_counts.h5ad'
 with h5py.File(f,'r') as h:
  obs=h['obs'];v=obs[obs.attrs['_index']].asstr()[:]
 bc=[pat.search(s).group(1) if pat.search(s) else None for s in v];bb=sorted(set(b for b in bc if b));out=D/(r.sample+'.txt.gz')
 with gzip.open(out,'wt') as stream:stream.write('\n'.join(bb)+'\n')
 clinical=cm.loc[r.sample]
 rows.append(dict(sample=r.sample,dataset=r.dataset,donor_record=clinical.donor_record,site=clinical.strict_site_group,cells=len(v),parsed_barcode_cells=sum(x is not None for x in bc),unique_barcodes=len(bb),barcode_list=str(out),barcode_sha256=hashlib.sha256(('\n'.join(bb)+'\n').encode()).hexdigest()))
pd.DataFrame(rows).to_csv(O/'tables/ORIGINAL_SOURCE_BARCODE_INDEX.tsv',sep='\t',index=False);(O/'resources/ORIGINAL_BARCODE_INDEX_COMPLETE.json').write_text(json.dumps({'samples':len(rows),'status':'complete','cells':sum(r['cells'] for r in rows)},indent=2));print('ORIGINAL_BARCODE_INDEX',len(rows),sum(r['cells'] for r in rows))
