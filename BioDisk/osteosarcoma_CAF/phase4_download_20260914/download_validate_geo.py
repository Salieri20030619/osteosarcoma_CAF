"""SYNTAX_CHECK: python -m py_compile download_validate_geo.py
RUN_COMMAND: python download_validate_geo.py --root /mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase4_download_20260914
EXPECTED_SUCCESS_MARKER: NON_PDX_DOWNLOAD_VALIDATION_COMPLETE
"""
import os
os.environ.setdefault('OMP_NUM_THREADS','2')
os.environ.setdefault('OPENBLAS_NUM_THREADS','2')
from pathlib import Path
import argparse,json,csv,gzip,tarfile,shutil,hashlib,subprocess,sys,time,platform
from concurrent.futures import ThreadPoolExecutor
import numpy as np
from scipy.io import mmread

def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()

def lines(p):
 op=gzip.open if p.suffix=='.gz' else open
 with op(p,'rt',encoding='utf-8-sig') as f:return [x.rstrip('\r\n').split('\t') for x in f if x.strip()]

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);a=ap.parse_args();root=a.root.resolve()
 (root/'downloads').mkdir(parents=True,exist_ok=True);(root/'matrices').mkdir(exist_ok=True)
 manifest=json.loads((root/'download_manifest.json').read_text())
 assert len(manifest)==17
 assert all(x['direct_patient_verified'] and x['is_xenograft_verified'] is False for x in manifest)
 (root/'python.pid').write_text(str(os.getpid()))
 (root/'environment.json').write_text(json.dumps({'python':sys.version,'executable':sys.executable,'platform':platform.platform(),'numpy':np.__version__},indent=2))
 jobs=[(x,u) for x in manifest for u in x['urls']]
 def download(job):
  x,url=job;p=root/'downloads'/x['dataset']/url.rsplit('/',1)[-1];p.parent.mkdir(exist_ok=True)
  side=p.with_name(p.name+'.verified.json')
  if p.exists() and side.exists() and json.loads(side.read_text())['sha256']==sha(p):return json.loads(side.read_text())
  part=p.with_name(p.name+'.geo.part')
  transport='https://www.ncbi.nlm.nih.gov/geo/download/?acc='+p.name.split('_')[0]+'&format=file&file='+p.name
  print('DOWNLOAD_START',p.name,flush=True)
  subprocess.run(['curl','--fail','--location','--retry','5','--retry-delay','5','--connect-timeout','30','--speed-limit','1024','--speed-time','180','--output',str(part),'--silent','--show-error',transport],check=True)
  assert part.stat().st_size>0
  # Consume gzip/tar streams to verify integrity before considering a download complete.
  if url.endswith('.tar.gz'):
   with tarfile.open(part,'r:gz') as t:
    for m in t:
     if m.isfile():
      with t.extractfile(m) as f:
       while f.read(4*1024*1024):pass
  elif url.endswith('.gz'):
   with gzip.open(part,'rb') as f:
    while f.read(4*1024*1024):pass
  part.replace(p)
  rec={'dataset':x['dataset'],'sample':x['sample'],'file':str(p),'url':url,'transport_url':transport,'bytes':p.stat().st_size,'sha256':sha(p),'status':'download_integrity_pass'}
  side.write_text(json.dumps(rec,indent=2));print('DOWNLOAD_OK',p.name,rec['bytes'],flush=True);return rec
 with ThreadPoolExecutor(max_workers=3) as pool:download_results=list(pool.map(download,jobs))
 (root/'download_checks.json').write_text(json.dumps(download_results,indent=2))
 summary=[]
 for x in manifest:
  folder=root/'matrices'/x['dataset']/x['sample'];folder.mkdir(parents=True,exist_ok=True)
  for url in x['urls']:
   p=root/'downloads'/x['dataset']/url.rsplit('/',1)[-1]
   if p.name.endswith('.tar.gz'):
    with tarfile.open(p,'r:gz') as t:
     for m in t:
      if not m.isfile():continue
      name=Path(m.name).name
      if not any(k in name for k in ['matrix.mtx','genes.tsv','features.tsv','barcodes.tsv']):continue
      dst=folder/name
      if dst.exists():assert dst.stat().st_size==m.size
      else:
       with t.extractfile(m) as src,dst.open('wb') as target:shutil.copyfileobj(src,target)
   else:
    dst=folder/p.name
    if not dst.exists():os.link(p,dst)
  fs=list(folder.iterdir())
  matrix=[p for p in fs if 'matrix.mtx' in p.name];genes=[p for p in fs if 'genes.tsv' in p.name or 'features.tsv' in p.name];barcodes=[p for p in fs if 'barcodes.tsv' in p.name]
  assert len(matrix)==len(genes)==len(barcodes)==1,(x,fs)
  g=lines(genes[0]);bc=lines(barcodes[0]);m=mmread(str(matrix[0])).tocsr()
  assert m.shape==(len(g),len(bc)),(x,m.shape,len(g),len(bc))
  assert np.isfinite(m.data).all() and (m.data>=0).all() and np.equal(m.data,np.floor(m.data)).all()
  assert len({r[0] for r in bc})==len(bc)
  human=sum(r[0].startswith('ENSG') for r in g);mouse=sum(r[0].startswith('ENSMUSG') for r in g)
  assert mouse==0,('MOUSE_GENE_IDS_SOURCE_REVIEW_REQUIRED',x['sample'])
  counts=np.asarray(m.sum(axis=0)).ravel();ng=np.asarray((m>0).sum(axis=0)).ravel()
  rec={k:x[k] for k in ['dataset','sample','patient_id','direct_patient_verified','is_xenograft_verified','clinical_group','source_evidence_url']}
  rec.update({'genes':m.shape[0],'input_barcodes':m.shape[1],'nnz':m.nnz,'zero_count_cells':int((counts==0).sum()),'median_umi':float(np.median(counts)),'median_genes':float(np.median(ng)),'human_ensembl_ids':human,'mouse_ensembl_ids':mouse,'duplicate_gene_ids':len(g)-len({r[0] for r in g}),'duplicate_gene_symbols':len(g)-len({r[1] if len(r)>1 else r[0] for r in g}),'matrix_dir':str(folder),'validation':'PASS','eligible_for_QC':True,'CAF_confirmed':False})
  summary.append(rec);print('MATRIX_VALIDATED',x['sample'],m.shape,flush=True)
  (root/'matrix_validation.json').write_text(json.dumps(summary,indent=2))
  del m
 with (root/'source_and_matrix_validation.tsv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(summary[0]),delimiter='\t');w.writeheader();w.writerows(summary)
 assert len(summary)==17
 (root/'COMPLETE').write_text('NON_PDX_DOWNLOAD_VALIDATION_COMPLETE\n')
 print('NON_PDX_DOWNLOAD_VALIDATION_COMPLETE',len(summary),sum(r['input_barcodes'] for r in summary),flush=True)
if __name__=='__main__':main()
