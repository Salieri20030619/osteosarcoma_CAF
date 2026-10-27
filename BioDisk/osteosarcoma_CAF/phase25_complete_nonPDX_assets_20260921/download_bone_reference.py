from pathlib import Path
import subprocess,json,hashlib,tarfile,h5py,shutil,sys,datetime
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase25_complete_nonPDX_assets_20260921/references/GSE217792';O.mkdir(parents=True,exist_ok=True)
url='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217792/suppl/GSE217792_RAW.tar';p=O/'GSE217792_RAW.tar';part=O/'GSE217792_RAW.tar.part'
with (O/'download.log').open('a') as log:
 rc=subprocess.run(['curl','-q','--fail','--location','--silent','--show-error','--connect-timeout','30','--speed-limit','1024','--speed-time','180','--retry','3','--retry-delay','10','--continue-at','-','--proxy','http://127.0.0.1:7897','--noproxy','','--output',str(part),url],stdout=log,stderr=subprocess.STDOUT).returncode
(O/'download.exit_code').write_text(str(rc)+'\n')
if rc:sys.exit(rc)
rows=[]
with tarfile.open(part) as t:
 for m in t.getmembers():
  if not m.isfile():continue
  target=(O/'matrices'/m.name).resolve();target.relative_to((O/'matrices').resolve());target.parent.mkdir(parents=True,exist_ok=True)
  with t.extractfile(m) as src,target.open('wb') as dst:shutil.copyfileobj(src,dst)
  if target.suffix=='.h5':
   with h5py.File(target) as f:
    group=f['matrix'] if 'matrix' in f else next(v for v in f.values() if isinstance(v,h5py.Group) and 'shape' in v)
    shape=group['shape'][:];rows.append({'file':str(target),'genes':int(shape[0]),'barcodes':int(shape[1]),'nonzero_entries':len(group['data']),'h5_readable':True})
assert len(rows)==9,'expected 9 official sample matrices'
part.replace(p)
h=hashlib.sha256()
with p.open('rb') as f:
 for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
q={'source':url,'sha256':h.hexdigest(),'bytes':p.stat().st_size,'samples':rows,'note':'9 direct-human-bone reference samples include 4 controls and 5 nonunions; do not treat all 9 as healthy controls','formal_analysis_started':False,'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
(O/'QA.json').write_text(json.dumps(q,indent=2));(O/'COMPLETE').write_text('REFERENCE_MATRICES_READABLE\n');print('GSE217792_REFERENCE_DOWNLOAD_COMPLETE',len(rows))
