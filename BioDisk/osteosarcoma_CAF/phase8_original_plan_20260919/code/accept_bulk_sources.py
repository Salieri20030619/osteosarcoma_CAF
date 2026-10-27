from pathlib import Path
import gzip,io,json,hashlib,os,signal,ctypes,select
import pandas as pd
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase8_original_plan_20260919')
p=O/'bulk/GPL10295_family.soft.gz.part'
if not p.exists():p=O/'bulk/GPL10295_family.soft.gz'
lines=[];inside=False;complete=False
with gzip.open(p,'rt') as f:
 for line in f:
  if line.startswith('!platform_table_begin'):inside=True;continue
  if line.startswith('!platform_table_end'):complete=True;break
  if inside:lines.append(line)
assert complete
t=pd.read_csv(io.StringIO(''.join(lines)),sep='\t',dtype=str)
assert len(t)==48701 and t.ID.is_unique and {'ID','Symbol'}<=set(t.columns)
dest=O/'bulk/GPL10295_platform_table.tsv';dest.write_text(''.join(lines))
targets=[]
want=str(O/'code/prepare.py')
for q in Path('/proc').iterdir():
 if not q.name.isdigit():continue
 try:args=(q/'cmdline').read_bytes().decode().split('\0')
 except (OSError,UnicodeDecodeError):continue
 if want in args and 'python' in Path(args[0]).name:targets.append(int(q.name))
assert len(targets)<=1,targets
for pid in targets:
 libc=ctypes.CDLL(None,use_errno=True);fd=libc.syscall(434,pid,0)
 print('STOP_OWN_REDUNDANT_FAMILY_DOWNLOAD',pid,flush=True);os.kill(pid,signal.SIGTERM)
 if fd>=0:select.select([fd],[],[]);os.close(fd)
# Matrix is a complete previously verified gzip, not a partial source.
matrix=O/'bulk/GSE21257_series_matrix.txt.gz'
with gzip.open(matrix,'rb') as f:
 while f.read(1024*1024):pass
meta={'matrix_complete_gzip':True,'platform_table_complete':True,'platform_rows':len(t),'expected_platform_rows':48701,'annotation_strategy':'complete delimited platform section extracted from official family SOFT; subsequent sample records not required','source_family_download_finished':p.suffix=='.gz','extra_family_download_stopped':bool(targets),'source_partial_retained':str(p),'platform_table_sha256':hashlib.sha256(dest.read_bytes()).hexdigest(),'source_url':'https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL10nnn/GPL10295/soft/GPL10295_family.soft.gz'}
(O/'bulk/data_acceptance.json').write_text(json.dumps(meta,indent=2));(O/'bulk/DOWNLOAD_COMPLETE').write_text('REQUIRED_MATRIX_AND_COMPLETE_PLATFORM_TABLE_ACCEPTED\n')
print(json.dumps(meta),flush=True)
