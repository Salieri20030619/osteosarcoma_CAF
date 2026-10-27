import pathlib,subprocess,concurrent.futures,hashlib,json,gzip
p=pathlib.Path(__file__).resolve().parents[1]
tasks=[('bone_msc_verified.rds.gz','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE253nnn/GSE253355/suppl/GSE253355_MSC_Subset_Seurat.rds.gz',401085038),('lung_counts_verified.tsv.gz','https://cells.ucsc.edu/stanford-czb-hlca/droplet/exprMatrix.tsv.gz',132732799),('lung_meta_verified.tsv','https://cells.ucsc.edu/stanford-czb-hlca/droplet/meta.tsv',9616483)]
parts=p/'inputs/parts';mini=p/'inputs/mini';mini.mkdir(exist_ok=True);jobs=[]
for n,u,size in tasks:
 for start in range(0,size,4194304):
  end=min(start+4194303,size-1);f=parts/(n+'.'+str(start))
  if f.exists() and f.stat().st_size==end-start+1:continue
  for a in range(start,end+1,1048576):jobs.append((n,u,a,min(a+1048575,end)))
def fetch(j):
 n,u,s,e=j;f=mini/(n+'.'+str(s))
 if f.exists() and f.stat().st_size==e-s+1:return True
 r=subprocess.run(['curl','-fLsS','--connect-timeout','10','--max-time','120','--retry','12','--retry-all-errors','--retry-delay','3','-r',str(s)+'-'+str(e),'-o',str(f),u],capture_output=True,text=True)
 good=r.returncode==0 and f.exists() and f.stat().st_size==e-s+1
 if not good:print('FAILED',j,r.stderr,flush=True)
 return good
print('small chunks needed',len(jobs),flush=True)
with concurrent.futures.ThreadPoolExecutor(6) as pool:results=list(pool.map(fetch,jobs))
assert all(results),'incomplete download'
for n,u,size in tasks:
 for start in range(0,size,4194304):
  end=min(start+4194303,size-1);f=parts/(n+'.'+str(start))
  if f.exists() and f.stat().st_size==end-start+1:continue
  with f.open('wb') as w:
   for a in range(start,end+1,1048576):w.write((mini/(n+'.'+str(a))).read_bytes())
manifest=[]
for n,u,size in tasks:
 f=p/'inputs'/n
 with f.open('wb') as w:
  for s in range(0,size,4194304):w.write((parts/(n+'.'+str(s))).read_bytes())
 assert f.stat().st_size==size
 if n.endswith('.gz'):
  with gzip.open(f,'rb') as z:
   while z.read(8388608):pass
 manifest.append(dict(file=str(f),url=u,bytes=size,sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
(p/'resources/VERIFIED_DOWNLOADS.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
print('DOWNLOADS_VERIFIED',flush=True)
