from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='2'
import requests,hashlib,json,base64,zlib,struct,time
import cellxgene_census,tiledbsoma
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase122_HLCA_supplement_20261026')
url='https://cellxgene-census-public-us-west-2.s3.us-west-2.amazonaws.com/cell-census/2025-11-08/soma/census_data/homo_sapiens/obs/__schema/__enumerations/__0149b1dec606ee8d30c591c3fd6f86e2_0'
r=requests.get(url,headers={'Range':'bytes=0-182','x-amz-checksum-mode':'ENABLED'},timeout=25);r.raise_for_status();body=r.content;etag=r.headers.get('ETag','').strip('"');checksum=r.headers.get('x-amz-checksum-crc32');owncrc=base64.b64encode(struct.pack('!I',zlib.crc32(body)&0xffffffff)).decode();record=dict(status=r.status_code,bytes=len(body),etag=etag,md5=hashlib.md5(body).hexdigest(),server_crc32=checksum,computed_crc32=owncrc);assert len(body)==183 and record['md5']==etag
if checksum:assert checksum==owncrc
(O/'resources/CENSUS_FAILED_OBJECT_INDEPENDENT_CHECK.json').write_text(json.dumps(record,indent=2),encoding='utf-8');print('OBJECT_VERIFIED',json.dumps(record),flush=True)
config={'vfs.s3.region':'us-west-2','vfs.s3.no_sign_request':'true','sm.compute_concurrency_level':'2','sm.io_concurrency_level':'2','vfs.s3.max_parallel_ops':'1'};ctx=tiledbsoma.SOMATileDBContext(tiledb_config=config)
with cellxgene_census.open_soma(census_version='2025-11-08',context=ctx) as census:
 obs=census['census_data']['homo_sapiens'].obs;print('SCHEMA',obs.schema,flush=True)
 fl="dataset_id == '9f222629-9e39-47d0-b83f-e08d610c7479' and (cell_type == 'fibroblast' or cell_type == 'myofibroblast cell' or cell_type == 'pericyte' or cell_type == 'smooth muscle cell')"
 z=obs.read(value_filter=fl).concat().to_pandas();z.to_csv(O/'tables/census_stromal_obs.tsv.gz',sep='\t',index=False);print('STROMA',len(z),z.groupby(['cell_type','disease']).size().to_dict(),flush=True)
print('CENSUS_PROBE_COMPLETE',flush=True)
