import fsspec,h5py,json,time
from pathlib import Path
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase122_HLCA_supplement_20261026');url='https://datasets.cellxgene.cziscience.com/5f863718-02a6-46cb-83e7-52be988b915b.h5ad'
t=time.time()
with fsspec.open(url,'rb',block_size=1048576,cache_type='blockcache') as f:
 with h5py.File(f,'r') as h:
  print('ROOT',list(h),flush=True)
  for key in ['X','raw','layers','obs','var']:
   g=h.get(key)
   if g is None:continue
   print(key,dict(g.attrs),list(g)[:60] if isinstance(g,h5py.Group) else g.shape,flush=True)
  print('TIME',time.time()-t,flush=True)
