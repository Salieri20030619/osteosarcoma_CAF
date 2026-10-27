from probe_sources import ROOT,fetch
import concurrent.futures,json,xml.etree.ElementTree as ET
jobs=[]
for n in (1,2,3):
 ext={1:'docx',2:'xls',3:'xlsx'}[n]
 jobs.append((f'shen_sup{n}',f'https://static-content.springer-cdn.com/esm/art%3A10.1007%2Fs13402-026-01162-0/MediaObjects/13402_2026_1162_MOESM{n}_ESM.{ext}','proxy'))
jobs += [('hra_supplement','https://image.oaes.cc/published/article/82f5bd4c4ceff0b525e5e945aafb3d02/cdr9009-SupplementaryMaterials.pdf','proxy'),('ncc_methods','https://static-content.springer-cdn.com/esm/art%3A10.1186%2Fs13045-024-01598-7/MediaObjects/13045_2024_1598_MOESM4_ESM.docx','proxy'),('hra_omix_search','https://ngdc.cncb.ac.cn/omix/search?query=HRA013407','direct'),('hra_ngdc_root','https://ngdc.cncb.ac.cn/','direct'),('hra_human_download_root','https://download.cncb.ac.cn/','direct'),('hra_geo_related','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE317512&targ=self&form=text&view=full','proxy')]
with concurrent.futures.ThreadPoolExecutor(max_workers=5) as ex:rows=list(ex.map(fetch,jobs))
(ROOT/'supplement_routes.json').write_text(json.dumps(rows,indent=2),encoding='utf8')
root=ET.parse(ROOT/'evidence/shen_paper_proxy.txt').getroot()
texts=[' '.join(e.itertext()) for e in root.iter() if e.tag in ('p','fn') and any(w in ' '.join(e.itertext()).lower() for w in ['gene expression omnibus','data avail','deposited','sequencing','patients'])]
(ROOT/'evidence/shen_relevant_text.txt').write_text('\n\n'.join(texts),encoding='utf8')
print('\nSHEN RELEVANT TEXT\n'+'\n'.join(texts)[:8000])
