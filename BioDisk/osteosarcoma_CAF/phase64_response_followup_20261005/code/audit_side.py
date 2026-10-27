from pathlib import Path
import requests,json,re,concurrent.futures
B=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=B/'phase64_response_followup_20261005'
print('statefiles',[x.name for x in B.iterdir() if any(k in x.name.upper() for k in ['CURRENT','STATUS','AGENTS'])])
for f in B.glob('**/AGENTS.md'):
 print('AGENTS',f, f.read_text()[:10000])
 break
for g in ['GPL14951','GPL10295']:
 p=O/'sources'/f'{g}_platform_only.soft'
 if p.exists():
  t=p.read_text();i=t.index('!platform_table_begin'); print(g,len(t),t[i:i+1800])
def get(url,name):
 try:
  r=requests.get(url,timeout=(15,30));r.raise_for_status();(O/'sources'/name).write_text(r.text);print(name,len(r.content),r.text[:100])
 except Exception as e:print(name,str(e))
with concurrent.futures.ThreadPoolExecutor(2) as ex:
 list(ex.map(lambda v:get(*v),[
 ('https://zenodo.org/api/records/19369637/files/Supplementary_Code.Rmd/content','GSE293065_Supplementary_Code.Rmd'),
 ('https://codeocean.com/capsule/9535428/tree/v1','BS3_codeocean_landing.html')]))
