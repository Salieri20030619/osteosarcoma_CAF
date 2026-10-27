from pathlib import Path
import requests, json, concurrent.futures, datetime, hashlib
ROOT=Path(__file__).resolve().parent
ROOT.mkdir(exist_ok=True)
OUT=ROOT/'evidence'; OUT.mkdir(exist_ok=True)
jobs={
 'ncc_repo_tree':'https://api.github.com/repos/zhengxj1/A-Single-Cell-and-Spatially-Resolved-Atlas-of-Human-Osteosarcomas/git/trees/main?recursive=1',
 'ncc_releases':'https://api.github.com/repos/zhengxj1/A-Single-Cell-and-Spatially-Resolved-Atlas-of-Human-Osteosarcomas/releases',
 'ncc_readme':'https://raw.githubusercontent.com/zhengxj1/A-Single-Cell-and-Spatially-Resolved-Atlas-of-Human-Osteosarcomas/main/README.md',
 'ncc_capsule':'https://codeocean.com/capsule/9535428/tree/v1',
 'hra_current':'https://ngdc.cncb.ac.cn/gsa-human/browse/HRA013407',
 'hra_ngdc_old':'https://ngdc.cncb.ac.cn/gsa-human/browse/HRA013407/',
 'hra_bigd':'https://bigd.big.ac.cn/gsa-human/browse/HRA013407',
 'hra_paper':'https://www.oaepublish.com/articles/cdr.2026.09',
 'shen_paper':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12816095/fullTextXML',
 'hra_pmc':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13174192/fullTextXML',
 'ncc_pmc':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11337771/fullTextXML',
}
def fetch(item):
 name,url,mode=item
 s=requests.Session(); s.trust_env=False
 if mode=='proxy': s.proxies={'http':'http://127.0.0.1:7897','https':'http://127.0.0.1:7897'}
 row={'name':name,'url':url,'route':mode,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
 try:
  r=s.get(url,timeout=(12,35)); p=OUT/(name+'_'+mode+'.txt'); p.write_bytes(r.content)
  row.update(status=r.status_code,bytes=len(r.content),final_url=r.url,file=p.name,sha256=hashlib.sha256(r.content).hexdigest())
 except Exception as e: row['error']=str(e)
 print(json.dumps(row,ensure_ascii=False),flush=True)
 return row
if __name__=='__main__':
 items=[(k,v,m) for k,v in jobs.items() for m in (['direct','proxy'] if k.startswith(('hra_current','hra_bigd','ncc_capsule')) else ['proxy'])]
 with concurrent.futures.ThreadPoolExecutor(max_workers=6) as ex: rows=list(ex.map(fetch,items))
 (ROOT/'route_probes.json').write_text(json.dumps(rows,ensure_ascii=False,indent=2),encoding='utf-8')
