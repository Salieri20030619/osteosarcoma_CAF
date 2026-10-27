from probe_sources import ROOT,fetch
import concurrent.futures,json
jobs=[('hra_osteo_all22','https://ngdc.cncb.ac.cn/search/api/specific?db=hra&q=osteosarcoma&sort=desc&length=100&start=0','direct'),('hra_project032564','https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA032564','direct'),('hra_project_search','https://ngdc.cncb.ac.cn/search/api/specific?db=bioproject&q=PRJCA032564','direct'),('hra_sameproject_omix','https://ngdc.cncb.ac.cn/search/api/specific?db=omix&q=PRJCA032564','direct')]
for p in ['PRJCA057037','PRJCA024817','PRJCA028332','PRJCA025405']:
 jobs.append((p,'https://ngdc.cncb.ac.cn/bioproject/browse/'+p,'direct'))
with concurrent.futures.ThreadPoolExecutor(max_workers=5) as ex:rows=list(ex.map(fetch,jobs))
(ROOT/'hra_crosslink_routes.json').write_text(json.dumps(rows,indent=2),encoding='utf8')
