from probe_sources import ROOT,fetch
import concurrent.futures,json
jobs=[
('hra_commonjs','https://ngdc.cncb.ac.cn/search/static/js/my/common.js','direct'),
('hra_omix_list','https://ngdc.cncb.ac.cn/omix/releaseList','direct'),
('hra_direct_directory','https://download.cncb.ac.cn/gsa-human/HRA013407/','direct'),
('shen_supp_epmc_direct','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12816095/supplementaryFiles','direct'),
('ncc_supp_epmc_direct','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11337771/supplementaryFiles','direct'),
('shen_springer_page','https://link.springer.com/article/10.1007/s13402-026-01162-0','proxy'),
('hra_geo_direct','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE317512&targ=self&form=text&view=full','direct'),
('shen_geo_esearch','https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=osteosarcoma%20AND%20%28Shen%20OR%20aging%20OR%20Fujian%29&retmode=json&retmax=100','direct')]
with concurrent.futures.ThreadPoolExecutor(max_workers=5) as ex:rows=list(ex.map(fetch,jobs))
(ROOT/'next_routes.json').write_text(json.dumps(rows,indent=2),encoding='utf8')
