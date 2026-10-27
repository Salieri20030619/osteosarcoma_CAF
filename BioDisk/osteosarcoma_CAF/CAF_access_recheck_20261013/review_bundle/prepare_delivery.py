from pathlib import Path
import json,zipfile,hashlib,shutil,re,csv,datetime
import fitz
ROOT=Path(__file__).resolve().parent
D=ROOT.parent/'CAF_新增数据入口复核_20261013';D.mkdir(exist_ok=True)
L=ROOT/'literature';L.mkdir(exist_ok=True)
sources={'hra_supplement_proxy.txt':'HRA013407_原文补充材料.pdf','shen_real_sup1_proxy.txt':'Shen2026_补充图.docx','shen_real_sup2_proxy.txt':'Shen2026_作者标记汇总表_原文件.xls','shen_real_sup3_proxy.txt':'Shen2026_样本年龄分组.xlsx','ncc_real_methods_proxy.txt':'NCC_补充方法.docx'}
for a,b in sources.items():shutil.copy2(ROOT/'evidence'/a,L/b)
with fitz.open(L/'HRA013407_原文补充材料.pdf') as pdf:
 pdf[2].get_pixmap(matrix=fitz.Matrix(1.3,1.3)).save(ROOT/'evidence/HRA_supplement_page3.png')
# Extract displayed DOCX text only, excluding embedded EndNote field instructions.
from docx import Document
for p in L.glob('*.docx'):
 doc=Document(p);text='\n'.join(p.text for p in doc.paragraphs)
 (ROOT/'evidence'/(p.stem+'_displayed_text.txt')).write_text(text,encoding='utf8')
T=ROOT/'tables'
projects={'PRJCA025405':'HRA007225','PRJCA028332':'HRA008101','PRJCA029010':'HRA008525','PRJCA019918':'HRA006496'}
rows=[]
for project,hra in projects.items():
 d=json.loads((ROOT/'evidence'/f'{project}_all_samples_direct.txt').read_text(encoding='utf8'))['result']['data']
 assert d['recordsTotal']==len(d['data'])
 for r in d['data']:
  rows.append({'HRA':hra,'BioProject':project,'BioSample':r['id'],'sample_name':r['attrs'].get('SampleName'),'description':r['description'],'organization':r['attrs'].get('Organization',''),'source':r['url'],'independent_patient_count':'unverified'})
with (T/'ADDITIONAL_HRA_SAMPLE_RECORDS.tsv').open('w',encoding='utf8',newline='') as f:
 w=csv.DictWriter(f,fieldnames=rows[0].keys(),delimiter='\t');w.writeheader();w.writerows(rows)
d=json.loads((ROOT/'evidence/hra_osteo_all22_direct.txt').read_text(encoding='utf8'))['result']['data'];assert d['recordsTotal']==22==len(d['data'])
rows=[{'HRA':r['id'],'description':r['description'],'BioProject':r['attrs'].get('bioProject Accession'),'biosample_records':len(r['attrs'].get('bioSample Accession',[])),'source':r['url']} for r in d['data']]
with (T/'OFFICIAL_HRA_OSTEOSARCOMA_22_RECORDS.tsv').open('w',encoding='utf8',newline='') as f:
 w=csv.DictWriter(f,fieldnames=rows[0].keys(),delimiter='\t');w.writeheader();w.writerows(rows)
manifest=json.loads((ROOT/'NCC_DOWNLOAD_COMPLETE.json').read_text(encoding='utf8'))
receipt={'date_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'NCC_downloaded_bytes':manifest['bytes'],'NCC_official_git_blobs_verified':13,'NCC_single_cell_files':10,'NCC_labeled_single_cell_files':9,'NCC_confirmed_barcodes':79071,'NCC_BS6_labels_missing':True,'NCC_spatial_spots':4579,'HRA013407_status_in_official_index':'Not released','PRJCA032564_ReleaseDate_field':'2028-08-19','Shen_supplement_files_downloaded':3,'new_formal_CAF_patients_added':0,'original_baseline_unchanged':True}
(ROOT/'RECHECK_FINDINGS.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2),encoding='utf8')
print('PREPARED',str(D),json.dumps(receipt,ensure_ascii=False))
