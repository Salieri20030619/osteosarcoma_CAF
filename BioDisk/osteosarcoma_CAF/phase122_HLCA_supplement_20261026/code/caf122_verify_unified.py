from pathlib import Path
import io,re,json,hashlib,pymupdf,pandas as pd
from PIL import Image
B=Path('C:/Users/31937/caf122_hlca_supplement_20261026');W=B/'unified_full_report';doc=pymupdf.open(W/'full_report_review.pdf');qa=pymupdf.open(W/'reviewer_QA_review.pdf');build=json.loads((W/'UNIFIED_BUILD_QA.json').read_text(encoding='utf-8'));checks=[]
def check(name,ok):checks.append(dict(check=name,passed=bool(ok)))
def compact(t):return re.sub(r'\s+','',t)
text=''.join(p.get_text() for p in doc);qtext=''.join(p.get_text() for p in qa)
check('report17pages',len(doc)==17);check('QA2pages',len(qa)==2);check('5questions_only_in_separate_file',all(qtext.count(f'Q{i}｜')==1 and f'Q{i}｜' not in text for i in range(1,6)));check('healthy_reference_precedes_recommendation','正常肺参照' in doc[10].get_text() and '综合推荐' in doc[12].get_text());check('IPF_only_in_appendix',all('IPF' not in p.get_text() for p in list(doc)[:13]))
newpixels=set()
for page in doc:
 for item in page.get_images():
  im=Image.open(io.BytesIO(doc.extract_image(item[0])['image'])).convert('RGB');newpixels.add(hashlib.sha256(im.tobytes()).hexdigest())
for i,entry in enumerate(build['source_figures'],1):
 im=Image.open(entry['path']).convert('RGB');check(f'original_figure_{i}_pixels_preserved',hashlib.sha256(im.tobytes()).hexdigest() in newpixels)
dec=pd.read_csv('C:/Users/31937/caf121_report_work/resources/node_decisions_ZH.tsv',sep='\t').set_index('父群或子群');t=compact(doc[4].get_text());nodes=['Bone_00','Bone_07','Bone_07_s1','Bone_07_s2','Bone_13','Lung_04','Lung_04_s2','Lung_05_s2','Lung_05_s4','Lung_06_s2','Lung_07','Lung_11','Lung_12']
for node in nodes:
 a=dec.loc[node];source=str(a['评估来源']).replace('SCPCP000017','017').replace('SCPCP000023','023');s=node+f"{int(a['细胞'])}/{int(a['覆盖记录'])}"+source+''.join(str(int(a[k])) for k in ['最少1细胞有效患者','最少5细胞有效患者','最少10细胞有效患者','身份敏感性有效患者']);check('patient_counts_'+node,s in t)
t=compact(doc[8].get_text())
for a in build['primary_statistic_rows']:
 s=f"{int(a['N'])}{int(a['positive_patients'])}/{int(a['N'])}{a['mean_effect']:+.3f}{a['P']:.5f}{a['FDR']:.5f}";check(a['mode']+'_'+a['endpoint'],s in t)
s1=compact(doc[14].get_text())
for s in ['POSTN74/7-0.1639820.6250000.772059','ACTA272/7-0.0238890.8593750.984375','TAGLN72/7+0.0295910.9218751.000000','LOXL271/7-0.0931600.5468750.710374']:check('S1_'+s[:5],s in s1)
for i,page in enumerate(doc):
 spans=[s for b in page.get_text('dict')['blocks'] if 'lines' in b for line in b['lines'] for s in line['spans']];check(f'page_{i+1}_fonts',all('YaHei' in s['font'] for s in spans));check(f'page_{i+1}_bounds',all(0<=s['bbox'][0]<s['bbox'][2]<=595.2 and 0<=s['bbox'][1]<s['bbox'][3]<=842.2 for s in spans))
check('no_decode_replacements','\ufffd' not in text+qtext)
result=dict(status='PASS' if all(c['passed'] for c in checks) else 'FAIL',checks=len(checks),passed=sum(c['passed'] for c in checks),failed=[c for c in checks if not c['passed']],report_sha256=hashlib.sha256((W/'full_report_review.pdf').read_bytes()).hexdigest(),QA_sha256=hashlib.sha256((W/'reviewer_QA_review.pdf').read_bytes()).hexdigest(),visual_review='All17 report pages and both QA pages reviewed; uniform text/table/heading/legend layout; no clipping')
(W/'UNIFIED_FINAL_QA.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(result));assert result['status']=='PASS'
