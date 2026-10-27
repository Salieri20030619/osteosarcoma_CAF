from pathlib import Path
import re,json,xml.etree.ElementTree as ET
B=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=B/'phase64_response_followup_20261005'
p=O/'sources/GSE293065_Supplementary_Code.Rmd'
lines=p.read_text().splitlines();hits=[]
for i,l in enumerate(lines):
 if re.search(r'mask|patholog|annotat|boundary|necros|tumor|tumour|ROI|region',l,re.I):
  hits.append({'line':i+1,'text':'\n'.join(lines[max(0,i-2):min(len(lines),i+3)])})
(O/'resources/SPATIAL_CODE_AUDIT_EXCERPTS.json').write_text(json.dumps(hits,indent=2))
print('CODE_HITS',len(hits))
for h in hits[:35]:print(h)
print('CURRENT', (B/'CAF_NEXT_STAGE_CURRENT.json').read_text()[:6500])
root=ET.fromstring((O/'sources/GSE152048_paper.xml').read_text())
for x in root.iter('supplementary-material'):
 print('SUPPLEMENT',' '.join(x.itertext())[:1800])
