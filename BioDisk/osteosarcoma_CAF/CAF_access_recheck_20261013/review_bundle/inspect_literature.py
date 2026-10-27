from pathlib import Path
import xml.etree.ElementTree as ET,re,json
ROOT=Path(__file__).resolve().parent
for name in ['shen_paper','hra_pmc','ncc_pmc']:
 r=ET.parse(ROOT/'evidence'/f'{name}_proxy.txt').getroot()
 print('\nARTICLE',name)
 for e in r.iter():
  if e.tag in ['supplementary-material','ext-link','media']:
   print(e.tag,e.attrib,' '.join(e.itertext())[:350])
  if e.tag=='sec' and e.find('title') is not None and re.search('availability|sequencing|cohort', ''.join(e.find('title').itertext()),re.I):
   print('SECTION',' '.join(e.itertext())[:2000])
 print('ACCESSIONS',sorted(set(re.findall(r'(?:GSE|GSM|HRA|CRA|PRJNA|PRJCA|SRP|SRR)\d+',ET.tostring(r,encoding='unicode')))))
