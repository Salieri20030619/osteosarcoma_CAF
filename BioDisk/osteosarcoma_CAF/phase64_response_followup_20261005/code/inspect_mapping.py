from pathlib import Path
import json
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase64_response_followup_20261005')
for g in ['gpl14951_data','gpl10295_data']:
 s=(O/'sources'/f'probe_{g}.txt').read_text();i=s.find('!platform_table_begin');print(g,s[i:i+2500])
s=(O/'sources/GSE293065_Supplementary_Code.Rmd').read_text().splitlines()
print('CANCER_SELECTION','\n'.join(s[180:212]))
p=O.parent/'phase63_matrix_resistance_20261005/resources/PROGRAMS.json'
print('PROGRAMS',[(k,len(v)) for k,v in json.loads(p.read_text())['sets'].items()])
