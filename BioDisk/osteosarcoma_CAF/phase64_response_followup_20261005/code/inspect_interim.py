from pathlib import Path
import pandas as pd,json
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase64_response_followup_20261005')
p=O/'code/analyze_response.py';s=p.read_text().replace("programs_evaluable=len(scores)","programs_evaluable=len(scores.columns)");p.write_text(s)
p=O/'sources/GPL14951_incomplete_attempt1.txt'
s=p.read_text(errors='replace');lines=s.splitlines()
print('INCOMPLETE_DIAGNOSTIC','bytes',p.stat().st_size,'rows',len(lines),'tail',lines[-1][:120])
print(pd.read_csv(O/'interim_GSE21257/tables/program_coverage.tsv',sep='\t')[['program','n_expected','n_measured','missing','evaluable']].to_string(index=False))
print(pd.read_csv(O/'interim_GSE21257/tables/program_response_statistics.tsv',sep='\t').query('model=="unadjusted_Welch"')[['program','effect','p','fdr_secondary13']].to_string(index=False))
