from pathlib import Path
import json, pandas as pd, numpy as np
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase43_full_transcriptome_discovery_20261001'
assert (O/'full_gene_statistics.exit_code').read_text().strip()=='0'
contexts=json.loads((O/'resources/EFFECT_CONTEXTS.json').read_text());rows=[];cached=None;table=None
symbols=['COL1A1','LOX','MMP14','RUNX2','ANGPT1','SMAD3','CCN2','LDHA','CA9','ACTA2','CXCL12','TGFB1']
features=pd.read_csv(P/'phase41_full_transcriptome_setup_20261001/tables/feature_presence_gene_index.tsv.gz',sep='\t').feature.to_numpy()
for ctx in contexts:
 if ctx['gene_result_file']!=cached:table=pd.read_csv(O/ctx['gene_result_file'],sep='\t');cached=ctx['gene_result_file']
 d=table[table.contrast.eq(ctx['contrast'])];selected=d[d.symbol.isin(symbols)]
 valid=d[d.P.notna()];extra=valid.sample(min(5,len(valid)),random_state=20261001)
 chosen=pd.concat([selected,extra]).drop_duplicates('feature')
 ef=np.load(O/ctx['effect_file']);lookup={features[j]:i for i,j in enumerate(ef['feature_indices'])}
 for row in chosen.itertuples(index=False):
  effects=ef['effects'][:,lookup[row.feature]]
  for donor,value in zip(ef['donors'],effects):rows.append({'context':ctx['context_id'],'feature':row.feature,'symbol':row.symbol,'donor':donor,'effect':value,'expected_N':row.N,'expected_P':row.P,'expected_mean':row.mean_delta,'expected_median':row.median_delta})
pd.DataFrame(rows).to_csv(O/'tables/independent_R_patient_effect_QA_input.tsv.gz',sep='\t',index=False)
print('INDEPENDENT_R_QA_INPUT_READY',len(rows))