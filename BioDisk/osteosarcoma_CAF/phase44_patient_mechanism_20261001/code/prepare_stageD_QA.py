from pathlib import Path
import json,numpy as np,pandas as pd
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');C=P/'phase43_full_transcriptome_discovery_20261001';O=P/'phase44_patient_mechanism_20261001'
assert (O/'robustness_tf.exit_code').read_text().strip()=='0' and (O/'clinical.exit_code').read_text().strip()=='0'
contexts=json.loads((C/'resources/EFFECT_CONTEXTS.json').read_text());cross=pd.read_csv(C/'tables/HGNC_full_gene_symbol_crosswalk.tsv.gz',sep='\t',keep_default_na=False);features=cross.feature.to_numpy();symbols=cross.canonical_symbol.to_numpy();quality=cross.mapping_method.eq('HGNC_Ensembl_approved').map({True:0,False:1}).to_numpy();net=pd.read_csv(O/'tables/signed_TF_target_network_used.tsv.gz',sep='\t');selected=net[net.regulator.isin(['SMAD3','TEAD1','JUN','STAT3','RUNX2'])];rows=[]
for ctx in contexts:
 ef=np.load(C/ctx['effect_file']);gi=ef['feature_indices'];pref={}
 for k,j in enumerate(gi):
  key=(quality[j],features[j]);sy=symbols[j]
  if sy and (sy not in pref or key<pref[sy][0]):pref[sy]=(key,k)
 for regulator,g in selected.groupby('regulator'):
  for target,weight in g[['target_symbol','weight']].itertuples(index=False,name=None):
   if target not in pref:continue
   values=ef['effects'][:,pref[target][1]]
   for donor,value in zip(ef['donors'],values):rows.append({'context':ctx['context_id'],'regulator':regulator,'donor':donor,'target':target,'weight':weight,'gene_effect':value})
pd.DataFrame(rows).to_csv(O/'tables/independent_TF_target_reconstruction_input.tsv.gz',sep='\t',index=False)
print('STAGE_D_QA_INPUT_READY',len(rows))