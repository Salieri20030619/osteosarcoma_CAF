from pathlib import Path
import pandas as pd,json
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase62_F1F5_internal_20261005');jobs=json.loads((O/'resources/GSEA_JOBS.json').read_text());names=json.loads((O/'resources/WORKING_NAMES.json').read_text());focus={
'F1_originalC10':['NABA_CORE_MATRISOME','HALLMARK_INFLAMMATORY_RESPONSE','HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION'],
'F1_originalC4':['HALLMARK_TNFA_SIGNALING_VIA_NFKB','HALLMARK_G2M_CHECKPOINT'],
'F1_originalC5':['HALLMARK_BILE_ACID_METABOLISM','HALLMARK_FATTY_ACID_METABOLISM'],
'F2_C0':['NABA_PROTEOGLYCANS'],'F2_C1':['REACTOME_SMOOTH_MUSCLE_CONTRACTION','HALLMARK_MYOGENESIS'],
'F2_C2':['NABA_SECRETED_FACTORS'],'F2_C3':['NABA_COLLAGENS','REACTOME_COLLAGEN_FORMATION'],
'F2_C4':['HALLMARK_E2F_TARGETS','HALLMARK_G2M_CHECKPOINT'],'F2_C5':['HALLMARK_OXIDATIVE_PHOSPHORYLATION','HALLMARK_MYC_TARGETS_V1'],
'F3_C0':['HALLMARK_IL6_JAK_STAT3_SIGNALING'],'F3_C1':['NABA_CORE_MATRISOME','NABA_ECM_GLYCOPROTEINS'],'F3_C2':['HALLMARK_MYC_TARGETS_V1','REACTOME_TRANSLATION'],
'F4_C1':['NABA_CORE_MATRISOME','REACTOME_COLLAGEN_FORMATION'],'F4_C2':['HALLMARK_UNFOLDED_PROTEIN_RESPONSE'],
'F4_C5':['HALLMARK_INFLAMMATORY_RESPONSE','NABA_BASEMENT_MEMBRANES'],'F4_C7':['HALLMARK_E2F_TARGETS','HALLMARK_G2M_CHECKPOINT'],
'F4_C8':['NABA_PROTEOGLYCANS','NABA_CORE_MATRISOME'],'F4_C9':['NABA_CORE_MATRISOME'],
'F4_C10':['REACTOME_RUNX2_REGULATES_OSTEOBLAST_DIFFERENTIATION','NABA_ECM_REGULATORS'],
'F4_C13':['HALLMARK_INTERFERON_ALPHA_RESPONSE'],'F4_C14':['HALLMARK_MYC_TARGETS_V1','HALLMARK_OXIDATIVE_PHOSPHORYLATION'],
'F4_C16':['HALLMARK_E2F_TARGETS','HALLMARK_G2M_CHECKPOINT'],
'F5_originalC1':['NABA_CORE_MATRISOME'],'F5_originalC12':['HALLMARK_INFLAMMATORY_RESPONSE','HALLMARK_MYC_TARGETS_V1'],
'F5_originalC6':['REACTOME_RUNX2_REGULATES_OSTEOBLAST_DIFFERENTIATION']}
frames=[];done=0
for j in jobs:
 if not Path(j['output']).exists():continue
 done+=1;g=pd.read_csv(j['output'],sep='\t');v=g[g.Term.isin(focus.get(j['contrast'],[]))];frames.append(v)
out=pd.concat(frames);out.to_csv(O/'tables/functional_name_evidence.tsv',sep='\t',index=False);print('JOBS',done,len(jobs))
for st,g in out.groupby('contrast'):
 print(st)
 for r in g.itertuples():print(r.scope,r.Term,round(r.NES,3),'P',round(r.GSEA_model_P,5),'BH',round(r.model_BH_FDR,5))
for j in jobs:
 if j['scope']!='all' or j['library']!='hallmark' or j['contrast'] not in ['F1_originalC4','F4_C1','F4_C11','F4_C12','F4_C15','F4_C16','F5_originalC1','F5_originalC23','F5_originalC6']:continue
 if not Path(j['output']).exists():continue
 g=pd.read_csv(j['output'],sep='\t');g=g[g.NES.gt(0)&g.GSEA_model_P.lt(.05)].sort_values('NES',ascending=False).head(5);print('ADDITIONAL',j['contrast'],'; '.join(f'{r.Term}:{r.NES:.2f},P={r.GSEA_model_P:.4g},BH={r.model_BH_FDR:.4g}' for r in g.itertuples()))
