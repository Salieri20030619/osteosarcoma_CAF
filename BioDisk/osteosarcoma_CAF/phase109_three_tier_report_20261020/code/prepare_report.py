import pathlib,json,hashlib,re,zipfile,xml.etree.ElementTree as ET
import numpy as np,pandas as pd
R=pathlib.Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');P=R/'phase109_three_tier_report_20261020'
for d in ['code','resources','tables','figures','report','render']:(P/d).mkdir(parents=True,exist_ok=True)
lit=R/'phase96_all14_literature_20261016/revision_volcano'
content=json.loads((lit/'content.json').read_text())
six=['PNKD','TOR1AIP1','GRID2','TAF10','BEST1','RGS2-AS1']
bio=[g for g in content['genes'] if g['gene'] in six]
refs=[r for r in content['references'] if r['id'] in [7,8,9,10,11,12,15,19,20,27]]
with zipfile.ZipFile(lit/'CAF_14gene_literature_20261016.docx') as z:txt=''.join(ET.fromstring(z.read('word/document.xml')).itertext())
assert all(g in txt for g in six)
(P/'resources/previous_six_gene_biology.json').write_text(json.dumps(dict(genes=bio,references=refs,source=str(lit/'CAF_14gene_literature_20261016.docx')),ensure_ascii=False,indent=2),encoding='utf-8')
old=pd.read_csv(R/'phase81_allCAF_bone_lung_20261012/tables/all_genes_statistics.tsv.gz',sep='\t')
old=old[old.layer.eq('expanded')&old.sensitivity.eq('all_cells')&old.design.eq('independent')&old.dataset.isin(['SCPCP000017','SCPCP000023'])]
g=pd.read_csv(R/'phase108_each_source_significant_20261020/tables/all_genes_each_source_gate.tsv.gz',sep='\t')
registry=json.loads((lit/'project_registry.json').read_text())
order=[v['gene'] for v in registry];mapping={v['feature_ids']:v['gene'] for v in registry}
st1=g.P_17.lt(.05)&g.P_23.lt(.05)&g.source17_effect.mul(g.source23_effect).gt(0)
st2=st1&g.baseline_selected
st3=st2&g.remain_after_both_corrections
assert (int(st1.sum()),int(st2.sum()),int(st3.sum()))==(14,8,6)
assert set(g.loc[st3,'gene'])==set(six)
c=g[st1].set_index('gene').loc[order].reset_index();c['tier1']=True;c['tier2']=c.baseline_selected;c['tier3']=c.remain_after_both_corrections
c.to_csv(P/'tables/tier_candidate_registry.tsv',sep='\t',index=False,encoding='utf-8',na_rep='NA')
g[st3].set_index('gene').loc[six].reset_index().to_csv(P/'tables/final_six_complete_statistics.tsv',sep='\t',index=False,encoding='utf-8',na_rep='NA')
specs=[]
for src,n in [('SCPCP000017','tier1_source17'),('SCPCP000023','tier1_source23')]:
 z=old[old.dataset.eq(src)].copy().rename(columns={'mean_lung_minus_bone':'effect'})
 z['gene']=[mapping.get(k,s) for k,s in zip(z.feature,z.symbol)]
 z[['feature','gene','effect','P','FDR']].to_csv(P/'tables'/(n+'_volcano_source.tsv.gz'),sep='\t',index=False,encoding='utf-8',na_rep='NA')
 specs.append(dict(key=n,stage=1,source=src,valid=int(z.P.notna().sum()),nominal=int(z.P.lt(.05).sum()),FDR=int(z.FDR.lt(.05).sum())))
for key,stage,file,eff in [
 ('tier2_pooled',2,R/'phase106_pooled_test_original_directions_20261020/tables/patient_permutation_allgenes_with_directions.tsv.gz','mean_lung_minus_bone'),
 ('tier2_NCI',2,R/'phase106_pooled_test_original_directions_20261020/tables/third_NCI_vs_Bone_allgenes.tsv.gz','mean_lung_minus_bone'),
 ('tier3_pooled',3,R/'phase107_third_NCI_normal_background_20261020/tables/pooled_normal_allgenes.tsv.gz','effect'),
 ('tier3_NCI',3,R/'phase107_third_NCI_normal_background_20261020/tables/NCI_normal_allgenes.tsv.gz','effect')]:
 z=pd.read_csv(file,sep='\t').rename(columns={eff:'effect'}) if eff!='effect' else pd.read_csv(file,sep='\t')
 z[['feature','gene','effect','P','FDR']].to_csv(P/'tables'/(key+'_volcano_source.tsv.gz'),sep='\t',index=False,encoding='utf-8',na_rep='NA')
 specs.append(dict(key=key,stage=stage,source=str(file),valid=int(z.P.notna().sum()),nominal=int(z.P.lt(.05).sum()),FDR=int(z.FDR.lt(.05).sum())))
(P/'resources/FIGURE_SPECS.json').write_text(json.dumps(specs,indent=2),encoding='utf-8')
inputs=[lit/'content.json',lit/'CAF_14gene_literature_20261016.docx',R/'phase108_each_source_significant_20261020/tables/all_genes_each_source_gate.tsv.gz']
(P/'resources/INPUT_MANIFEST.json').write_text(json.dumps([dict(path=str(f),sha256=hashlib.sha256(f.read_bytes()).hexdigest()) for f in inputs],indent=2),encoding='utf-8')
(P/'resources/REPORT_SCOPE.json').write_text(json.dumps(dict(format='PDF',tiers=[14,8,6],volcano_count=6,normal_min_cells=20,final_six=six,full_gene_FDR_preserved=True,no_new_significance_tests=True,source_adjusted_model_not_in_this_report=True,biology_source='20261016 revision_volcano content, six sections retained with evidence levels',web_rechecked_primary_articles=['10.1155/humu/8430177','10.1038/s41556-022-01042-3','10.1007/s11010-025-05276-1','PMC8085875','10.1002/ctm2.1153','10.1074/jbc.M703758200','10.1038/s41413-022-00197-x','10.1002/advs.202306472']),ensure_ascii=False,indent=2),encoding='utf-8')
print('PREPARED',specs)
