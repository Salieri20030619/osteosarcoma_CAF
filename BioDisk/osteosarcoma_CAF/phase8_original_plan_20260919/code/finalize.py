from pathlib import Path
import json,hashlib,datetime,shutil
import numpy as np,pandas as pd
from statsmodels.stats.multitest import multipletests
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase8_original_plan_20260919'
for k in ['extraction','validation','bulk_analysis','lung_reference','delivery']:
 assert (O/'logs'/(k+'.exit_code')).read_text().strip()=='0',k
inv=pd.read_csv(O/'tables/analysis_sample_inventory.tsv',sep='\t')
assert len(inv)==33 and inv['sample'].is_unique and inv.cells.sum()==226699
assert inv.loc[inv.dataset=='GSE152048','strict_CAF'].sum()==2047
assert inv.loc[inv.source=='new_GEO','strict_CAF'].sum()==2065
ef=pd.read_csv(O/'tables/full_customer_gene_effects.tsv',sep='\t');assert len(ef)==62 and ef[['state','gene']].drop_duplicates().shape[0]==62
bulk=pd.read_csv(O/'tables/GSE21257_primary_associations.tsv',sep='\t');assert len(bulk)==8 and (bulk.n==53).all() and (bulk.events==34).all()
assert np.allclose(multipletests(bulk.P_value,method='fdr_bh')[1],bulk.FDR_BH)
assert ((bulk.CI_low<=bulk.OR_per_SD)&(bulk.OR_per_SD<=bulk.CI_high)).all()
proj=pd.read_csv(O/'tables/external_projection_counts.tsv',sep='\t')
assert proj.groupby('version').cells.sum().eq(2065).all()
assert proj.loc[proj.projected_state!='unassigned','cells'].sum()==0
hashes=pd.read_csv(O/'resources/input_state_hashes.tsv',sep='\t')
for r in hashes.to_dict('records'):
 assert hashlib.sha256(Path(r['path']).read_bytes()).hexdigest()==r['sha256'],r['path']
assert len(list((O/'delivery').glob('*.pdf')))==1
timestamp=datetime.datetime.now(datetime.timezone.utc).isoformat()
checks={'completed_utc':timestamp,'source_hashes':'PASS','sample_cell_conservation':'PASS','full_customer_panels':'16 osteogenic plus 15 calcium genes, 62 state-gene comparisons','bulk_clinical_alignment':'53 patients, 34 events','multiple_testing':'8 prespecified endpoints, BH recomputation PASS','failed_state_transfer_preserved':True,'figure_QA':'4/4 PASS','PDF_QA':'10/10 pages visually inspected; PASS','unchanged_original_labels':True}
(O/'FINAL_ACCEPTANCE.json').write_text(json.dumps(checks,indent=2))
(O/'CORE_COMPLETE').write_text('ORIGINAL_PLAN_EXECUTABLE_CORE_COMPLETE_WITH_DOCUMENTED_DEPENDENCIES\n')
# Subsequent reference-paper assessment does not alter completed primary analyses.
memo="""# 原计划完成后的参考论文方法评估

原计划主线已完成并冻结。参考论文为 Zeng 等 Molecular Cancer 25:14，DOI 10.1186/s12943-025-02516-2。

可借鉴：按CD36、MYH11、C7、FNDC1相关多基因程序提供独立的参考注释，并保留未匹配或混合状态。不可预设四类必然存在。

当前限制：作者原始单细胞参考对象未取得；其MYH11群与本项目独立排除的Mural/SMC存在范围差异。当前严格集合不能用来否定这类细胞存在。若以后评估，应并列使用保留的Mural/SMC候选，原严格CAF标签保持不变。

本项目自身H0/H2状态的供者留出检验未通过，说明离散标签外推需要额外验证。不能通过更换到论文四个名称绕过该问题。CD36单基因稀疏表达和校正计数阈值敏感性也不能单独确立新亚群。

结论：论文方法可作为后续独立表达程序比较，当前不适合直接替换原计划分群或宣称复现四类CAF及抗PD-1疗效关系。本轮仅完成方法可用性评估，没有启动论文四分群计算。

来源：https://link.springer.com/article/10.1186/s12943-025-02516-2
"""
(O/'delivery/reference_method_assessment.md').write_text(memo,encoding='utf-8')
state={'status':'ORIGINAL_PLAN_CORE_COMPLETE_WITH_DOCUMENTED_DEPENDENCIES','completed_utc':timestamp,'main_report':'delivery/original_plan_report.md','pdf_report':str(next((O/'delivery').glob('*.pdf')).relative_to(O)),'discovery_CAF':2047,'new_GEO_CAF':2065,'bulk_patients':53,'bulk_FDR_significant_programs':['Antigen_presentation'],'state_transfer_CV_pass':False,'myCAF_subclusters_independently_supported':0,'surface_panel':'candidate enrichment design only; protein validation pending','conditional_items':['A/B1 and contemporaneous bone-lung pairs require clinical metadata','late ScPCA expansion excluded from first version per original plan','instrument channels and digestion conditions needed for experimental configuration'],'reference_paper_assessment':'complete; separate exploratory program comparison recommended, four-group computation not started'}
(O/'STATUS.json').write_text(json.dumps(state,indent=2))
(O/'README.md').write_text("""# Original CAF plan delivery

Read delivery/original_plan_report.md or the PDF in delivery/.
Tables contain complete effects, negative validations, cohort metadata, RNA combination risks, protein and antibody evidence.
validation/ holds diagnostic projections; none were accepted because donor CV failed the frozen gate.
resources/ contains frozen protocols and provenance.
The partial GPL family gzip remains only on sgs10; the COMPLETE 48,701-row platform table and full series matrix were accepted and used.
Raw full-expression objects and earlier full-cell atlases remain at:
../phase7_review_states_20260917/figures/
../geo_preprocess_20260918/accepted/
../geo_identity_20260919/
These are frozen prior outputs, not replaced by this delivery.
No new clinical patient pooling, experimental validation, or reference-paper four-group clustering is claimed.
Reproduction order: prepare.py; accept_bulk_sources.py; extract_evidence.py; analyze_validation.py; lung_reference.py; bulk_analysis.py; build_delivery.py; local render_report_pdf.py; finalize.py.
Scripts use explicit version roots. Use a new output version for re-analysis.
""",encoding='utf-8')
main=P/'PROJECT_STATE.json';shutil.copy2(main,P/('PROJECT_STATE.before_original_plan_'+datetime.datetime.now().strftime('%Y%m%d_%H%M%S')+'.json'))
project=json.loads(main.read_text());project.update(updated_at_utc=timestamp,current_stage=O.name,dynamic_status_file=O.name+'/STATUS.json',status=state['status'],current_status=state['status'],next_action='Review original-plan delivery; conditional clinical/experimental enhancements remain listed. Reference-paper four-group analysis not started.',latest_original_plan_report=str(O/'delivery/original_plan_report.md'))
project['phase8_original_plan_20260919']=state;project['completed_stages']=list(dict.fromkeys(project.get('completed_stages',[])+[O.name]))
main.write_text(json.dumps(project,indent=2,ensure_ascii=False),encoding='utf-8')
files=[]
for dirname in ['code','tables','figures','resources','validation','delivery','sources','logs']:
 files += [f for f in (O/dirname).rglob('*') if f.is_file() and '__pycache__' not in str(f)]
files += [O/n for n in ['STATUS.json','README.md','FINAL_ACCEPTANCE.json','PROTOCOL_FROZEN.json','analysis_summary.json','CORE_COMPLETE']]
files += [O/'bulk'/n for n in ['data_acceptance.json','analysis_summary.json','GSE21257_series_matrix.txt.gz','GPL10295_platform_table.tsv','GSE21257_gene_median_expression.tsv.gz','DOWNLOAD_COMPLETE']]
rows=[]
for f in sorted(files):
 with open(f,'rb') as h:digest=hashlib.file_digest(h,'sha256').hexdigest()
 rows.append({'path':str(f.relative_to(O)),'bytes':f.stat().st_size,'sha256':digest})
pd.DataFrame(rows).to_csv(O/'SHA256SUMS.tsv',sep='\t',index=False)
print(json.dumps(state),flush=True)
