from pathlib import Path
import json,csv,datetime,subprocess,hashlib,shutil
import pandas as pd,anndata as ad,numpy as np
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase31_customer_completion_20260926';B=P/'phase26_client_reanalysis_20260922'
for d in ['code','resources','tables','logs','delivery','monitoring']:(O/d).mkdir(parents=True,exist_ok=True)
now=datetime.datetime.now().astimezone().isoformat()
ev=pd.read_csv(P/'phase30_two_layer_classification_20260926/tables/expanded_classified_cells.tsv.gz',sep='\t',low_memory=False).set_index('cell_id')
need=[]
for sm,g in ev.groupby('sample'):
 a=ad.read_h5ad(B/'inputs'/sm/'full_called_counts.h5ad',backed='r');x=a.obs.loc[g.index]
 sing=x.get('scDblFinder_class',pd.Series('',index=x.index)).astype(str).str.lower().eq('singlet')
 scored=pd.to_numeric(x.get('doublet_score',pd.Series(np.nan,index=x.index)),errors='coerce').notna()
 if (~(sing|scored)).any():need.append(sm)
 a.file.close()
manifest=pd.read_csv(B/'tables/analysis_source_manifest.tsv',sep='\t');assert len(manifest)==107
# Interleave studies; every input remains in the queue.
groups=[list(g['sample']) for _,g in manifest.groupby('dataset')];order=[]
while any(groups):
 for g in groups:
  if g:order.append(g.pop(0))
jobs=[{'method':'decontx','sample':sm} for sm in order]+[{'method':'doublet','sample':sm} for sm in need]
(O/'METHOD_JOBS.json').write_text(json.dumps(jobs,indent=2))
registry={
'created_at':now,'user_authorization':'Customize next plan and start; continue without asking while user is away; do not add/remove customer thresholds.',
'customer_explicit':{'human_nonPDX':True,'lymph_node_not_B2':True,'CNV_crosscheck_methods_at_least':2,'flow_positive_surface_markers_at_least':2,'stress_gene_set_size':512,'myCAF_genes':['ACTA2','TAGLN','MYH11','MYL9','TPM2','CNN1'],'iCAF_genes':['IL6','CXCL12','PDGFRA','DPT','HAS1'],'apCAF_genes':['CD74','HLA-DRA','HLA-DRB1'],'osteogenic_panel_size':16,'mineral_panel_size':15},
'user_scope':{'supported_primary_cells':1224,'expanded_primary_cells':8174,'dual114_reference_only':True,'nominal_P_threshold':0.05,'FDR_recorded_not_entry_gate':True,'external_bulk_not_required':True},
'not_customer_thresholds_and_not_to_add':{'UMI_min':None,'genes_min':None,'MT_max':None,'marker_count_for_CAF_entry':None,'signature_coverage_min':None,'normal_reference_min':None,'heldout_reference_min':None,'heldout_pass_fraction':None,'contamination_exclusion_cutoff':None,'atlas_prediction_confidence_cutoff':None},
'software_parameters':'Use registered existing analysis settings or documented software defaults/recommendations; record exact versions/formals; never label them customer thresholds.',
'resource_controls_only':{'R_jobs':2,'threads_per_R_job':4,'available_memory_reserve_GiB':50,'working_memory_admission_estimate_GiB':32},
'no_automatic_membership_mutation':'All correction/sensitivity outputs align to frozen cell IDs; keep source layers and raw counts. Evaluate changes as new version, not auto-delete candidates.',
'completion':'Exit code + stage markers + full scheduled outcomes + cell/feature alignment; failed jobs remain explicit, do not silently weaken criteria to pass.'}
(O/'THRESHOLD_REGISTRY.json').write_text(json.dumps(registry,ensure_ascii=False,indent=2))
plan={'created_at':now,'objective':'Complete audited customer work within supported and expanded layers without extra biological thresholds.',
'work_packages':[
 {'id':'W1','name':'Clinical comparisons and myCAF internal patient-level DE','acceptance':'A/B1/B2/any-metastasis/site/treatment comparisons with source-valid labels; paired bone-lung effects; full patient-level genes/programs; unavailable comparisons explicitly listed'},
 {'id':'W2','name':'Current-scope quality and identity checks','acceptance':'107 DecontX outcomes and candidate marker before/after; targeted missing doublet evidence; full512 stress and formal S/G2M; sample-color plots, normal-reference and mural checks'},
 {'id':'W3','name':'Lung reference phenotype and reference sensitivity','acceptance':'HLCA AdvF/AlvF mapping and Gong program with source coverage; paired/sampled phenotype proportions; alternate same-sample immune reference sensitivity without imposed reference counts'},
 {'id':'W4','name':'Flow validation design and integrated report','acceptance':'New candidate identity/protein/antibody evidence; negative exclusion plus multimerker combinations and controls; complete audit table; no false metastasis or protein-validation claim'}],
'methods_auto_queue':{'decontx':107,'doublet_gap_samples':need},
'dependencies':'Independent W1 and scoring work proceeds while native correction runs. Interpretation that depends on correction/ref mapping waits for actual verified outputs; no waiting for user preference or S525.',
'status':'EXECUTION_STARTED_METHODS'}
(O/'PLAN.json').write_text(json.dumps(plan,ensure_ascii=False,indent=2))
for name in ['stress_core_author.csv','seurat_cc_updated2019.tsv']:
 src=P/'phase7_review_states_20260917/remediation_20260917/resources'/name;shutil.copy2(src,O/'resources'/name)
shutil.copy2(P/'phase27_identity_review_20260926/resources/program_panels.json',O/'resources/customer_program_panels.json')
shutil.copy2(P/'phase27_identity_review_20260926/tables/clinical_mapping_v3.tsv',O/'resources/clinical_mapping_frozen.tsv')
shutil.copy2(P/'phase27_identity_review_20260926/tables/bone_lung_pairs_v3.tsv',O/'resources/bone_lung_pairs_frozen.tsv')
log=(O/'logs/methods_controller.log').open('x')
proc=subprocess.Popen(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python','-u',str(O/'code/methods_controller.py')],stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True,cwd=O)
(O/'methods.pid').write_text(str(proc.pid)+chr(10));print('METHODS_STARTED',proc.pid,'DECONTX',107,'DOUBLET_GAP_SAMPLES',need)
