from pathlib import Path
import sys,os,subprocess,json,datetime,hashlib,select,gzip,io,traceback
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');Q=P/'phase48_clinical_programs_delivery_20261001';D=Q/'delivery/client';O=P/'phase46_spatial_functional_20261001';E=P/'phase47_external_candidate_integration_20261001'
def put(path,obj):path.write_text(json.dumps(obj,ensure_ascii=False,indent=2,default=str))
def state(s,**kw):put(Q/'FINAL_DELIVERY_STATUS.json',{'state':s,'at':datetime.datetime.now().astimezone().isoformat(),**kw});print(s,kw,flush=True)
def run(name,cmd,allow_failure=False):
 with (Q/'logs'/f'{name}.log').open('w') as f:
  p=subprocess.Popen(cmd,stdout=f,stderr=subprocess.STDOUT);(Q/f'{name}.pid').write_text(str(p.pid));r=p.wait()
 (Q/f'{name}.exit_code').write_text(str(r));assert r==0 or allow_failure,(name,r)
if '--start' in sys.argv:
 f=(Q/'logs/finish_delivery.log').open('w');p=subprocess.Popen([sys.executable,__file__],stdout=f,stderr=subprocess.STDOUT,start_new_session=True,stdin=subprocess.DEVNULL);(Q/'finish_delivery.pid').write_text(str(p.pid));print('NATIVE_DELIVERY_CONTROLLER_STARTED',p.pid)
else:
 rc=1
 try:
  state('waiting_for_native_background_image_download')
  # Native file completion events work on this older Linux/Python build.
  import ctypes
  lib=ctypes.CDLL(None,use_errno=True);fd=lib.inotify_init1(0);assert fd>=0
  watch=lib.inotify_add_watch(fd,os.fsencode(Q),8);assert watch>=0
  while True:
   try:
    st=json.loads((Q/'HIRES_DOWNLOAD_STATUS.json').read_text())
    if (Q/'hires_parallel_native.exit_code').exists() and (Q/'hires_parallel_native.exit_code').stat().st_mtime>=(Q/'hires_parallel_worker.pid').stat().st_mtime:break
   except (FileNotFoundError,json.JSONDecodeError):pass
   os.read(fd,65536)
  os.close(fd)
  # A missing process is never treated as successful; use exit status and outcomes.
  for attempt in range(3):
   outcomes=json.loads((Q/'HIRES_DOWNLOAD_STATUS.json').read_text());code=(Q/'hires_download.exit_code').read_text().strip()
   if code=='0' and outcomes['state']=='complete':break
   state('retry_failed_image_resources',attempt=attempt+1);run('hires_retry_'+str(attempt+1),[sys.executable,str(Q/'code/download_hires_parallel.py'),'--worker'],allow_failure=True)
  outcomes=json.loads((Q/'HIRES_DOWNLOAD_STATUS.json').read_text());assert outcomes['state']=='complete' and len(outcomes['outcomes'])==26 and (Q/'hires_download.exit_code').read_text().strip()=='0'
  from PIL import Image
  import pandas as pd,anndata as ad,numpy as np
  rows=[]
  for r in outcomes['outcomes']:
   f=Path(r['path']);raw=f.read_bytes();assert hashlib.sha256(raw).hexdigest()==r['SHA256'];im=Image.open(io.BytesIO(gzip.decompress(raw)));im.load();a=ad.read_h5ad(P/'phase35_external_validation_prep_20260927/objects/GSE293065'/f'{r["sample"]}_counts_spatial.h5ad',backed='r');xy=np.asarray(a.obsm['spatial'])*float(a.uns['spatial_scalefactors']['tissue_hires_scalef']);ok=(xy[:,0]>=0)&(xy[:,0]<im.width)&(xy[:,1]>=0)&(xy[:,1]<im.height);assert ok.all();rows.append({**r,'spots':a.n_obs,'in_image':int(ok.sum()),'independent_tumor_mask':False});a.file.close()
  pd.DataFrame(rows).to_csv(Q/'tables/hires_image_coordinate_audit.tsv',sep='\t',index=False);put(Q/'delivery/HIRES_FINAL_ACCEPTANCE.json',{'images':26,'coordinates_in_bounds':sum(r['in_image'] for r in rows),'all_bytes_SHA_gzip_PNG_checked':True,'proxy':False,'geometry_still_unassessable_independent_mask_absent':True})
  state('image_resources_accepted; waiting_for_report_review_gate')
  # The review gate is an explicit artefact generated only after all final pages are inspected.
  while not (D/'REPORT_VISUAL_QA.json').exists():
   import time;time.sleep(15)
  review=json.loads((D/'REPORT_VISUAL_QA.json').read_text());assert review['all_pages_visual_passed'] and review['pages']==16
  markers=[O/'SPATIAL_COMPLETE',O/'INDEPENDENT_R_QA_COMPLETE',O/'MODEL_SENSITIVITY_REPAIR_QA_COMPLETE',E/'EXTERNAL_COMPLETE',E/'SOURCE_PROJECTION_COMPLETE',E/'INTEGRATION_COMPLETE',Q/'CLINICAL_PROGRAM_COMPLETE',Q/'INDEPENDENT_R_QA_COMPLETE',Q/'SURFACE_COMPLETE',Q/'SURFACE_INDEPENDENT_R_QA_COMPLETE',Q/'HISTOLOGY_AUDIT_COMPLETE',Q/'HIRES_DOWNLOAD_COMPLETE']
  assert all(f.read_text().strip()=='SUCCESS' for f in markers)
  checks=[O/'controller.exit_code',E/'controller.exit_code',Q/'QA_controller.exit_code',Q/'surface_controller.exit_code',Q/'surface_R_QA.exit_code',Q/'histology_audit.exit_code',Q/'fraction_HC3_audit_v3.exit_code']
  assert all(f.read_text().strip()=='0' for f in checks)
  conditions=[('tumor_boundary_geometry','unassessable','independent tumor/region mask absent;26lowres and26hires images available'),('real_immune_inhibition_or_protection','not_established','functional experiments not part of computational completion'),('protein_sorting_gate','unvalidated','RNA-based two-marker candidates and negative controls proposed'),('cross_source_global_patient_dedup','audit_complete_correspondence_unresolved','no unknown source pooling or independent confirmation claims'),('OsteoCAR_Budhathoki','unobtained_historical_source_gap','403/access gaps retained; available cohorts analyzed'),('computational_perturbation','condition_not_triggered','optional in frozen plan; no validated intervention model answering key matrix mechanism')]
  pd.DataFrame(conditions,columns=['endpoint','status','reason']).to_csv(D/'CONDITIONAL_ENDPOINTS.tsv',sep='\t',index=False)
  acceptance={'at':datetime.datetime.now().astimezone().isoformat(),'state':'all_evaluable_computational_modules_and_report_accepted','scope':'full detectable RNA transcriptome, source-scoped patient tests, functional/TF/LR discovery, spatial fallback, external expression migration, marker candidates and report','success_markers':[str(f) for f in markers],'exit_codes':{str(f):0 for f in checks},'report_QA':review,'hires_images':26,'lowres_images':26,'histology_boundary_accepted':False,'source_patient_pooling':False,'conditional_endpoints':conditions,'preferred_fraction_models':str(Q/'tables/D_fraction_HC3_current_statistics.tsv'),'preferred_continuous_clinical':str(Q/'tables/clinical_program_statistics.tsv'),'historical_failures_retained':True}
  put(D/'FINAL_ACCEPTANCE.json',acceptance);put(Q/'delivery/FINAL_ACCEPTANCE.json',acceptance)
  state('building_and_checking_full_statistics_archive');run('statistical_package',[sys.executable,str(Q/'code/build_statistical_package.py')]);integrity=json.loads((D/'PACKAGE_INTEGRITY.json').read_text());assert integrity['CRC_passed'] and (Q/'STATISTICAL_PACKAGE_COMPLETE').read_text()=='SUCCESS'
  previous=P/'CAF_NEXT_STAGE_CURRENT.json'
  if previous.exists():(Q/'delivery/CURRENT_BEFORE_FINAL_DELIVERY.json').write_bytes(previous.read_bytes())
  current={'at':acceptance['at'],'state':'complete_evaluable_next_stage_analysis','phase':'phase41_to_phase48; client delivery20261002','report':str(D/'CAF客户结果报告.pdf'),'statistics_archive':integrity,'acceptance':str(D/'FINAL_ACCEPTANCE.json'),'conditional_endpoints':conditions,'no_proxy':True,'cross_source_patient_independence_not_fabricated':True}
  put(previous,current);put(Q/'delivery/FINAL_PROJECT_PROGRESS.json',current);(Q/'FINAL_DELIVERY_COMPLETE').write_text('SUCCESS');state('complete',archive=integrity);rc=0
 except Exception as e:traceback.print_exc();state('failed',error=repr(e))
 (Q/'finish_delivery.exit_code').write_text(str(rc));print('FINAL_DELIVERY_EXIT',rc,flush=True)
