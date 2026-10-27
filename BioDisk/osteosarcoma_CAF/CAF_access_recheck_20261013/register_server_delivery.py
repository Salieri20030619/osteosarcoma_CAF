from pathlib import Path
import json,hashlib,zipfile,datetime
R=Path(__file__).resolve().parent
delivery=R/'delivery';receipt=json.loads((delivery/'PACKAGE_RECEIPT.json').read_text(encoding='utf8'))
p=delivery/receipt['package'];assert p.stat().st_size==receipt['bytes'] and hashlib.sha256(p.read_bytes()).hexdigest()==receipt['sha256']
target=R/'review_bundle';target.mkdir(exist_ok=True)
with zipfile.ZipFile(p) as z:
 for name in z.namelist():assert (target/name).resolve().is_relative_to(target.resolve())
 z.extractall(target)
verified=json.loads((R/'SERVER_ACCEPTANCE.json').read_text());assert verified['status']=='PASS'
pointer={'status':'completed_available_sources_and_access_audit','updated_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'task_directory':str(R),'report':str(target/'CAF_NCC_HRA及相关文献复核报告_20261013.txt'),'delivery_archive':str(p),'NCC_downloaded_bytes':487493626,'NCC_single_cell_files':10,'NCC_labeled_single_cell_files':9,'HRA013407_official_release_state':'Not released','HRA013407_bioproject':'PRJCA032564','project_ReleaseDate_field':'2028-08-19','Shen_expression_accession':'unresolved','new_formal_CAF_patients':0,'baseline_unchanged':True}
(R.parent/'CAF_SOURCE_ACCESS_RECHECK_CURRENT.json').write_text(json.dumps(pointer,ensure_ascii=False,indent=2),encoding='utf8')
(delivery/'SERVER_DELIVERY_COMPLETE.json').write_text(json.dumps({'status':'PASS','sha256':receipt['sha256'],'bytes':receipt['bytes'],'pointer':str(R.parent/'CAF_SOURCE_ACCESS_RECHECK_CURRENT.json')},indent=2))
print(json.dumps(pointer,ensure_ascii=False))
