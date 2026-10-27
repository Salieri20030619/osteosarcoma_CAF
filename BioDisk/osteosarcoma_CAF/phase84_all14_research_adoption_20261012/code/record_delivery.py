from pathlib import Path
import json
root=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF')
receipt={'local_delivery_verified': True, 'local_archive': 'C:\\codex_test\\CAF_all14_research_20261012\\delivery\\CAF_all14_research_adoption_20261012.zip', 'archive_sha256': '4d0c91a596d2e5a3d6c7ce47f1c9bbbdfd50a6bde2e5f78373a2b23636bf2f4d', 'local_files': ['C:\\codex_test\\CAF_all14_research_20261012\\delivery\\ALL14_NEXT_STAGE_RESEARCH_REGISTRY.json', 'C:\\codex_test\\CAF_all14_research_20261012\\delivery\\ALL14_NEXT_STAGE_RESEARCH_REGISTRY.tsv', 'C:\\codex_test\\CAF_all14_research_20261012\\delivery\\DATA_SOURCE_ROLES.tsv', 'C:\\codex_test\\CAF_all14_research_20261012\\delivery\\RESEARCH_SELECTION_POLICY.json', 'C:\\codex_test\\CAF_all14_research_20261012\\delivery\\14个CAF候选_全部纳入下一阶段研究_20261012.txt'], 'verified_at': '2026-10-12T07:26:37.128525+00:00'}
names=['CAF_ALL14_RESEARCH_CURRENT.json','CAF_NEXT_STAGE_CURRENT.json','CAF_RESUME_STATUS_CURRENT.json','CAF_CLASSIFICATION_DIFFERENCE_RESUME_CURRENT.json','CAF_ALLCAF_CROSSDATASET_MARKERS_CURRENT.json','CAF_UPSTREAM_DIRECTION_LITERATURE_CURRENT.json']
for name in names:
    path=root/name;obj=json.loads(path.read_text(encoding='utf-8'))
    rec=obj if name=='CAF_ALL14_RESEARCH_CURRENT.json' else obj['all14_research_adoption']
    assert rec['retained_count']==14 and rec['excluded_count']==0
    assert rec['archive_sha256']==receipt['archive_sha256']
    assert rec['selection_policy']['record_only_sources_affect_ranking'] is False
    rec.update(receipt)
    if name in ['CAF_NEXT_STAGE_CURRENT.json','CAF_RESUME_STATUS_CURRENT.json']:
        for key in ['latest_task','last_completed_task','last_completed_analysis','latest_research_results']:obj[key].update(receipt)
    path.write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print(json.dumps(dict(status='completed_verified',retained_count=14,local_delivery_verified=True,updated_records=len(names))))
