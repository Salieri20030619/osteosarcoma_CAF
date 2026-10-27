from pathlib import Path
import json
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase83_upstream_and_direction_20261012'
delivery=dict(local_delivery_verified=True,archive_sha256='78f6fdacaae9f07ac79ce8209538a595123d4c4db52b80ad66b3b070f5677d5a',local_archive='C:/codex_test/CAF_upstream_review_20261012/delivery/CAF_upstream_direction_literature_20261012.zip',local_literature_report='C:/codex_test/CAF_upstream_review_20261012/delivery/14个CAF候选marker_文献与方向联合评估_20261012.txt',local_upstream_report='C:/codex_test/CAF_upstream_review_20261012/delivery/CAF上游核查与跨来源方向复核_20261012.txt')
def write(f,d):
    tmp=f.with_suffix('.json.tmp');tmp.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');tmp.replace(f)
for f in [P/'CAF_UPSTREAM_DIRECTION_LITERATURE_CURRENT.json',O/'COMPLETE.json']:
    d=json.loads(f.read_text());assert d['archive_sha256']==delivery['archive_sha256'];d.update(delivery);write(f,d)
for name in ['CAF_NEXT_STAGE_CURRENT.json','CAF_RESUME_STATUS_CURRENT.json']:
    f=P/name;d=json.loads(f.read_text());assert d['current_phase']==O.name
    for key in ['latest_task','last_completed_task','last_completed_analysis','latest_research_results']:
        assert d[key]['archive_sha256']==delivery['archive_sha256'];d[key].update(delivery)
    d.update(local_delivery_verified=True,local_upstream_direction_delivery=delivery['local_archive']);write(f,d)
write(O/'delivery/LOCAL_DELIVERY.json',delivery)
print('Phase83 upstream, direction and literature review delivered and verified.')
