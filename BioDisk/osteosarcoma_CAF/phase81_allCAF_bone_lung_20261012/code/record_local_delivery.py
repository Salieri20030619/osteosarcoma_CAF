from pathlib import Path
import json,hashlib
from datetime import datetime,timezone
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase81_allCAF_bone_lung_20261012'
expected='5c934a605ace244fb3ce3fe9d923797e743c4df6f384db218a7694c75e90a631'
delivery=dict(local_delivery_verified=True,local_archive='C:/codex_test/CAF_allCAF_bone_lung_20261012/delivery/CAF_allCAF_bone_lung_20261012.zip',local_report='C:/codex_test/CAF_allCAF_bone_lung_20261012/delivery/全部CAF候选_骨肺病灶比较_20261012.txt',archive_sha256=expected,local_verified_at=datetime.now(timezone.utc).isoformat())
def write(f,d):
    tmp=f.with_suffix('.json.tmp');tmp.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');tmp.replace(f)
assert hashlib.sha256((O/'delivery/CAF_allCAF_bone_lung_20261012.zip').read_bytes()).hexdigest()==expected
for f in [P/'CAF_ALLCAF_BONE_LUNG_CURRENT.json',P/'CAF_SCHEME1_CURRENT.json',O/'COMPLETE.json']:
    d=json.loads(f.read_text());assert d['archive_sha256']==expected;d.update(delivery);write(f,d)
for name in ['CAF_NEXT_STAGE_CURRENT.json','CAF_RESUME_STATUS_CURRENT.json']:
    f=P/name;d=json.loads(f.read_text());assert d['current_phase']==O.name
    for key in ['latest_task','last_completed_task','last_completed_analysis','latest_research_results']:
        assert d[key]['archive_sha256']==expected;d[key].update(delivery)
    d.update(local_delivery_verified=True,local_allCAF_delivery=delivery['local_archive']);write(f,d)
write(O/'delivery/LOCAL_DELIVERY.json',delivery)
print('All-CAF analysis and verified local delivery recorded.')
