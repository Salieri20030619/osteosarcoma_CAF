from pathlib import Path
import json
from datetime import datetime,timezone
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase86_functional_bone_lung_20261012'
expected='90e65ad6b9742407ba4e0641609a292cc4b904d9989feb2f83d188e9a7c383e5'
receipt=dict(local_delivery_verified=True,local_archive='C:\\codex_test\\CAF_functional_bone_lung_20261012\\delivery\\CAF_functional_bone_lung_20261012_verified.zip',local_archive_sha256=expected,delivery_verified_at=datetime.now(timezone.utc).isoformat())
for name in ['CAF_FUNCTIONAL_BONE_LUNG_CURRENT.json','CAF_NEXT_STAGE_CURRENT.json','CAF_RESUME_STATUS_CURRENT.json','CAF_FUNCTIONAL_STATE_OBJECTIVE_CURRENT.json','CAF_ALL14_RESEARCH_CURRENT.json']:
 p=P/name;obj=json.loads(p.read_text());rec=obj if name=='CAF_FUNCTIONAL_BONE_LUNG_CURRENT.json' else obj['functional_bone_lung_analysis']
 assert rec['archive_sha256']==expected and rec['status']=='completed_verified'
 rec.update(receipt)
 if name in ['CAF_NEXT_STAGE_CURRENT.json','CAF_RESUME_STATUS_CURRENT.json']:
  for key in ['latest_task','last_completed_task','last_completed_analysis','latest_research_results']:obj[key].update(receipt)
 p.write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
done=json.loads((O/'COMPLETED.json').read_text());done.update(receipt);(O/'COMPLETED.json').write_text(json.dumps(done,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print(json.dumps(dict(status='completed_verified_and_delivered',phase=O.name,all14_retained=True,clinical_framework='Bone versus Lung',archive_sha256=expected)))
