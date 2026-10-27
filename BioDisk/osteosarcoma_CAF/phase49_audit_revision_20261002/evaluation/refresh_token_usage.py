from pathlib import Path
import json,datetime
paths=[Path('C:/Users/31937/.codex/sessions/2026/09/30/rollout-2026-09-30T21-18-20-01a0f50a-7b58-7780-8a17-c92abc352d45.jsonl'),Path('C:/Users/31937/.codex/sessions/2026/10/02/rollout-2026-10-02T09-01-18-01a0fa21-3f57-7e03-9d9a-3c15b85585e8.jsonl')]
records=[]
for p in paths:
 turns={};cur=None
 for line in p.read_text(encoding='utf-8').splitlines():
  try:o=json.loads(line)
  except json.JSONDecodeError:continue
  payload=o.get('payload',{})
  if o.get('type')=='turn_context':
   cur=payload.get('turn_id')
   if cur not in turns:turns[cur]={'turn_id':cur,'model':payload.get('model'),'effort':payload.get('effort'),'start_utc':o['timestamp']}
  if cur and o.get('type')=='event_msg' and payload.get('type')=='token_count' and (payload.get('info') or {}).get('total_token_usage'):
   turns[cur]['cumulative']=payload['info']['total_token_usage'];turns[cur]['last_snapshot_utc']=o['timestamp']
 prev={}
 for r in turns.values():
  if 'cumulative' not in r:continue
  r['turn_usage']={k:v-prev.get(k,0) for k,v in r['cumulative'].items()};r['turn_usage']['uncached_input_tokens']=r['turn_usage']['input_tokens']-r['turn_usage']['cached_input_tokens'];r['counter_nonnegative']=all(v>=0 for v in r['turn_usage'].values());assert r['counter_nonnegative'];u=r['turn_usage'];assert u['input_tokens']+u['output_tokens']==u['total_tokens'];assert 0<=u['cached_input_tokens']<=u['input_tokens'];assert 0<=u['reasoning_output_tokens']<=u['output_tokens'];r['source_file']=str(p);prev=r['cumulative'];records.append(r)
out={'captured_at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope':'main conversation counters; auto-review omitted; current repair turn is interim','output_includes_reasoning':True,'cached_is_subset_of_input':True,'records':records};Path('C:/Users/31937/CAF_revision_20261002/TOKEN_USAGE_BY_TURN.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8');print('TOKEN_USAGE_REFRESHED',out['captured_at'])
