from pathlib import Path
import json
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase82_allCAF_crossdataset_markers_20261012'
delivery=dict(local_delivery_verified=True,local_archive='C:/codex_test/CAF_crossdataset_markers_20261012/delivery/CAF_allCAF_crossdataset_markers_20261012.zip',local_report='C:/codex_test/CAF_crossdataset_markers_20261012/delivery/整体CAF_跨数据集marker筛选_20261012.txt',local_marker_table='C:/codex_test/CAF_crossdataset_markers_20261012/delivery/MARKER_EVIDENCE_REGISTRY.tsv',archive_sha256='cbebbb0125cd68b3c6beec2b4e98f4b0556fa7fc7b1b6f73004dac043e08a181')
def write(f,d):
    tmp=f.with_suffix('.json.tmp');tmp.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');tmp.replace(f)
for f in [P/'CAF_ALLCAF_CROSSDATASET_MARKERS_CURRENT.json',O/'COMPLETE.json']:
    d=json.loads(f.read_text());assert d['archive_sha256']==delivery['archive_sha256'];d.update(delivery);write(f,d)
for name in ['CAF_NEXT_STAGE_CURRENT.json','CAF_RESUME_STATUS_CURRENT.json']:
    f=P/name;d=json.loads(f.read_text());assert d['current_phase']==O.name
    for key in ['latest_task','last_completed_task','last_completed_analysis','latest_research_results']:
        assert d[key]['archive_sha256']==delivery['archive_sha256'];d[key].update(delivery)
    d.update(local_delivery_verified=True,local_crossdataset_marker_delivery=delivery['local_archive']);write(f,d)
write(O/'delivery/LOCAL_DELIVERY.json',delivery)
print('Phase82 marker analysis and verified local delivery recorded.')
