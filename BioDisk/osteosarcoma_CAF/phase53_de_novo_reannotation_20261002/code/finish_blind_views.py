from pathlib import Path
import json,hashlib,hmac,gc
import pandas as pd,numpy as np
R=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase53_de_novo_reannotation_20261002');B=R/'blind_review';U=R/'unblinding'
maps=json.loads((U/'identifier_maps.json').read_text());cfg=json.loads((U/'anonymization_secret.json').read_text());salt=bytes.fromhex(cfg['salt'])
def cid(x):return 'C_'+hmac.new(salt,str(x).encode(),hashlib.sha256).hexdigest()[:24]
out=B/'regional_CNV';out.mkdir(exist_ok=True)
for f in sorted((R/'objects/infercnv_regions').glob('*')):
 matches=[sm for sm in maps['samples'] if f.name.startswith(sm+'__')];assert len(matches)==1
 sample=matches[0];tail=f.name[len(sample)+2:];samplecode=maps['samples'][sample]
 context=tail.replace('_groups.tsv.gz','').replace('_regions.tsv.gz','')
 t=pd.read_csv(f,sep='\t',float_precision='round_trip')
 if 'cell_id' in t:t['cell_id']=t.cell_id.map(cid)
 if 'cell_group_name' in t:
  t['group_id']=t.cell_group_name.map(lambda v:'G_'+hmac.new(salt,(sample+'|'+context+'|'+str(v)).encode(),hashlib.sha256).hexdigest()[:12]);t=t.drop(columns=['cell_group_name'])
 if 'cnv_name' in t:t=t.drop(columns=['cnv_name'])
 t.to_csv(out/(samplecode+'__'+tail),sep='\t',index=False)
cm=pd.read_csv(R/'tables/copykat_cell_continuous_metrics.tsv.gz',sep='\t')
cm['cell_id']=cm.cell_id.map(cid);cm['sample_id']=cm['sample'].map(maps['samples']);cm=cm.drop(columns='sample')
cm.to_csv(B/'copykat_cell_continuous_metrics.tsv.gz',sep='\t',index=False)
# Supply full library totals for independent normalization of the cluster count matrices.
key=pd.read_csv(U/'cluster_identifier_key.tsv',sep='\t',dtype={'original_cluster':str})
for ds,source_id in maps['sources'].items():
 obs=pd.read_csv(R/'tables'/ds/'cell_structure.tsv.gz',sep='\t',usecols=['r10','full_library'],dtype={'r10':str})
 lib=obs.groupby('r10').full_library.sum();km=key[(key.dataset==ds)&(key.resolution=='r10')].set_index('cluster_id').original_cluster
 f=B/source_id/'clusters.tsv';t=pd.read_csv(f,sep='\t');t['full_library']=t.cluster_id.map(km).map(lib)
 assert t.full_library.notna().all();t.to_csv(f,sep='\t',index=False)
figs=sorted(str(f.relative_to(B)) for d in maps['sources'].values() for f in (B/d).glob('*.png'))
pd.DataFrame({'figure':figs}).to_csv(B/'FIGURE_INDEX.tsv',sep='\t',index=False)
roster=pd.read_csv(B/'anonymous_cell_cluster_roster.tsv.gz',sep='\t',usecols=['cell_id','source_id'])
template=pd.read_csv(B/'reviewer_response_template.tsv',sep='\t')
assert len(roster)==876648 and roster.cell_id.is_unique and len(template)==307
qa={'primary_clusters':len(template),'cells':len(roster),'sources':len(maps['sources']),'figures':len(figs),'proposed_labels_not_read_by_exporter':True,'original_identifiers_mapped':True,'sample_names_with_double_underscores_handled':True,'genomic_group_ids_scoped_by_sample_and_method':True,'full_library_metadata_included':True}
(R/'delivery/BLIND_VIEW_QA.json').write_text(json.dumps(qa,indent=2));(R/'BLIND_VIEWS_COMPLETE').write_text('SUCCESS');print(json.dumps(qa),flush=True)
