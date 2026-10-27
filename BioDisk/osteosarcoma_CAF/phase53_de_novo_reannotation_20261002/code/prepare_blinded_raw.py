from pathlib import Path
import os,json,secrets,hmac,hashlib,random,gc
import numpy as np,pandas as pd,anndata as ad,h5py
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');R=P/'phase53_de_novo_reannotation_20261002'
m=pd.read_csv(P/'phase41_full_transcriptome_setup_20261001/tables/clinical_sample_registry.tsv',sep='\t',keep_default_na=False)
B=R/'blind_review';U=R/'unblinding';(B/'raw_counts').mkdir(parents=True,exist_ok=True)
secret=U/'anonymization_secret.json'
if secret.exists():cfg=json.loads(secret.read_text())
else:
 cfg={'salt':secrets.token_hex(32),'seed':secrets.randbits(64)};secret.write_text(json.dumps(cfg))
rng=random.Random(cfg['seed'])
def keymap(values,prefix,width):
 v=sorted(set(values));out=[prefix+str(i+1).zfill(width) for i in range(len(v))];rng.shuffle(out);return dict(zip(v,out))
sources=keymap(m.dataset,'D',2);samples=keymap(m['sample'],'S',3);patients=keymap(m.donor_record,'P',3)
maps={'sources':sources,'samples':samples,'patients':patients};(U/'identifier_maps.json').write_text(json.dumps(maps,indent=2))
def cellid(x):return 'C_'+hmac.new(bytes.fromhex(cfg['salt']),str(x).encode(),hashlib.sha256).hexdigest()[:24]
technical=[];keys=[];manifest=[];seen=set()
cw=pd.read_csv(P/'phase43_full_transcriptome_discovery_20261001/tables/HGNC_full_gene_symbol_crosswalk.tsv.gz',sep='\t',keep_default_na=False);sym=dict(zip(cw.feature,cw.canonical_symbol))
for row in m.itertuples(index=False):
 src=P/'phase26_client_reanalysis_20260922/inputs'/row.sample/'full_called_counts.h5ad';a=ad.read_h5ad(src);original=a.obs_names.astype(str);ids=np.array([cellid(x) for x in original])
 assert len(set(ids))==len(ids) and not set(ids)&seen;seen.update(ids)
 obs=pd.DataFrame(index=pd.Index(ids,name='cell_id'));obs['source_id']=sources[row.dataset];obs['sample_id']=samples[row.sample];obs['patient_id']=patients[row.donor_record];obs['tissue']=row.site;obs['technology']=row.technology
 obs['library_counts']=np.asarray(a.X.sum(1)).ravel();obs['detected_features']=a.X.getnnz(axis=1)
 var=pd.DataFrame(index=pd.Index(a.var_names.astype(str),name='feature_id'));var['gene_symbol']=[sym.get(str(f),'') or str(a.var.loc[f,'symbol']) for f in a.var_names]
 b=ad.AnnData(a.X.tocsr(),obs=obs,var=var);b.uns['content']='Provider-called raw counts; identifiers anonymized; no proposed or historical lineage labels'
 dest=B/'raw_counts'/f'{samples[row.sample]}.h5ad'
 b.write_h5ad(dest,compression='gzip')
 with h5py.File(dest,'r') as h:
  assert np.array_equal(h['X/data'][:],a.X.data)
  assert np.array_equal(h['X/indices'][:],a.X.indices)
  assert np.array_equal(h['X/indptr'][:],a.X.indptr)
  assert set(h['obs'].keys())=={'cell_id','source_id','sample_id','patient_id','tissue','technology','library_counts','detected_features'}
 keys.append(pd.DataFrame({'original_cell_id':original,'cell_id':ids,'original_sample':row.sample,'sample_id':samples[row.sample],'original_dataset':row.dataset,'source_id':sources[row.dataset],'original_patient':row.donor_record,'patient_id':patients[row.donor_record]}))
 technical.append({'source_id':sources[row.dataset],'sample_id':samples[row.sample],'patient_id':patients[row.donor_record],'tissue':row.site,'technology':row.technology,'cells':a.n_obs,'features':a.n_vars})
 manifest.append({'relative_path':str(dest.relative_to(B)),'bytes':dest.stat().st_size,'SHA256':hashlib.sha256(dest.read_bytes()).hexdigest(),'cells':a.n_obs,'features':a.n_vars,'raw_count_arrays_unchanged':True})
 print('BLINDED_RAW_COMPLETE',samples[row.sample],a.n_obs,flush=True);del a,b;gc.collect()
pd.concat(keys,ignore_index=True).to_csv(U/'cell_identifier_key.tsv.gz',sep='\t',index=False)
pd.DataFrame(technical).sort_values('sample_id').to_csv(B/'sample_technical_metadata.tsv',sep='\t',index=False)
pd.DataFrame(manifest).to_csv(B/'RAW_COUNT_MANIFEST.tsv',sep='\t',index=False)
assert sum(x['cells'] for x in manifest)==876648 and len(manifest)==107
(R/'delivery/RAW_BLINDING_QA.json').write_text(json.dumps({'samples':107,'cells':876648,'no_original_identifiers_in_raw_obs':True,'all_count_arrays_verified_unchanged':True,'identifier_key_separate':True,'full_raw_counts_server_only_compact_archive_has_manifest':True},indent=2))
(R/'RAW_BLINDING_COMPLETE').write_text('SUCCESS');print('RAW_BLINDING_COMPLETE',flush=True)
