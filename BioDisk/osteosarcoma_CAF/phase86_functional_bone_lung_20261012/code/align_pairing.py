from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='2'
import json,importlib.util,shutil
import numpy as np,pandas as pd
from statsmodels.stats.multitest import multipletests
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase86_functional_bone_lung_20261012'
spec=importlib.util.spec_from_file_location('core',P/'phase73_F134_bone_lung_20261009/strict_v2/code/analyze_bone_lung_v2.py');core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core);core.SEED=20261012
old=O/'tables/all_functional_bone_lung_statistics.tsv'
archive=O/'resources/clinical_manifest_pair_exclusion_sensitivity'
archive.mkdir(exist_ok=True)
for rel in ['tables/all_functional_bone_lung_statistics.tsv','tables/nominal_functional_bone_lung_statistics.tsv','QA/INDEPENDENT_QA.json','COMPUTATION_COMPLETE.json','resources/PROTOCOL_FROZEN.json']:
 f=O/rel;shutil.copy2(f,archive/f.name)
amendment=dict(reason='Consistency audit found original comparison excluded a clinically paired donor even when only one site had evaluable CAF. Primary analysis must reproduce the existing Phase81 endpoint-available pairing rule and scheme1 inclusion; prior stricter clinical-pair exclusion retained as sensitivity, not discarded.',primary='At each source/scope/matrix/minimum-cell endpoint, paired donors have available cells at both Bone and Lung. These paired observations are removed from independent comparison; one-site-evaluable donors contribute once to independent group.',no_duplicate_patient=True,pre_existing_policy_source='phase81_allCAF_bone_lung_20261012/resources/PROTOCOL_FROZEN.json',no_program_or_threshold_change=True)
(O/'resources/PAIRING_CONTINUITY_AMENDMENT.json').write_text(json.dumps(amendment,ensure_ascii=False,indent=2),encoding='utf-8')
p=pd.read_csv(O/'tables/patient_site_programs.tsv',sep='\t');rows=[]
for (scope,ds,matrix,pr),g in p.groupby(['scope','dataset','matrix','program'],sort=False):
 for minimum in [1,5]:
  gg=g[g.cells>=minimum];b=gg[gg.site=='Bone'].set_index('donor_record');l=gg[gg.site=='Lung'].set_index('donor_record');pair=sorted(set(b.index)&set(l.index))
  for design in ['independent','paired']:
   ib=pair if design=='paired' else sorted(set(b.index)-set(pair));il=pair if design=='paired' else sorted(set(l.index)-set(pair))
   for endpoint in ['mean_score','high_fraction']:
    x=b.loc[ib,endpoint].to_numpy().reshape(-1,1);y=l.loc[il,endpoint].to_numpy().reshape(-1,1);st=core.stats(x,y,design=='paired').iloc[0].to_dict()
    if design=='paired':st['mean_lung_minus_bone']=core.mean(y-x)[0]
    effect=st['mean_lung_minus_bone']
    if design=='paired':loo=[float(np.mean(np.delete(y-x,i,axis=0))) for i in range(len(x)) if len(x)>1]
    else:loo=[float(np.mean(y)-np.mean(np.delete(x,i,axis=0))) for i in range(len(x)) if len(x)>1 and len(y)]+[float(np.mean(np.delete(y,i,axis=0))-np.mean(x)) for i in range(len(y)) if len(y)>1 and len(x)]
    st.update(scope=scope,dataset=ds,matrix=matrix,program=pr,minimum_cells=minimum,design=design,endpoint=endpoint,bone_mean=float(x.mean()) if x.size else np.nan,lung_mean=float(y.mean()) if y.size else np.nan,bone_cells=int(b.loc[ib,'cells'].sum()),lung_cells=int(l.loc[il,'cells'].sum()),bone_patients=';'.join(ib),lung_patients=';'.join(il),LOPO_min=min(loo) if loo else np.nan,LOPO_max=max(loo) if loo else np.nan,LOPO_direction_retained=bool(all(v*effect>0 for v in loo)) if loo else None,family=f'{scope}|{ds}|{matrix}|{minimum}|{design}|{endpoint}',pairing_policy='Phase81_endpoint_available')
    rows.append(st)
f=pd.DataFrame(rows);f['FDR']=np.nan
for _,g in f.groupby('family'):
 ix=g.index[g.P.notna()]
 if len(ix):f.loc[ix,'FDR']=multipletests(f.loc[ix,'P'],method='fdr_bh')[1]
f['nominal']=f.P.lt(.05);f.to_csv(old,sep='\t',index=False,na_rep='NA');f[f.nominal].to_csv(O/'tables/nominal_functional_bone_lung_statistics.tsv',sep='\t',index=False,na_rep='NA')
q=O/'resources/PROTOCOL_FROZEN.json';v=json.loads(q.read_text());v['independence']=amendment['primary'];v['pairing_amendment']=str(O/'resources/PAIRING_CONTINUITY_AMENDMENT.json');q.write_text(json.dumps(v,ensure_ascii=False,indent=2),encoding='utf-8')
print('PRIMARY_PAIRING_ALIGNED',f[f.nominal][['scope','dataset','matrix','program','minimum_cells','design','endpoint','bone_N','lung_N','mean_lung_minus_bone','P','FDR']].to_string(index=False))
