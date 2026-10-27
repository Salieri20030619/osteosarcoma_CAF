import pathlib,json,itertools,math
import numpy as np,pandas as pd
from statsmodels.stats.multitest import multipletests
R=pathlib.Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');P=R/'phase105_original_rule_NCI_20261020'
a=np.load(R/'phase104_full_gene_four_models_20261020/objects/tumor_common_patient_expression.npz')
Y=a['expression'];src=a['sources'];site=a['sites'];features=a['features'];fi={v:i for i,v in enumerate(features)}
qa=[]
for strategy in ['augmented','NCI_only']:
 for source in ['SCPCP000017','SCPCP000023']:
  f=P/'tables'/(strategy+'_'+source+'_allgenes.tsv.gz');d=pd.read_csv(f,sep='\t').set_index('feature')
  b=Y[(src==source)&(site=='Bone')];mask=(src=='NCI')|((src==source)&(site=='Lung')) if strategy=='augmented' else src=='NCI';l=Y[mask]
  assert np.allclose(d.loc[features,'mean_lung_minus_bone'],l.mean(0)-b.mean(0),atol=1e-12)
  valid=d[d.P.notna()]
  assert np.allclose(multipletests(valid.P,method='fdr_bh')[1],valid.FDR,atol=1e-12)
  selected=valid.sort_values('P').iloc[np.unique(np.r_[np.arange(4),np.linspace(0,len(valid)-1,7,dtype=int)])].index
  idx=[fi[k] for k in selected];v=np.vstack([b[:,idx],l[:,idx]]);nb=len(b);nl=len(l);obs=abs(l[:,idx].mean(0)-b[:,idx].mean(0))
  exact=math.comb(nb+nl,nl)<=65536
  if exact:
   combos=np.array(list(itertools.combinations(range(nb+nl),nb)),dtype=int)
   sb=v[combos].sum(1);null=abs((v.sum(0)-sb)/nl-sb/nb)
   p=(null>=obs-1e-12).mean(0)
  else:
   hits=np.zeros(len(idx),int)
   for j in range(65536):
    ix=np.random.default_rng(20261012+j).choice(nb+nl,nl,replace=False)
    sl=v[ix].sum(0);null=abs(sl/nl-(v.sum(0)-sl)/nb);hits+=null>=obs-1e-12
   p=(hits+1)/65537
  error=float(np.max(abs(p-d.loc[selected,'P'].to_numpy())))
  assert error<1e-12,(strategy,source,error)
  boundary=valid[(valid.P-.05).abs()<.002]
  qa.append(dict(strategy=strategy,source=source,independent_recalculated_genes=len(idx),max_P_error=error,exact=exact,assignments=math.comb(nb+nl,nl) if exact else 65536,near_P_threshold_within_0_002=len(boundary)))
  boundary.reset_index().to_csv(P/'tables'/(strategy+'_'+source+'_P_boundary.tsv'),sep='\t',index=False,encoding='utf-8')
  print('QA',qa[-1],flush=True)
for strategy in ['augmented','NCI_only']:
 z=pd.read_csv(P/'tables'/(strategy+'_intersection_allgenes.tsv.gz'),sep='\t')
 independent=(z.P_17<.05)&(z.P_23<.05)&(z.mean_lung_minus_bone_17*z.mean_lung_minus_bone_23>0)
 assert np.array_equal(independent,z.selected)
(P/'resources/INDEPENDENT_QA.json').write_text(json.dumps(dict(status='PASS',checks=qa,original_seed=20261012,all_intersections_verified=True),indent=2),encoding='utf-8')
