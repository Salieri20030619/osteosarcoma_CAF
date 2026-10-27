from pathlib import Path
import os,sys,json,hashlib,itertools
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
import numpy as np,pandas as pd,anndata as ad
from scipy.stats import false_discovery_control
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase63_matrix_resistance_20261005'
def read(p):return pd.read_csv(p,sep='\t',low_memory=False)
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for z in iter(lambda:f.read(8388608),b''):h.update(z)
 return h.hexdigest()
for stage in ['freeze','mechanism','spatial','structure','integrate']:assert (O/f'{stage}.exit_code').read_text().strip()=='0'
orig=read(P/'phase62_F1F5_internal_20261005/tables/F1F5_internal_cell_labels.tsv.gz').set_index('cell_id');new=read(O/'tables/cell_mechanism_scores_with_confidence.tsv.gz').set_index('cell_id');assert len(new)==8174 and new.index.is_unique and orig.index.tolist()==new.index.tolist()
for c in ['class_id','round1_state','subgroup_id','working_name','lineage_confidence','functional_interpretation_confidence','expression_confidence','cross_patient_support','CNV_evidence_current','candidate_tier','mixed_program_flag','doublet_positive_any']:assert orig[c].fillna('NA').equals(new[c].fillna('NA')),c
for r in json.loads((O/'resources/INPUT_MANIFEST.json').read_text()):assert sha(Path(r['path']))==r['sha256']
# Reconstruct all program statistics from per-patient effects, including absence of P when no independent repeats.
checks=[];pchecks=0;fdrfamilies=0
def exact(v):
 n=len(v);signs=np.array(list(itertools.product([-1.,1.],repeat=n)));return float((np.abs(signs@v/n)>=abs(np.mean(v))-1e-12).mean())
for statfile,efffile,keys in [
('patient_program_statistics.tsv.gz','patient_program_effects.tsv.gz',['dataset','subgroup_id','matrix','scope','minimum_cells','program']),
('posthoc_collagen_structure_statistics.tsv','posthoc_collagen_structure_patient_effects.tsv',['dataset','subgroup_id','matrix','scope','minimum_cells','program']),
('spatial_patient_statistics.tsv','spatial_patient_effects.tsv',['version','program','endpoint','role'])]:
 s=read(O/'tables'/statfile);ef=read(O/'tables'/efffile);groups={k:g.effect.dropna().to_numpy() for k,g in ef.groupby(keys,dropna=False)};maxerr=0
 for r in s.itertuples():
  key=tuple(getattr(r,k) for k in keys);v=groups.get(key,np.array([]));assert len(v)==r.N
  if len(v):
   err=abs(v.mean()-r.mean_effect);assert err<1e-11;maxerr=max(maxerr,err);assert int((v>0).sum())==r.positive_patients
  else:assert pd.isna(r.mean_effect)
  if len(v)<2 or not np.any(v):assert pd.isna(r.P)
 sample=s[s.N.between(2,12)&s.P.notna()].iloc[::max(1,len(s)//30)].head(30)
 for r in sample.itertuples():
  v=groups[tuple(getattr(r,k) for k in keys)];assert abs(exact(v)-r.P)<1e-12;pchecks+=1
 for fam,g in s.groupby('family'):
  valid=g.P.notna();q=false_discovery_control(g.P.fillna(1).to_numpy(),method='bh');assert np.allclose(q[valid],g.loc[valid,'FDR'],atol=1e-12);fdrfamilies+=1
 checks.append({'table':statfile,'rows':len(s),'maximum_mean_reconstruction_error':maxerr})
# Spatial exact sign-flip checks for N=19 use independent bit enumerations.
sp=read(O/'tables/spatial_patient_statistics.tsv');pat=read(O/'tables/spatial_patient_effects.tsv')
for p,ep in [('collagen_biosynthesis','tumor_local_adjusted'),('collagen_biosynthesis','CD8_local_adjusted'),('collagen_crosslinking','tumor_local_adjusted')]:
 row=sp[(sp.version=='full_program')&(sp.program==p)&(sp.endpoint==ep)].iloc[0];v=pat[(pat.version=='full_program')&(pat.program==p)&(pat.endpoint==ep)].effect.dropna().to_numpy();n=len(v);hits=0
 for start in range(0,2**n,8192):
  nums=np.arange(start,min(start+8192,2**n),dtype=np.uint64);signs=2*((nums[:,None]>>np.arange(n,dtype=np.uint64))&1).astype(float)-1;hits+=int((np.abs(signs@v/n)>=abs(v.mean())-1e-12).sum())
 assert abs(hits/(2**n)-row.P)<1e-12;pchecks+=1
g=read(O/'tables/targeted_mechanism_gene_statistics.tsv.gz')
for fam,fr in g.groupby('family'):
 valid=fr.P.notna();q=false_discovery_control(fr.P.fillna(1).to_numpy(),method='bh');assert np.allclose(q[valid],fr.loc[valid,'FDR'],atol=1e-12);fdrfamilies+=1
assert g.loc[g.N.lt(2),'P'].isna().all()
# Direct count reconstruction of target genes in source/patient/site comparisons.
z=np.load(O/'objects/mechanism_counts.npz');genes=z['genes'];rr=read(O/'tables/paired_patient_site_roster.tsv.gz');gene_probes=0;countsmax=0
for ds,st,mode in [('SCPCP000017','F2_C1','corrected'),('SCPCP000017','F2_C3','corrected'),('SCPCP000017','F2_C3','raw'),('SCPCP000017','F3_C1','corrected')]:
 cached=np.load(O/'objects'/f'{ds}__{st}__{mode}__patient_gene_effects.npz');ix=[int(np.flatnonzero(genes==x)[0]) for x in ['COL1A1','COL4A1','LOXL2','PCOLCE','PXDN']];source=rr[(rr.dataset==ds)&(rr.subgroup_id==st)&(rr.matrix==mode)&(rr.scope=='all')&(rr.minimum_cells==1)]
 for k,donor in enumerate(cached['donors']):
  vals=[]
  for row in source[source.donor_record==donor].itertuples():
   mask=(orig.dataset.eq(ds)&orig.donor_record.eq(donor)&orig.site.eq(row.site)).to_numpy();left=np.flatnonzero(mask&orig.subgroup_id.eq(st).to_numpy());right=np.flatnonzero(mask&orig.subgroup_id.ne(st).to_numpy());lib=z[mode+'_library'];x=z[mode];v=np.log2(x[np.ix_(left,ix)].sum(0)/lib[left].sum()*1e4+.1)-np.log2(x[np.ix_(right,ix)].sum(0)/lib[right].sum()*1e4+.1);present=z['assayed'][np.r_[left,right]][:,ix].all(0);v[~present]=np.nan;vals.append(v)
  observed=np.nanmean(np.array(vals),axis=0);expected=cached['effects'][k,ix];assert np.allclose(observed,expected,atol=1e-12,equal_nan=True);gene_probes+=len(ix);countsmax=max(countsmax,float(np.nanmax(np.abs(observed-expected))))
# Independently re-aggregate spatial event/site/patient effects.
sf=read(O/'tables/spatial_sample_effects.tsv.gz');keys=['version','program','endpoint','role'];rebuilt=sf.groupby(['patient_id','site_group','event_month']+keys,dropna=False).effect.mean().groupby(['patient_id','site_group']+keys).mean().groupby(['patient_id']+keys).mean().reset_index();merged=rebuilt.merge(pat,on=['patient_id']+keys,validate='one_to_one',suffixes=('_new','_saved'));assert np.allclose(merged.effect_new,merged.effect_saved,atol=1e-12,equal_nan=True)
cov=read(O/'tables/spatial_program_coverage.tsv');assert cov['sample'].nunique()==26
spots=read(O/'tables/spatial_spot_scores.tsv.gz');assert len(spots)==53309 and spots.patient_id.nunique()==19 and not spots[['sample','spot_id']].duplicated().any()
# Raw spatial count checks use the previously accepted canonical feature mapping.
sys.path.insert(0,str(P/'phase46_spatial_functional_20261001/code'));from common import feature_map
sets=json.loads((O/'resources/PROGRAMS.json').read_text())['sets'];spatialprobes=0
for sm in ['PT1','MR1','LR17']:
 a=ad.read_h5ad(P/'phase35_external_validation_prep_20260927/objects/GSE293065'/f'{sm}_counts_spatial.h5ad');t=spots[spots['sample'].eq(sm)];assert a.obs_names.tolist()==t.spot_id.tolist();gm=feature_map(a);lib=np.asarray(a.X.sum(1)).ravel()
 for name in ['collagen_biosynthesis','collagen_crosslinking']:
  have=[gm[g] for g in sets[name] if g in gm];score=np.log1p(a.X[:,have].toarray()/lib[:,None]*1e4).mean(1);assert np.allclose(score,t['full_program__'+name],atol=1e-12);spatialprobes+=len(score)
files=[P/'phase41_full_transcriptome_setup_20261001/objects/all_candidate_full_genes.h5ad',P/'phase58_F1F5_round1_20261004/objects/expanded_corrected_full_genes.h5ad',P/'phase41_full_transcriptome_setup_20261001/objects/sample_feature_presence.npz',P/'phase58_F1F5_round1_20261004/resources/source_canonical_feature_map.tsv']
(O/'resources/ADDITIONAL_INPUT_HASHES.json').write_text(json.dumps([{'path':str(p),'sha256':sha(p)} for p in files],indent=2))
qa={'status':'PASS','cells':8174,'subgroups':33,'original_labels_and_confidence_preserved':True,'all_program_effects_reconstructed':checks,'independent_P_checks':pchecks,'BH_families_checked':fdrfamilies,'gene_count_effect_probes':gene_probes,'maximum_count_reconstruction_error':countsmax,'raw_spatial_score_probes':spatialprobes,'spatial_patient_aggregation_verified':True,'spatial_spots':53309,'spatial_patients':19,'spatial_slices':26,'direct_response_P_generated':False,'shell_geometry_claimed':False}
(O/'QA/ACCEPTANCE.json').write_text(json.dumps(qa,indent=2));(O/'VERIFIED').write_text('SUCCESS\n');print(json.dumps(qa,indent=2),flush=True)
st=read(O/'tables/posthoc_collagen_structure_statistics.tsv');print('F2_FIBRILLAR_SENSITIVITY');print(st[(st.dataset=='SCPCP000017')&st.subgroup_id.isin(['F2_C1','F2_C3'])&st.program.eq('fibrillar_I_III_V')][['subgroup_id','matrix','scope','minimum_cells','N','mean_effect','P','FDR']].to_string(index=False))
