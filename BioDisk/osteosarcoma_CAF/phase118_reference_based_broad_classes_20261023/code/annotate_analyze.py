from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
import numpy as np,pandas as pd,json,hashlib,datetime,itertools
from scipy.stats import false_discovery_control,spearmanr
R=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=R/'phase118_reference_based_broad_classes_20261023';B=R/'phase117_joint_restratification_20261023'
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def dump(n,x):(O/n).write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str)+'\n',encoding='utf-8')
def save(d,n):d.to_csv(O/'tables'/n,sep='\t',index=False,encoding='utf-8')
def read(n):return pd.read_csv(O/'tables'/n,sep='\t',low_memory=False)
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
def event(stage,status,**kw):
 rec=dict(utc=now(),stage=stage,status=status,**kw)
 with (O/'events.jsonl').open('a',encoding='utf-8') as f:f.write(json.dumps(rec,ensure_ascii=False)+'\n')
 print(json.dumps(rec,ensure_ascii=False),flush=True);dump('STATUS.json',rec)
event('recovery','profile_verified',profile_exit_code=int((O/'logs/profile.exit_code').read_text()),duplicate_job_started=False)
assert json.loads((O/'PROFILE_COMPLETE.json').read_text())['status']=='PASS_PROFILE_READY_FOR_ANNOTATION'
m=read('cell_registry.tsv.gz');scores=read('cell_rank_scores_1500.tsv.gz');assert scores.cell_id.equals(m.cell_id)
panels=json.loads((O/'resources/PANELS_FROZEN.json').read_text());counts=m[~m.doublet].groupby(['source','donor','joint_cluster']).size().rename('target_cells').reset_index()
patient=read('patient_panel_scores_1500.tsv').merge(counts,on=['source','donor','joint_cluster'],validate='one_to_one');adequate=patient[patient.target_cells.ge(5)];ap=adequate.groupby(['source','joint_cluster'])[list(panels)].mean().groupby('joint_cluster').mean().reset_index();save(ap,'balanced_cluster_scores_min5.tsv')
dump('resources/ANNOTATION_PROFILE_SENSITIVITY.json',dict(created_at=now(),reason='Source/donor equal weighting can give single-cell source contributions disproportionate influence; inspect min5 target-cell donor profiles alongside all-cell balanced profiles',threshold=5,no_cell_removed_from_master=True,original_all_patient_profiles_preserved=True,role='Annotation sensitivity, not patient admission or P-value gate'))
classdefs={
 'L1':dict(name='炎症及外膜样表达候选',reference='DOCX 2.1/5.1；论文C7型对照',identity='纤维细胞表达背景；不推断骨膜或外膜组织来源'),
 'L2':dict(name='收缩及活化基质候选',reference='DOCX 2.2/5.1；论文MYH11型对照',identity='保留基质重塑及壁细胞重叠；非已验证纯myCAF'),
 'L3':dict(name='基质合成重塑及应答候选',reference='DOCX 2.3/5.1；论文FNDC1型对照',identity='区分强基质型与应答过渡型，保留原亚群来源'),
 'L4':dict(name='成骨与基质重叠表达层',reference='DOCX 2.4/5.1；论文成骨细胞与CAF区室区分',identity='成骨谱系/肿瘤竞争未闭合，不能统称已确认osteogenic CAF'),
 'U':dict(name='混合竞争谱系或功能未定保留层',reference='不强制套入四种CAF功能类别',identity='不是第五种CAF功能类型；细胞全部保留')}
ann=[
 ('J00','U','mixed_osteogenic_unresolved','low','原谱系混合且未定细胞多；四核心签名均弱，不按最高小分数强行命名。'),
 ('J01','L4','osteogenic_source_dominated','provisional','成骨核心及成骨谱系程序支持，原成骨表达细胞占多数；来源高度集中，CAF身份未闭合。'),
 ('J02','U','cycling_mixed','low','多来源增殖/混合表达，炎症、收缩、基质和成骨程序未形成清晰单一主型。'),
 ('J03','L3','response_matrix_transition','provisional','纤维细胞背景明确，mCAF与FNDC1面板中等支持；炎症核心不占优势，作为基质应答过渡群而非强基质执行群。'),
 ('J04','U','stress_secretory_unresolved','low','纤维/成骨混合且应激明显；MYH11论文面板分数受非专属成员影响，原myCAF核心及mCAF核心较弱。'),
 ('J05','L4','osteochondral_matrix_overlap','provisional','成骨核心、软骨和基质程序共同升高；原成骨与纤维谱系重叠，不命名为纯CAF。'),
 ('J06','L4','osteogenic_low_matrix','provisional','成骨核心/成骨谱系相对突出，基质重塑与收缩核心弱；保留身份未定边界。'),
 ('J07','L2','contractile_matrix_mural_overlap','supported_expression','myCAF收缩核心与mCAF基质核心同时强，MYH11/壁细胞程序叠加；按收缩活化基质归纳并保留mCAF次级特征。'),
 ('J08','L1','inflammatory_adventitial_matrix_overlap','supported_expression','iCAF原核心及C7面板在22群中突出，纤维细胞背景支持；同时有FNDC1/基质和收缩成分，不推断组织来源。'),
 ('J09','L4','osteogenic_matrix_uncertain','provisional','原成骨及未定谱系居多，成骨背景高于炎症/收缩核心；四分类中的成骨重叠研究层，身份仍待定。'),
 ('J10','U','immune_stromal_competition','low','NCI来源PTPRC高检出、MYH11和FABP4少见；不能以CD36或非专属面板分数命名MYH11型/CD36型CAF。'),
 ('J11','L4','osteogenic_matrix_overlap','provisional','原成骨谱系占多数并有成骨/基质程序；脂质面板部分信号不构成完整CD36型。'),
 ('J12','U','immune_mural_vascular_competition','low','NCI来源PTPRC、PECAM1与RGS5等竞争信号并存；FABP4少见，不确认为CD36脂质CAF或纯收缩CAF。'),
 ('J13','L4','osteogenic_mesenchymal_uncertain','provisional','成骨背景与间充质未定成分占主导，原四类中仅支持成骨重叠表达层；核心证据弱于J05等。'),
 ('J14','L4','osteogenic_stress_overlay','provisional','原成骨表达谱系高度集中，保留成骨大类及应激附加状态；原成骨CAF核心弱，不能称已证实骨化CAF。'),
 ('J15','L4','osteogenic_mineralization_overlap','provisional','成骨核心及成骨谱系程序支持，伴基质成分；保留CAF/恶性成骨竞争解释。'),
 ('J16','L4','osteogenic_matrix_overlap','provisional','基质与成骨程序重叠，原成骨/未定谱系居多；不因mCAF分数而直接归纯基质CAF。'),
 ('J17','L4','osteogenic_stress_lipid_partial','provisional','成骨核心/成骨背景突出，原成骨谱系占绝大多数；脂质面板仅部分支持，不单列CD36型。'),
 ('J18','U','cycling_stress_transition','low','跨来源增殖与应激叠加，基质/成骨/炎症部分成分并存；暂不强制定为某一CAF功能大类。'),
 ('J19','L3','matrix_synthesis_remodeling','supported_expression','mCAF原核心最强之一，FNDC1对照支持，纤维细胞背景为主；收缩及成骨背景相对低。'),
 ('J20','U','endothelial_stromal_competition','low','NCI来源PECAM1高检出且FABP4稀少，CD36可见但不足以确认为论文脂质CAF；保留血管竞争。'),
 ('J21','L3','matrix_remodeling','supported_expression','mCAF/FNDC1原面板共同支持，纤维细胞背景居多；来源和患者覆盖单独保留。')]
annotation=pd.DataFrame(ann,columns=['joint_cluster','broad_class','secondary_state','annotation_support','reason_ZH'])
for k in ['name','reference','identity']:annotation[k]=annotation.broad_class.map({c:v[k] for c,v in classdefs.items()})
assert len(annotation)==22 and set(annotation.joint_cluster)==set(m.joint_cluster)
save(annotation,'J_to_broad_class_FROZEN.tsv');dump('resources/CLASS_DEFINITIONS.json',classdefs)
dump('ANNOTATION_FROZEN.json',dict(frozen_at=now(),file=str(O/'tables/J_to_broad_class_FROZEN.tsv'),sha256=sha(O/'tables/J_to_broad_class_FROZEN.tsv'),functional_classes=4,reserve_layer='U is not a fifth CAF type',annotation_uses_inferential_P=False,manual_evidence_annotation=True,no_forced_argmax=True,no_independent_CD36_class=True,explanation='Four literature-informed functional axes plus retained unresolved/competing layer. L4 is an expression compartment with unresolved CAF identity. Statistical results will not silently change the frozen annotation.'))
event('annotation','frozen',clusters=22,classes=4,reserve_layer='U')
mp=annotation.set_index('joint_cluster');m['broad_class']=m.joint_cluster.map(mp.broad_class);m['broad_name']=m.joint_cluster.map(mp['name']);m['secondary_state']=m.joint_cluster.map(mp.secondary_state);save(m,'cell_broad_class_registry.tsv.gz')
summary=m.groupby(['broad_class','broad_name']).agg(all_cells=('cell_id','size'),J_groups=('joint_cluster','nunique'),software_positive=('doublet','sum')).reset_index();fit=m[~m.doublet].groupby('broad_class').size();summary['main_cells']=summary.broad_class.map(fit);save(summary,'broad_class_summary.tsv')
save(m[~m.doublet].groupby(['broad_class','source','donor']).size().rename('cells').reset_index(),'broad_source_patient_counts.tsv')
for oldcol in ['original_class','original_state']:
 ct=pd.crosstab(m.loc[m.cohort.eq('baseline'),oldcol],m.loc[m.cohort.eq('baseline'),'broad_class']);save(ct.reset_index(),oldcol+'_to_broad_counts.tsv')
pan_names=list(panels);raw=np.load(O/'objects/target_gene_values.npz');assert np.array_equal(raw['cell_id'],m.cell_id.to_numpy());target=raw['genes'].tolist();values=np.column_stack([scores[pan_names].to_numpy(),raw['logCP10k']]);names=pan_names+['gene__'+g for g in target];families=['signature']*len(pan_names)+['gene']*len(target)
rows=[];effects=[];scache={}
def ptest(v,key):
 n=len(v)
 if n<2 or np.all(abs(v)<1e-12):return np.nan
 if n<=16:
  if n not in scache:scache[n]=1-2*((np.arange(2**n,dtype=np.uint32)[:,None]>>np.arange(n))&1).astype(np.int8)
  return float(np.mean(abs(scache[n]@v)>=abs(v.sum())-1e-10))
 seed=int(hashlib.sha256(str(key).encode()).hexdigest()[:8],16);a=np.random.default_rng(seed).choice([-1,1],size=(65536,n));return float((np.sum(abs(a@v)>=abs(v.sum())-1e-10)+1)/65537)
for src,q in m.groupby('source'):
 for arm in ['all','singlet','singlet_min5']:
  base=q if arm=='all' else q[~q.doublet];minimum=5 if arm.endswith('min5') else 1
  for cls in sorted(base.broad_class.unique()):
   per={}
   for (donor,site),g in base.groupby(['donor','site'],dropna=False):
    a=g.index[g.broad_class.eq(cls)].to_numpy();b=g.index[g.broad_class.ne(cls)].to_numpy()
    if min(len(a),len(b))<minimum:continue
    # A missing endpoint in either complete comparison arm remains unassessed for that patient/site.
    v=values[a].mean(0)-values[b].mean(0);per.setdefault(donor,[]).append((v,len(a),len(b)))
   donors=sorted(per);e=np.array([np.mean([t[0] for t in per[d]],axis=0) for d in donors]) if donors else np.empty((0,len(names)))
   for j,name in enumerate(names):
    ok=np.isfinite(e[:,j]);v=e[ok,j];ds=np.array(donors)[ok];n=len(v);p=ptest(v,(src,arm,cls,name))
    for donor,delta in zip(ds,v):effects.append(dict(source=src,arm=arm,broad_class=cls,endpoint=name,donor=donor,effect=delta))
    rows.append(dict(source=src,arm=arm,broad_class=cls,endpoint=name,family=families[j],cells=int(base.broad_class.eq(cls).sum()),eligible_pairs_before_coverage=len(donors),N=n,mean_effect=float(v.mean()) if n else np.nan,positive_patients=int((v>1e-12).sum()),P=p,test='insufficient_N' if n<2 else 'all_zero_not_tested' if np.all(abs(v)<1e-12) else 'exact_signflip' if n<=16 else 'MC65536_plus_observed'))
 event('statistics','source_completed',source=src)
df=pd.DataFrame(rows);df['FDR']=np.nan
for _,idx in df.groupby(['source','arm','family']).groups.items():
 ii=df.loc[idx].index[df.loc[idx].P.notna()];df.loc[ii,'FDR']=false_discovery_control(df.loc[ii,'P'].to_numpy())
df['nominal']=df.P.lt(.05);df['FDR_significant']=df.FDR.lt(.05);save(df,'broad_class_statistics.tsv');save(pd.DataFrame(effects),'broad_patient_effects.tsv.gz')
# Primary and sensitivity broad profiles, preserving all marker and competing-lineage axes.
for mr in [1000,1500,2000]:
 f=read('cell_rank_scores_'+str(mr)+'.tsv.gz');f['broad_class']=m.broad_class;pp=f[~f.doublet].groupby(['source','donor','broad_class'])[pan_names].mean().reset_index();save(pp,'broad_patient_panel_scores_'+str(mr)+'.tsv');sp=pp.groupby(['source','broad_class'])[pan_names].mean().reset_index();save(sp,'broad_source_panel_scores_'+str(mr)+'.tsv');save(sp.groupby('broad_class')[pan_names].mean().reset_index(),'broad_balanced_panel_scores_'+str(mr)+'.tsv')
robust=[];main=read('balanced_cluster_scores_1500.tsv').set_index('joint_cluster')
for mr in [1000,2000]:
 q=read('balanced_cluster_scores_'+str(mr)+'.tsv').set_index('joint_cluster').loc[main.index]
 for pn in pan_names:robust.append(dict(maxRank=mr,panel=pn,cluster_profile_spearman=float(spearmanr(main[pn],q[pn],nan_policy='omit').statistic)))
save(pd.DataFrame(robust),'maxRank_sensitivity.tsv')
# Co-detection is descriptive RNA evidence, not flow/protein positivity.
gi={g:i for i,g in enumerate(target)};rv=raw['raw'];cod=[]
for (src,donor,cl),g in m[~m.doublet].groupby(['source','donor','joint_cluster']):
 ix=g.index.to_numpy();co={'CD36_FABP4':['CD36','FABP4'],'CD36_PTPRC':['CD36','PTPRC'],'CD36_PECAM1':['CD36','PECAM1'],'COL1A1_PTPRC':['COL1A1','PTPRC']}
 for n,gs in co.items():
  vv=rv[ix][:,[gi[x] for x in gs]];ok=np.isfinite(vv).all(1);cod.append(dict(source=src,donor=donor,joint_cluster=cl,combination=n,measured_cells=int(ok.sum()),fraction=float((vv[ok]>0).all(1).mean()) if ok.any() else np.nan))
cod=pd.DataFrame(cod);save(cod,'patient_RNA_codetection.tsv');save(cod.groupby(['source','joint_cluster','combination']).fraction.mean().reset_index(),'source_RNA_codetection.tsv')
old=pd.read_csv(R/'phase66_doublet_repair_DSC1_20261005/tables/cells_repaired_with_confidence.tsv.gz',sep='\t',low_memory=False).set_index('cell_id');nci=pd.read_csv(R/'phase102_NCI_author_CAF_20261020/tables/author_CAF_cell_registry.tsv.gz',sep='\t').set_index('cell_id');ncc=pd.read_csv(R/'phase93_expanded_admission_20261016/tables/expanded_registry_8980.tsv.gz',sep='\t',low_memory=False);ncc=ncc[ncc.registry_origin.eq('NCC_addition_806')].set_index('cell_id')
m['CNV_record']=m.cell_id.map(old.CNV_evidence_current);m.loc[m.source.eq('NCI'),'CNV_record']=m.loc[m.source.eq('NCI'),'cell_id'].map(nci.native_copykat_call);m.loc[m.source.eq('NCC'),'CNV_record']=m.loc[m.source.eq('NCC'),'cell_id'].map(ncc.native_copykat_call)
save(m[~m.doublet].groupby(['broad_class','source','CNV_record'],dropna=False).size().rename('cells').reset_index(),'inherited_CNV_evidence_by_broad_class.tsv')
assert sha(O/'tables/J_to_broad_class_FROZEN.tsv')==json.loads((O/'ANNOTATION_FROZEN.json').read_text())['sha256']
dump('ANALYSIS_COMPLETE.json',dict(status='PASS_COMPUTATION',all_cells=len(m),main_cells=int((~m.doublet).sum()),J_groups=22,functional_classes=4,reserve_J_groups=int(annotation.broad_class.eq('U').sum()),test_rows=len(df),valid_tests=int(df.P.notna().sum()),annotation_frozen_before_P=True,annotation_not_changed_by_P=True,independent_CD36_class=False,old_F_J_labels_changed=False))
event('analysis','completed',summary=summary.to_dict('records'));print(summary.to_string(index=False))
