from pathlib import Path
import os,select,subprocess,json
import pandas as pd,numpy as np
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');B=P/'phase56_patient_programs_20261004';O=P/'phase57_downstream_update_20261004';Q=O/'independent_QA';Q.mkdir(exist_ok=True)
if not (O/'downstream.exit_code').exists():
 try:
  fd=os.pidfd_open(int((O/'downstream.pid').read_text()));select.select([fd],[],[]);os.close(fd)
 except ProcessLookupError:pass
assert (O/'downstream.exit_code').read_text().strip()=='0' and (O/'DOWNSTREAM_COMPLETE').is_file()
stat=pd.read_csv(B/'tables/all_patient_program_statistics.tsv.gz',sep='\t');contexts=pd.read_csv(B/'tables/executed_contexts.tsv',sep='\t')
for n in sorted(set(contexts.N)):
 if n>15:np.savetxt(Q/f'MC_signs_N{n}.tsv.gz',np.random.default_rng(20261004).choice([-1.,1.],size=(19999,n)),fmt='%d',delimiter='\t')
selected=[];effects=[]
for f in sorted((O/'tables/genes').glob('*.tsv.gz')):
 ctx=f.name.replace('.tsv.gz','');s=pd.read_csv(f,sep='\t');rng=np.random.default_rng(20261004);fixed=['COL1A1','LOX','MMP14','RUNX2','ANGPT1','SMAD3'];fixed+=rng.choice(s.gene.to_numpy(),size=min(12,len(s)),replace=False).tolist();t=s[s.gene.isin(fixed)].copy();t['context_id']=ctx;selected.append(t)
 z=np.load(B/'objects'/f'{ctx}_gene_effects.npz');ix={x:i for i,x in enumerate(z['genes'].astype(str))}
 for gene in t.gene:
  for i,donor in enumerate(z['donors'].astype(str)):effects.append({'context_id':ctx,'donor_record':donor,'gene':gene,'effect':z['effects'][i,ix[gene]]})
pd.concat(selected).to_csv(Q/'gene_P_check_statistics.tsv',sep='\t',index=False);pd.DataFrame(effects).to_csv(Q/'gene_P_check_effects.tsv',sep='\t',index=False)
shares=pd.read_csv(B/'tables/patient_site_target_gene_RNA_shares.tsv',sep='\t');assert shares.candidate_RNA_share.dropna().between(-1e-12,1+1e-12).all()
for _,g in stat.groupby('family'):
 assert len(g)==84 and g.program.is_unique
 assert (~(g.N.lt(2)&g.P.notna())).all()
assert stat.N.max()<=26 and stat.P.dropna().between(0,1).all() and stat.FDR.dropna().between(0,1).all()
env=os.environ.copy();env.update(OMP_NUM_THREADS='4',OPENBLAS_NUM_THREADS='4',MKL_NUM_THREADS='4')
with (Q/'independent_R_QA.log').open('w') as f:r=subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript','--vanilla',str(O/'code/independent_QA.R')],stdout=f,stderr=subprocess.STDOUT,env=env)
(Q/'independent_R_QA.exit_code').write_text(str(r.returncode));assert r.returncode==0
(Q/'PYTHON_OUTPUT_QA.json').write_text(json.dumps({'program_family_slots':84,'source_scoped_N_max':int(stat.N.max()),'gene_P_sampling':'Six prior target genes plus 12 reproducibly random genes per evaluable raw background context','all_program_P_checked':True,'all_gene_program_ORA_FDR_checked':True,'RNA_shares_bounded':True,'no_N_below_2_P_values':True},indent=2))
print('INDEPENDENT_STATISTICAL_QA_COMPLETE',flush=True)