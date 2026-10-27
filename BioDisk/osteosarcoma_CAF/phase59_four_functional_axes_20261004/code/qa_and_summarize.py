from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
import pandas as pd,numpy as np,json,subprocess,sys,hashlib
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase59_four_functional_axes_20261004';B=P/'phase58_F1F5_round1_20261004'
def read(p,**kw):return pd.read_csv(p,sep='\t',low_memory=False,**kw)
assert (O/'analysis.exit_code').read_text().strip()=='0' and (O/'ANALYSIS_COMPLETE').exists()
s=read(O/'tables/axis_statistics.tsv');g=read(O/'tables/marker_statistics.tsv.gz');e=read(O/'tables/patient_axis_effects.tsv.gz');ob=read(O/'tables/frozen_cells_with_identity.tsv.gz').set_index('cell_id');genes=json.loads((O/'resources/GENE_SETS.json').read_text());sets=genes['definitions'];zf=np.load(O/'objects/selected_gene_counts.npz');z={k:zf[k] for k in zf.files};zf.close();gi={g:i for i,g in enumerate(z['genes'])}
for n in sorted(set(s.N)|set(g.N)):
 if n>18:np.savetxt(O/f'QA/MC_signs_N{n}.tsv.gz',np.random.default_rng(20261004).choice([-1.,1.],size=(99999,n)),fmt='%d',delimiter='\t')
# Independent pseudobulk reconstruction for the four axes in main source, all donor/site pairs and both matrices.
trace=json.loads((O/'QA/main_pseudobulk_reconstruction_indices.json').read_text());checks=[];reconstructed=[]
for t in trace:
 mode=t['matrix'];il=np.array(t['left_indices']);ir=np.array(t['right_indices']);cl=t['contrast'].split('_')[0];assert ob.iloc[il].class_id.eq(cl).all() and ob.iloc[ir].class_id.ne(cl).all();assert not set(il)&set(ir)
 for axis in ['inflammation','contraction','matrix','osteogenic']:
  ix=np.array([gi[g] for g in sets[axis]]);valid=z['assayed'][np.r_[il,ir]][:,ix].all(0);use=ix[valid]
  if len(use)<max(min(3,len(ix)),int(np.ceil(.8*len(ix)))):value=np.nan
  else:
   lc=z[mode][np.ix_(il,use)].sum(axis=0);rc=z[mode][np.ix_(ir,use)].sum(axis=0);ld=z[mode+'_library'][il].sum();rd=z[mode+'_library'][ir].sum();value=float(np.mean(np.log2((lc/ld*10000+.1)/(rc/rd*10000+.1))))
  reconstructed.append({'matrix':mode,'contrast':t['contrast'],'program':axis,'donor_record':t['donor_record'],'effect':value})
r=pd.DataFrame(reconstructed).groupby(['matrix','contrast','program','donor_record'],as_index=False).effect.mean();exp=e[(e.dataset=='SCPCP000017')&(e.context=='main')&e.program.isin(['inflammation','contraction','matrix','osteogenic'])];cmp=r.merge(exp,on=['matrix','contrast','program','donor_record'],suffixes=('_reconstructed','_saved'),validate='one_to_one');assert len(cmp)==len(exp);assert np.allclose(cmp.effect_reconstructed,cmp.effect_saved,atol=1e-10,equal_nan=True);cmp.to_csv(O/'QA/main_patient_axis_reconstruction.tsv',sep='\t',index=False)
assert ob.class_id.value_counts().to_dict()=={'F4':5304,'F5':1699,'F1':695,'F3':251,'F2':225}
con=read(O/'tables/F2F4_lineage_RNA_contribution.tsv.gz');key=['dataset','class_id','donor_record','site','matrix','gene'];mass=con.groupby(key).agg(component_sum=('RNA_counts','sum'),total=('total_class_RNA_counts','first'),fraction_sum=('RNA_fraction','sum'));ev=mass.total.notna();assert np.allclose(mass.loc[ev,'component_sum'],mass.loc[ev,'total'],atol=1e-7);assert np.allclose(mass.loc[mass.total.gt(0),'fraction_sum'],1)
comp=read(O/'tables/patient_class_lineage_composition.tsv');assert np.allclose(comp.groupby(['dataset','class_id','donor_record','site']).cell_fraction.sum(),1)
if not (O/'QA/INDEPENDENT_R_COMPLETE').exists():
 with (O/'QA/independent_R.log').open('w') as f:q=subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript','--vanilla',str(O/'code/independent_QA.R')],stdout=f,stderr=subprocess.STDOUT)
 (O/'QA/independent_R.exit_code').write_text(str(q.returncode));assert q.returncode==0,(O/'QA/independent_R.log').read_text()[-4000:]
else:assert (O/'QA/independent_R.exit_code').read_text().strip()=='0'
qa={'original_cells_and_labels_preserved':True,'main_patient_axis_reconstructions':len(cmp),'max_reconstruction_error':float((cmp.effect_reconstructed-cmp.effect_saved).abs().max()),'RNA_mass_and_fraction_conservation':True,'R_QA':read(O/'QA/independent_R_summary.tsv').to_dict('records')};(O/'delivery/QA_ACCEPTANCE.json').write_text(json.dumps(qa,indent=2))
# Summaries keep sources distinct. Patient equal contributions average sites first, absent lineage kept as zero.
patcon=con.groupby(['dataset','class_id','donor_record','matrix','gene','lineage'],as_index=False).RNA_fraction.mean();cs=patcon.groupby(['dataset','class_id','matrix','gene','lineage'],as_index=False).agg(mean_patient_RNA_fraction=('RNA_fraction','mean'),N_evaluable=('RNA_fraction','count'));cs.to_csv(O/'tables/lineage_RNA_patient_equal_summary.tsv',sep='\t',index=False)
pool=con.groupby(['dataset','class_id','matrix','gene','lineage'],as_index=False).RNA_counts.sum();pool['pooled_RNA_fraction']=pool.RNA_counts/pool.groupby(['dataset','class_id','matrix','gene']).RNA_counts.transform('sum');pool.to_csv(O/'tables/lineage_RNA_pooled_summary.tsv',sep='\t',index=False)
four=['inflammation','contraction','matrix','osteogenic'];v=e[(e.context=='main')&e.program.isin(four)].rename(columns={'matrix':'count_matrix'}).pivot(index=['dataset','count_matrix','contrast','donor_record'],columns='program',values='effect');v.to_csv(O/'tables/patient_four_axis_vectors.tsv',sep='\t');v['n_positive_axes']=(v[four]>0).sum(1);v['all_four_evaluable']=v[four].notna().all(1);v.reset_index().to_csv(O/'tables/patient_four_axis_overlap.tsv',sep='\t',index=False)
paired=s[s.context.eq('main')&s.program.isin(four)].pivot(index=['dataset','contrast','program'],columns='matrix',values=['mean_delta','P','FDR','N']);paired.to_csv(O/'tables/main_raw_corrected_comparison.tsv',sep='\t')
print('QA_PASS',json.dumps(qa),flush=True)
print('MAIN_SOURCE',s[(s.dataset=='SCPCP000017')&(s.context=='main')&s.program.isin(four)][['matrix','contrast','program','N','mean_delta','positive_patients','P','FDR']].to_string(index=False),flush=True)
sel=s[(s.dataset=='SCPCP000017')&(s.matrix=='corrected')&(s.context.isin(['Fibro_only','conservative_Fibro','within_F2','within_F4','min10']))&s.program.isin(four+['contractile_markers','matrix_markers','osteogenic_effectors','Liu_C2_OMD_OGN','Liu_C5_ossification','mineralization'])];sel.to_csv(O/'delivery/FOCUS_RESULTS.tsv',sep='\t',index=False)
print('FOCUS',sel[sel.contrast.str.contains('F2|F4')|sel.context.str.startswith('within')][['context','contrast','program','N','mean_delta','P','FDR']].to_string(index=False),flush=True)
print('EXTERNAL',s[(s.dataset=='SCPCP000017')&s.context.str.startswith('external')&s.program.isin(four)][['contrast','program','N','mean_delta','P','FDR']].to_string(index=False),flush=True)
print('CONTRIBUTIONS',cs[(cs.dataset=='SCPCP000017')&(cs.matrix=='corrected')&cs.gene.isin(['ACTA2','MYL9','COL1A1','POSTN','RUNX2','ALPL','SP7','IBSP','OMD','OGN'])].to_string(index=False),flush=True)
(O/'QA_COMPLETE').write_text('SUCCESS\n')
