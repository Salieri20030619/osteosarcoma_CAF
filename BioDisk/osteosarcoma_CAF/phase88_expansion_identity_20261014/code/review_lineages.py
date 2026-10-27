"""SYNTAX_CHECK: python -m py_compile code/review_lineages.py
RUN_COMMAND: python code/review_lineages.py
EXPECTED_SUCCESS_MARKER: GEO_LINEAGE_REVIEW_READY
"""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='4'
from pathlib import Path
import ast,json,gzip,hashlib
import numpy as np,pandas as pd,anndata as ad
from threadpoolctl import threadpool_limits
threadpool_limits(limits=4)
r=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase88_expansion_identity_20261014/preprocess');out=r/'review';out.mkdir(exist_ok=True)
src=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase3_identity_20260914/code/review_inputs.py')
for n in ast.parse(src.read_text()).body:
 if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='PANELS' for t in n.targets):PANELS=ast.literal_eval(n.value)
PANELS.update({'Mast':['TPSAB1','TPSB2','CPA3','MS4A2','KIT'],'Erythroid':['HBA1','HBA2','HBB','AHSP','ALAS2','SLC4A1','KLF1'],'Myogenic':['MYOD1','MYOG','MYMX','ACTA1','MYH3','MYH8','MYLPF','TNNT1','TNNT3','KLHL41','CASQ2','MYBPC1'],'Stress_core':['FOS','FOSB','JUN','JUNB','ATF3','NR4A1','HSPA1A','HSPA1B','DDIT3']})
ANCHORS={k:list(v) for k,v in PANELS.items() if k in ['T_NK','B_Plasma','Myeloid','DC','Endo','Epi','Fib','Osteo','Mural','Mast','Erythroid','Myogenic']}
ANCHORS['Myeloid']+=['C1QA','C1QB','C1QC','ITGAM','ITGB2','MPEG1','FCGR3A','CD14','TREM1','MARCO','AIF1','LST1']
ANCHORS['Endo']+=['ADGRL4','TIE1','CLEC14A','ROBO4','ESAM']
ANCHORS['Mural']+=['ACTG2','LMOD1','HIGD1B','HEY1']
ANCHORS['Osteo']+=['SOST','PTH1R'];ANCHORS['T_NK']+=['CD3G','TRAC','TRBC1','TRBC2','LCK','IL7R']
# Include the extra anchors in quantitative evidence, not only marker-list matching.
for k,genes in ANCHORS.items():PANELS[k]=sorted(set(PANELS[k]+genes))
PANELS['OC']+=['SIGLEC15','SUCNR1','NFATC1']
manifest=pd.read_csv(r/'manifest.tsv',sep='\t').fillna('')
assert len(manifest)==12
packets=[];summaries=[];coverage=[];ambient=[]
for row in manifest.to_dict('records'):
 sample=row['sample'];p=r/'results'/sample;assert (p/'complete.json').exists(),sample
 a=ad.read_h5ad(p/'annotated.h5ad');x=a.layers['counts'];logx=a.X;o=a.obs.copy();o['cluster']=o.cluster.astype(str)
 assert len(o)==json.loads((p/'complete.json').read_text())['retained']
 assert np.allclose(np.asarray(x.sum(1)).ravel(),o.n_counts)
 symbols=a.var.symbol.astype(str).to_numpy();geneix={g:np.flatnonzero(symbols==g) for g in set(sum(PANELS.values(),[]))}
 for name,genes in PANELS.items():
  detected=np.zeros(a.n_obs,dtype=np.int16);expr=np.zeros(a.n_obs);umis=np.zeros(a.n_obs)
  for g in genes:
   ix=geneix[g]
   if len(ix):
    count=np.asarray(x[:,ix].sum(1)).ravel();detected+=(count>0);umis+=count;expr+=np.asarray(logx[:,ix].sum(1)).ravel()
  o['det_'+name+'_review']=detected;o['meanlog_'+name]=expr/max(1,len(genes));o['pctUMI_'+name]=umis/o.n_counts*100
  coverage.append({'sample':sample,'panel':name,'requested':len(genes),'covered':sum(len(geneix[g])>0 for g in genes),'missing':';'.join(g for g in genes if not len(geneix[g]))})
 mk=pd.read_csv(p/'cluster_markers.tsv',sep='\t',dtype={'group':str});labels={};reasons={}
 for cluster,cc in o.groupby('cluster',sort=False):
  top=mk.loc[mk.group==cluster,'symbol'].head(50).astype(str).tolist();tset=set(top)
  hits={k:sorted(set(v)&tset) for k,v in ANCHORS.items()}
  active=[k for k in ANCHORS if len(hits[k])>=2 and (cc['det_'+k+'_review']>=2).mean()>=.2]
  oc_special=len(set(['OCSTAMP','DCSTAMP','RUFY4','ATP6V0D2','SIGLEC15','SUCNR1','CALCR'])&tset)
  oc_core=len(set(['ACP5','CTSK'])&tset)
  oc=oc_special>=1 and oc_core>=1 and (cc.det_OC_review>=2).mean()>=.2
  immune=[k for k in active if k in ['T_NK','B_Plasma','Myeloid','DC','Mast','Erythroid']]
  structural=[k for k in active if k in ['Endo','Epi','Mural','Fib','Osteo','Myogenic']]
  reason=[];lab='Unresolved'
  if oc:lab='OC_like';reason.append('osteoclast_core_and_specialized_markers')
  elif immune and not structural:
   # DC and myeloid overlap is retained as a myeloid-family ambiguity when tied.
   scores={k:len(hits[k])*float((cc['det_'+k+'_review']>=2).mean()) for k in immune}
   order=sorted(scores,key=scores.get,reverse=True)
   lab=order[0] if len(order)==1 or scores[order[0]]>=1.5*scores[order[1]] else 'Immune_ambiguous'
  elif immune and structural:lab='Mixed'
  elif len(structural)==1:lab={'Fib':'FibLike','Osteo':'OsteoLike','Myogenic':'Myogenic_pending'}.get(structural[0],structural[0])
  elif len(structural)>1:
   lab='MesAmb' if set(structural)<=set(['Mural','Fib','Osteo']) else 'Mixed'
  else:
   # Retain prior quantitative lineage evidence as provisional, never as normal reference proof.
   lab='Unresolved'
  cta=len(tset&set(PANELS['CTA']))>=2
  cycle=len(tset&set(PANELS['Cycle']))>=2
  if cta:lab='MesAmb';reason.append('cancer_testis_program_identity_pending')
  if cycle:reason.append('cycling_program')
  reason.append('active='+','.join(active))
  labels[cluster]=lab;reasons[cluster]=';'.join(reason)
  z={'sample':sample,'short':row['short'],'cluster':cluster,'cells':len(cc),'original_label':str(cc.lineage_prelim.iloc[0]),'proposed_lineage':lab,'CTA_context':cta,'cycle_context':cycle,'reason':reasons[cluster],'top15':','.join(top[:15])}
  for k in ANCHORS:z[k+'_hits']=','.join(hits[k]);z[k+'_pct2']=round(float((cc['det_'+k+'_review']>=2).mean()*100),2)
  packets.append(z)
 o['lineage_review_proposed']=o.cluster.map(labels);o['lineage_review_reason']=o.cluster.map(reasons)
 o['phase2_label']=o.lineage_prelim.astype(str);o['lineage_phase3']=o.lineage_review_proposed
 mes=(o.det_Fib_review>=2)|(o.det_Osteo_review>=2)|(o.det_Mural_review>=2)
 immune_strong=(o.det_T_NK_review>=3)|(o.det_B_Plasma_review>=3)|(o.det_Myeloid_review>=4)
 mes_strong=(o.det_Fib_review>=4)|(o.det_Osteo_review>=3)|(o.det_Mural_review>=4)
 o['mixed_lineage_flag']=(immune_strong&mes_strong)|((o.det_Endo_review>=4)&(o.det_Osteo_review>=3))
 qscore=o.doublet_score.quantile(.9)
 scorepass=o.doublet_score.notna()&(o.doublet_score<=qscore)
 base=(o.n_genes>=500)&(o.n_counts>=750)&scorepass&~o.mixed_lineage_flag
 o['ref_eligible_immune']=False
 for k in ['T_NK','B_Plasma','Myeloid','DC']:
  o['ref_eligible_immune']|=(o.lineage_review_proposed==k)&(o['det_'+k+'_review']>=2)&base&~mes&(o.det_OC_review<2)&(o.det_Epi_review<2)&(o.det_Endo_review<2)
 o['ref_eligible_endo']=(o.lineage_review_proposed=='Endo')&(o.det_Endo_review>=3)&base&~mes&(o.det_OC_review<2)&(o.det_Epi_review<2)
 o['cnv_target_expanded']=o.lineage_review_proposed.isin(['FibLike','OsteoLike','Mural','MesAmb','Mixed','Unresolved','Myogenic_pending','OC_like'])|mes
 o['cell_id']=o.index;o['short']=row['short'];o['sample']=sample
 o.to_csv(out/(sample+'_cell_evidence.tsv.gz'),sep='\t',index=False)
 summary={'sample':sample,'retained':len(o),'proposed_fibroblast_cells':int((o.lineage_review_proposed=='FibLike').sum()),'broad_cnv_targets':int(o.cnv_target_expanded.sum()),'eligible_immune_reference':int(o.ref_eligible_immune.sum()),'eligible_endothelial_reference':int(o.ref_eligible_endo.sum()),'mixed_lineage_flags':int(o.mixed_lineage_flag.sum()),'status':'marker_review_pending_CNV'}
 summaries.append(summary)
 for group in ['T_NK','B_Plasma','Myeloid','Endo']:
  mask=(o.lineage_review_proposed==group).values
  if mask.sum():
   for gene in ['COL1A1','SPP1','IBSP']:
    ix=np.flatnonzero(symbols==gene);vals=np.asarray(x[:,ix].sum(1)).ravel()
    ambient.append({'sample':sample,'lineage':group,'cells':int(mask.sum()),'gene':gene,'pct_detected':float((vals[mask]>0).mean()*100),'median_UMI':float(np.median(vals[mask])),'interpretation':'cross_lineage_expression_screen; not a contamination estimate or correction'})
 print('REVIEW_PACKET_READY',json.dumps(summary),flush=True)
 del a,x,logx,o
pd.DataFrame(packets).to_csv(out/'cluster_review_packets.tsv',sep='\t',index=False)
pd.DataFrame(summaries).to_csv(out/'cnv_readiness.tsv',sep='\t',index=False)
pd.DataFrame(coverage).to_csv(out/'panel_coverage.tsv',sep='\t',index=False)
pd.DataFrame(ambient).to_csv(out/'cross_lineage_expression_screen.tsv',sep='\t',index=False)
(out/'lineage_rules.json').write_text(json.dumps({'panels':PANELS,'anchors':ANCHORS,'note':'Candidate proposals for evidence review, no CAF confirmation before CNV'},indent=2))
(out/'PACKETS_COMPLETE').write_text('GEO_LINEAGE_REVIEW_READY\n')
print('GEO_LINEAGE_REVIEW_READY',len(summaries),len(packets),flush=True)
