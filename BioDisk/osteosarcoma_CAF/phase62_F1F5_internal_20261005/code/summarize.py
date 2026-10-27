from pathlib import Path
import json,re
import numpy as np,pandas as pd
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase62_F1F5_internal_20261005')
def read(p,**kw):return pd.read_csv(p,sep='\t',low_memory=False,**kw)
def dump(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str),encoding='utf-8')
assert (O/'PREPARE_COMPLETE').exists()
roster=read(O/'tables/subgroup_roster.tsv');ob=read(O/'tables/cells_original_with_review.tsv.gz').set_index('cell_id');top=read(O/'tables/marker_display.tsv');cr=read(O/'tables/contrast_registry.tsv');rows=[];supports=[];markerchecks=[];brief=[]
def technical(g):return g.startswith(('MT-','RPL','RPS','MTCO','MTND','MTATP')) or g in ['MALAT1','NEAT1']
for rec in roster.to_dict('records'):
 st=rec['subgroup_id'];ds=rec['primary_source'];cell=ob[ob.subgroup_id.eq(st)];ctx=f'{ds}__{st}__corrected__all';t=top[top.context.eq(ctx)].copy();cor=read(O/'tables/genes'/f'{ctx}.tsv.gz').set_index('symbol');raw=read(O/'tables/genes'/f'{ds}__{st}__raw__all.tsv.gz').set_index('symbol');genes=t.symbol.tolist();preserved=float((raw.loc[genes,'mean_delta']>0).mean()) if genes else 0
 stable=[g for g in genes if cor.loc[g,'mean_delta']>=.5 and cor.loc[g,'detection_difference']>=.1 and raw.loc[g,'mean_delta']>0 and not technical(g)];expr='高' if len(stable)>=3 and preserved>=.8 else '中' if stable else '低';fractions=cell.lineage_stratum.value_counts(normalize=True);flags=[]
 if rec['max_donor_fraction']>=.9:flags.append('单患者主导')
 if cell.removed_fraction.median()>=.5:flags.append('去污染影响大')
 if cell.mitochondrial_fraction.median()>=.2:flags.append('线粒体比例偏高')
 for g in genes:markerchecks.append({'subgroup_id':st,'symbol':g,'primary_source':ds,'corrected_effect':cor.loc[g,'mean_delta'],'raw_effect':raw.loc[g,'mean_delta'],'corrected_detection':cor.loc[g,'target_detection'],'reference_detection':cor.loc[g,'reference_detection'],'stable_display_marker':g in stable,'technical_marker':technical(g)})
 assess=stable if stable else [g for g in genes if not technical(g)];best=None
 for path in sorted((O/'objects').glob(f'*__{st}__paired.npz')):
  z=np.load(path);gi={g:i for i,g in enumerate(z['symbols'])};available=[g for g in assess if g in gi];d=z['effects'][:,[gi[g] for g in available]];n=len(d);concord=float((d>0).mean()) if n and available else np.nan;source=path.name.split('__')[0];grade='高' if n>=10 and concord>=.7 and len(available)>=3 else '中' if n>=3 and concord>=.6 and len(available)>=1 else '有限';r={'subgroup_id':st,'dataset':source,'paired_donors':n,'markers_measured':len(available),'marker_direction_concordance':concord,'cross_patient_support':grade,'selection_context':'posthoc_same_data_marker_direction_support_not_independent_validation'};supports.append(r)
  if best is None or ({'高':2,'中':1,'有限':0}[grade],n,0 if np.isnan(concord) else concord)>({'高':2,'中':1,'有限':0}[best['cross_patient_support']],best['paired_donors'],0 if np.isnan(best['marker_direction_concordance']) else best['marker_direction_concordance']):best=r
 rep=best['cross_patient_support'];overall='高' if expr=='高' and rep=='高' else '中' if expr in ['高','中'] and rep in ['高','中'] else '探索性';fibctx=f'{ds}__{st}__corrected__Fibro';fp=O/'tables/genes'/f'{fibctx}.tsv.gz';fibpres=np.nan
 if fp.exists() and assess:
  f=read(fp).set_index('symbol');fibpres=float((f.loc[assess,'mean_delta']>0).mean())
  if fibpres<.5:flags.append('Fibro限定后标志改变')
 rec.update(expression_confidence=expr,cross_patient_support=rep,overall_evidence=overall,stable_marker_count=len(stable),raw_direction_retention=preserved,stable_markers=';'.join(stable),display_markers=';'.join(genes),Fibro_marker_direction_retention=fibpres,support_source=best['dataset'],support_paired_donors=best['paired_donors'],support_direction_fraction=best['marker_direction_concordance'],dominant_lineage=fractions.index[0],dominant_lineage_fraction=fractions.iloc[0],median_removed_fraction=cell.removed_fraction.median(),flags=';'.join(flags),confidence_reason=f'{len(stable)} stable descriptive markers; raw sign retention {preserved:.1%}; strongest source {best["dataset"]}, paired N={best["paired_donors"]}, marker concordance={best["marker_direction_concordance"]}')
 rows.append(rec);brief.append(st+' | '+str(rec['cells'])+' | '+expr+'/'+rep+' | '+','.join(stable)+' | '+rec['flags'])
result=pd.DataFrame(rows);result.to_csv(O/'tables/subgroup_confidence.tsv',sep='\t',index=False);pd.DataFrame(supports).to_csv(O/'tables/cross_patient_support.tsv',sep='\t',index=False);pd.DataFrame(markerchecks).to_csv(O/'tables/marker_raw_corrected_evidence.tsv',sep='\t',index=False)
if (O/'GSEA_COMPLETE').exists():
 g=read(O/'tables/all_GSEA_results.tsv.gz');show=[]
 for st in result.subgroup_id:
  for scope in ['all','Fibro']:
   z=g[g.contrast.eq(st)&g.scope.eq(scope)&g.NES.gt(0)&g.GSEA_model_P.lt(.05)]
   for lib in ['hallmark','reactome','matrisome']:
    v=z[z.library.eq(lib)].sort_values(['model_BH_FDR','NES'],ascending=[True,False]).head(5);show.append(v)
 pd.concat(show).to_csv(O/'tables/pathway_display.tsv',sep='\t',index=False)
 for st in result.subgroup_id:
  z=g[g.contrast.eq(st)&g.scope.eq('all')&g.NES.gt(0)&g.GSEA_model_P.lt(.05)]
  brief.append('\n'+st+' PATHWAYS')
  for lib in ['hallmark','matrisome','reactome']:
   v=z[z.library.eq(lib)].sort_values(['model_BH_FDR','NES'],ascending=[True,False]).head(4)
   brief.extend(f'{r.Term}: NES={r.NES:.2f} P={r.GSEA_model_P:.4g} BH={r.model_BH_FDR:.4g}' for r in v.itertuples())
(O/'QA/INTERPRETATION_BRIEF.txt').write_text('\n'.join(brief),encoding='utf-8');print('\n'.join(brief),flush=True)
