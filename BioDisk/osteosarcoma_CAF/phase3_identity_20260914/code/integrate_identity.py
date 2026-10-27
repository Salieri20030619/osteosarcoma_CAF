"""Integrate two actual CNV results, lineage and QC; retain discordant/unevaluable cells.
SYNTAX_CHECK: python -m py_compile integrate_identity.py
RUN_COMMAND: python integrate_identity.py --root <phase2> --phase3 <phase3>
EXPECTED_SUCCESS_MARKER: IDENTITY_INTEGRATION_COMPLETE
"""
import argparse,json,hashlib
from pathlib import Path
import pandas as pd,numpy as np,h5py
import scipy.sparse as sp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import anndata as ad
from takeover_audit import frame,sparse
SHORTS=['S0058','S0059','S0217','S0218','SC069','SC072','SC073']
MARKERS=['DCN','LUM','PDGFRA','FAP','COL1A1','COL11A1','ACTA2','TAGLN','MYH11','RGS5','RUNX2','ALPL','IBSP','MGP','ENPP1','SPP1','CD36','PECAM1','PTPRC']
IDS={'Fib_candidate':'FibCand','Malignant_OS':'MalOS','Mural_SMC':'Mural','Other':'Other','Identity_pending':'Unk'}
COL={'FibCand':'#196B31','MalOS':'#A51426','Mural':'#92278F','Other':'#B0B0B0','Unk':'#4B5590'}
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--phase3',type=Path,required=True);a=ap.parse_args();p=a.phase3
    for s in SHORTS:
        for method in ['infercnv','copykat']:
            if not (p/method/s/'COMPLETE').exists():raise RuntimeError(f'CNV not complete: {method} {s}')
    out=p/'identity';out.mkdir(exist_ok=True);(out/'figures').mkdir(exist_ok=True);(out/'objects').mkdir(exist_ok=True)
    final=[];toolqc=[];clstats=[];genesummary=[];strict_objs=[];extended_objs=[];summaries=[];selections=[]
    for short in SHORTS:
        o=pd.read_csv(p/'review'/f'{short}_cell_evidence.tsv.gz',sep='\t').set_index('cell_id',drop=False)
        gsm=o['sample'].iloc[0]
        marker_table=pd.read_csv(a.root/'results'/gsm/'cluster_markers.tsv',sep='\t')
        cta_genes={'MAGEA3','MAGEA10','MAGEA12','SSX1','CTAG2','PAGE2B','PAGE5','CSAG1','FMR1NB'}
        cta_clusters={int(c):';'.join(sorted(cta_genes&set(g.symbol))) for c,g in marker_table.groupby('group') if len(cta_genes&set(g.symbol))>=2}
        o['cluster_CTA_program']=o.cluster.astype(int).isin(cta_clusters)
        o['cluster_CTA_marker_genes']=o.cluster.astype(int).map(cta_clusters).fillna('')
        qualities={}
        for method,normal in [('infercnv','low_signal'),('copykat','diploid'),('copykat_balanced','diploid'),('scevan','normal')]:
            if not (p/method/short/'COMPLETE').exists():continue
            d=pd.read_csv(p/method/short/'cell_results.tsv',sep='\t').set_index('cell_id')
            o[method+'_call']=d['call'].reindex(o.index).fillna('not_run')
            o[method+'_signal']=d.signal_strength.reindex(o.index)
            o[method+'_role']=d.cnv_role.reindex(o.index).fillna('not_run')
            held=d[d.cnv_role.str.startswith('holdout_')];evaluated=held[held.call!='not_evaluable']
            fraction=float((evaluated.call==normal).mean()) if len(evaluated) else 0.
            reference_groups={}
            for role in ['ref_immune','ref_endo']:
                g=d[(d.cnv_role==role)&(d.call!='not_evaluable')]
                reference_groups[role]={'n':len(g),'normal_fraction':float((g.call==normal).mean()) if len(g) else 0.}
            ok=len(evaluated)>=10 and fraction>=.9 and all(g['n']>=10 and g['normal_fraction']>=.9 for g in reference_groups.values())
            qualities[method]=ok
            toolqc.append({'sample':gsm,'short':short,'method':method,'completed':True,'heldout_total':len(held),'heldout_evaluable':len(evaluated),'heldout_normal_fraction':fraction,'reference_immune_evaluable':reference_groups['ref_immune']['n'],'reference_immune_normal_fraction':reference_groups['ref_immune']['normal_fraction'],'reference_endo_evaluable':reference_groups['ref_endo']['n'],'reference_endo_normal_fraction':reference_groups['ref_endo']['normal_fraction'],'reference_check_pass':ok,'notes':'Reference QC checks each lineage; not proof of stromal specificity'})
        for extra in ['copykat_balanced','scevan']:
            if extra+'_call' not in o:
                fail=(p/extra/short/'FAILED').exists()
                ec=p/extra/short/'job.exit_code'
                fail=fail or (ec.exists() and ec.read_text().strip()!='0')
                o[extra+'_call']='failed' if fail else 'not_run';o[extra+'_signal']=np.nan;o[extra+'_role']='not_run'
                if fail:toolqc.append({'sample':gsm,'short':short,'method':extra,'completed':False,'reference_check_pass':False,'notes':'Tool failed, not CNV-negative'})
        priority=['scevan','copykat_balanced','copykat'] if short in ['S0058','SC072'] else ['copykat','scevan']
        selected=next((k for k in priority if qualities.get(k,False)),'copykat')
        selected_good=qualities.get(selected,False)
        normal2='normal' if selected=='scevan' else 'diploid';high2='tumor' if selected=='scevan' else 'aneuploid'
        o['decision_cnv2_source']=selected;o['decision_cnv2_call']=o[selected+'_call'];o['decision_cnv2_signal']=o[selected+'_signal'];o['decision_cnv2_reference_check']=selected_good
        o['decision_cnv2_role']=o[selected+'_role']
        if not selected_good:
            o['decision_cnv2_call']=np.where(o.cnv_role=='not_selected','not_run','failed_QC')
        selections.append({'short':short,'selected_source':selected,'reference_check_pass':selected_good,'rule':'S0058/SC072: independent SCEVAN preferred after original stromal/reference conflict; otherwise balanced or original only if reference QC passes. Other samples: original CopyKAT unless failed reference QC and SCEVAN available.'})
        ic=o.infercnv_call;ck=o.decision_cnv2_call
        low=(ic=='low_signal')&(ck==normal2);high=(ic=='elevated_signal')&(ck==high2)
        evalboth=ic.isin(['low_signal','elevated_signal','ambiguous_signal'])&ck.isin([normal2,high2])
        o['cnv_concordance']='not_jointly_evaluable'
        o.loc[low,'cnv_concordance']='both_low'
        o.loc[high,'cnv_concordance']='both_high'
        o.loc[evalboth&~low&~high,'cnv_concordance']='discordant_or_ambiguous'
        o['original_copykat_conflict']=((ic=='low_signal')&(o.copykat_call=='aneuploid'))|((ic=='elevated_signal')&(o.copykat_call=='diploid'))
        o['cluster_high_fraction']=0.
        for c,z in o.groupby('cluster'):
            ids=z.index;eligible=evalboth.loc[ids]
            fh=float(high.loc[ids][eligible].mean()) if eligible.any() else np.nan
            ic_valid=ic.loc[ids].isin(['low_signal','elevated_signal','ambiguous_signal'])
            risk=fh if eligible.any() else (float((ic.loc[ids][ic_valid]=='elevated_signal').mean()) if ic_valid.any() else np.nan)
            o.loc[ids,'cluster_high_fraction']=risk
            clstats.append({'sample':gsm,'short':short,'cluster':c,'lineage':z.lineage_phase3.iloc[0],'cells':len(z),'selected_cnv2':selected,'both_evaluable':int(eligible.sum()),'both_low':int(low.loc[ids].sum()),'both_high':int(high.loc[ids].sum()),'concordant_high_fraction':fh,'exclusion_risk_fraction':risk,'risk_basis':'two_tool_concordance' if eligible.any() else 'infercnv_only_secondary_unavailable','copykat_aneuploid':int((o.copykat_call.loc[ids]=='aneuploid').sum()),'selected_cnv2_altered':int((ck.loc[ids]==high2).sum()),'infercnv_elevated':int((ic.loc[ids]=='elevated_signal').sum())})
        qc=qualities['infercnv'] and selected_good
        clean=(~o.mixed_lineage_flag)&(~o.doublet_review_flag)&(o.det3_CTA<2)&(~o.cluster_CTA_program)
        nontumor_cluster=o.cluster_high_fraction.fillna(1)<.2
        fiblineage=(o.lineage_phase3=='FibLike')&(o.det3_Fib>=2)
        mesfib=(o.lineage_phase3=='MesAmb')&(o.det3_Fib>=3)&(o.det3_Osteo<2)&(o.det3_Mural<2)
        competing=(o.det3_Mural>=3)|(o.det3_OC>=2)|(o.det3_T_NK>=3)|(o.det3_Myeloid>=4)|(o.det3_Endo>=3)
        o['caf_strict']=fiblineage&low&clean&~competing&nontumor_cluster&(o.n_genes>=500)&qc
        nonmal_one=(ic=='low_signal')|(ck==normal2)
        no_opposing=~ic.isin(['elevated_signal'])&~ck.isin([high2])
        o['caf_extended']=o.caf_strict|((fiblineage|mesfib)&nonmal_one&no_opposing&clean&~competing&nontumor_cluster)
        o['final_identity']='Identity_pending';o['decision_reason']='Insufficient or discordant lineage/CNV evidence'
        other=o.lineage_phase3.isin(['T_NK','B_Plasma','Myeloid','DC','Mast','Endo','Epi','OC_like'])&~high&~o.mixed_lineage_flag&~o.doublet_review_flag
        o.loc[other,'final_identity']='Other';o.loc[other,'decision_reason']='Reviewed non-fibroblast lineage; CNV status retained separately'
        mural=(o.lineage_phase3=='Mural')&low&~o.mixed_lineage_flag&~o.doublet_review_flag&qc
        o.loc[mural,'final_identity']='Mural_SMC';o.loc[mural,'decision_reason']='Mural lineage and concordant low CNV signal'
        mal=high&o.lineage_phase3.isin(['FibLike','MesAmb','OsteoLike','Cycling','Mixed'])&~o.mixed_lineage_flag&~o.doublet_review_flag&qc
        o.loc[mal,'final_identity']='Malignant_OS';o.loc[mal,'decision_reason']='Mesenchymal tumor context and concordant CNV elevation'
        supported=o.caf_strict|(mesfib&low&o.caf_extended&qc)
        o.loc[supported,'final_identity']='Fib_candidate';o.loc[supported,'decision_reason']='Fibroblast lineage with concordant low CNV; non-malignant fibroblast candidate in tumor tissue'
        o['unresolved_issues']=''
        for mask,text in [(o.patient_id=='UNCONFIRMED','patient_mapping_unknown;'),(o.cnv_concordance=='discordant_or_ambiguous','CNV_discordance;'),(o.cnv_concordance=='not_jointly_evaluable','CNV_not_jointly_evaluable;'),(o.mixed_lineage_flag,'mixed_lineage;'),(o.doublet_review_flag,'doublet_review;'),(o.n_genes<500,'low_gene_complexity;')]:o.loc[mask,'unresolved_issues']+=text
        if not qc:o['unresolved_issues']+='reference_QC_not_passed;'
        o.loc[o.cluster_CTA_program&(o.final_identity=='Identity_pending'),'unresolved_issues']+='CTA_cluster_context_unresolved;'
        if selected!='copykat':
            raw_conflict=((ic=='low_signal')&(o.copykat_call=='aneuploid'))|((ic=='elevated_signal')&(o.copykat_call=='diploid'))
            o.loc[raw_conflict,'unresolved_issues']+='original_CopyKAT_conflict_retained;'
        o['dataset']='GSE270231';o['sample_id']=gsm;o['cell_barcode']=o.barcode
        o['reference_strategy']='same_sample_immune_and_endothelial; independent_holdout; selected_cnv2='+selected
        o['infercnv_reference_check']=qualities['infercnv'];o['copykat_reference_check']=qualities['copykat']
        o['caf_set_note']=np.where(o.caf_strict,'strict_candidate',np.where(o.caf_extended,'expanded_candidate_sensitivity_only','not_in_CAF_sets'))
        o.to_csv(out/f'{short}_identity.tsv.gz',sep='\t',index=False)
        final.append(o)
        sm={'sample':gsm,'short':short,'patient_id':o.patient_id.iloc[0],'retained':len(o),'cnv_input':int((o.cnv_role!='not_selected').sum()),'caf_strict':int(o.caf_strict.sum()),'caf_extended':int(o.caf_extended.sum()),'cnv_discordant':int((o.cnv_concordance=='discordant_or_ambiguous').sum()),**o.final_identity.value_counts().to_dict()};summaries.append(sm)
        with h5py.File(a.root/'results'/gsm/'annotated.h5ad') as h:
            original=frame(h['obs']);var=frame(h['var']);counts=sparse(h['layers/counts']);x=sparse(h['X']);umap=h['obsm/X_umap'][()]
        assert o.cell_id.tolist()==original.index.tolist()
        fig,axs=plt.subplots(1,3,figsize=(15,4.7))
        label=o.final_identity.map(IDS)
        for key,color in COL.items():
            z=label==key
            if z.any():axs[0].scatter(umap[z,0],umap[z,1],s=3,c=color,label=key,rasterized=True)
        axs[0].set_title(short+' identity');axs[0].legend(markerscale=3,fontsize=7)
        maps=[{'low_signal':'Low','elevated_signal':'High','ambiguous_signal':'Amb','not_evaluable':'NA','not_run':'NotRun'},{'diploid':'Diploid','aneuploid':'Aneuploid','normal':'Normal','tumor':'Tumor','not_evaluable':'NA','not_run':'NotRun','failed_QC':'BadQC'}]
        for ax,method,mp in zip(axs[1:],['infercnv',selected],maps):
            for key,color in [('NotRun','#DDDDDD'),('NA','#888888'),('BadQC','#888888'),('Amb','#996600'),('Low','#196B31'),('Diploid','#196B31'),('Normal','#196B31'),('High','#A51426'),('Aneuploid','#A51426'),('Tumor','#A51426')]:
                plotcall=o[method+'_call'] if method=='infercnv' else o.decision_cnv2_call
                z=plotcall.map(mp)==key
                if z.any():ax.scatter(umap[z,0],umap[z,1],s=3,c=color,label=key,rasterized=True)
            ax.set_title(short+' '+{'copykat_balanced':'CK_bal','scevan':'SCEVAN'}.get(method,method));ax.legend(markerscale=3,fontsize=7)
        for ax in axs:ax.set_xlabel('UMAP1');ax.set_ylabel('UMAP2')
        fig.tight_layout();fig.savefig(out/'figures'/f'{short}_identity_UMAP.png',dpi=180);plt.close(fig)
        groups=['FibCand','MalOS','Mural','Unk','Other'];points=[]
        for j,g in enumerate(MARKERS):
            ix=np.flatnonzero(var.symbol==g);ct=np.asarray(counts[:,ix].sum(1)).ravel();ex=np.asarray(x[:,ix].sum(1)).ravel()
            for i,lab in enumerate(groups):
                z=(label==lab).values
                if not z.any():continue
                pct=float((ct[z]>0).mean()*100);mean=float(ex[z].mean())
                points.append((j,i,pct,mean));genesummary.append({'sample':gsm,'short':short,'identity':lab,'gene':g,'cells':int(z.sum()),'pct_expressed':pct,'mean_log1p':mean})
        pts=np.asarray(points);fig,ax=plt.subplots(figsize=(12,3.8));im=ax.scatter(pts[:,0],pts[:,1],s=pts[:,2],c=pts[:,3],cmap='Reds',vmin=0)
        ax.set_xticks(range(len(MARKERS)),MARKERS,rotation=90);ax.set_yticks(range(len(groups)),groups);ax.invert_yaxis();ax.set_xlabel('Gene');ax.set_ylabel('Identity');ax.set_title(short+' markers');fig.colorbar(im,ax=ax,label='Mean log1p');fig.tight_layout();fig.savefig(out/'figures'/f'{short}_identity_markers.png',dpi=180);plt.close(fig)
        for flag,objects in [('caf_strict',strict_objs),('caf_extended',extended_objs)]:
            z=o[flag].values
            if not z.any():continue
            obs=o.loc[z].copy();obs.index=obs.cell_id
            # AnnData cannot persist arbitrary mixed object dtypes reliably; keep strings explicit.
            for c in obs.select_dtypes(include='object'):obs[c]=obs[c].fillna('').astype(str)
            obj=ad.AnnData(X=x[z],obs=obs,var=var[['symbol']].copy());obj.layers['counts']=counts[z];obj.obsm['X_umap_phase2']=umap[z];objects.append(obj)
        print('IDENTITY_SAMPLE_INTEGRATED',short,json.dumps(sm),flush=True)
    allcells=pd.concat(final,ignore_index=True);assert len(allcells)==28542 and not allcells.cell_id.duplicated().any()
    allcells.to_csv(out/'all_cell_identity_evidence.tsv.gz',sep='\t',index=False)
    for flag,objects in [('caf_strict',strict_objs),('caf_extended',extended_objs)]:
        allcells.loc[allcells[flag],['dataset','sample_id','short','patient_id','cell_barcode','cell_id','final_identity','caf_set_note','decision_reason','unresolved_issues']].to_csv(out/(flag+'_cells.tsv'),sep='\t',index=False)
        if objects:
            obj=ad.concat(objects,join='inner',merge='same');obj.uns['status']='identity_candidates; CAF reclustering not yet performed';obj.write_h5ad(out/'objects'/(flag+'.h5ad'),compression='gzip')
            check=ad.read_h5ad(out/'objects'/(flag+'.h5ad'),backed='r');assert check.n_obs==int(allcells[flag].sum());check.file.close()
    pd.DataFrame(summaries).fillna(0).to_csv(out/'sample_identity_summary.tsv',sep='\t',index=False)
    pd.DataFrame(toolqc).to_csv(out/'tool_reference_QC.tsv',sep='\t',index=False)
    pd.DataFrame(selections).to_csv(out/'method_selection.tsv',sep='\t',index=False)
    pd.DataFrame(clstats).to_csv(out/'cluster_CNV_concordance.tsv',sep='\t',index=False)
    pd.DataFrame(genesummary).to_csv(out/'identity_marker_expression.tsv',sep='\t',index=False)
    (out/'integration_parameters.json').write_text(json.dumps({'min_heldout_evaluable':10,'min_heldout_normal_fraction':.9,'min_each_reference_lineage_evaluable':10,'min_each_reference_lineage_normal_fraction':.9,'strict_min_genes':500,'exclude_cluster_concordant_high_fraction_at_least':.2,'CTA_cluster_exclusion':'At least 2 predefined cancer-testis genes in existing top50 cluster markers: excluded from both CAF sets, even when individual cells show low CNV signal','strict_definition':'phase3 FibLike + >=2 fibro genes + inferCNV and QC-selected independent method both low + no strong competing/mixed/CTA signal + reference QC pass','expanded_definition':'at least one low tool and no opposing elevation; lineage-supported sensitivity candidates; unresolved cells retained explicitly','method_selection_file':'method_selection.tsv','script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'no_patient_group_inference':True},indent=2),encoding='utf-8')
    (out/'COMPLETE').write_text('IDENTITY_INTEGRATION_COMPLETE\n')
    print('IDENTITY_INTEGRATION_COMPLETE')
if __name__=='__main__':main()
