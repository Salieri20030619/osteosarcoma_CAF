"""Review lineage and export full-RNA counts for independent CNV tools.
SYNTAX_CHECK: python -m py_compile review_inputs.py
RUN_COMMAND: python review_inputs.py --root <phase2> --out <phase3>
EXPECTED_SUCCESS_MARKER: CNV_INPUT_REVIEW_COMPLETE
"""
import argparse,json,gzip,hashlib
from pathlib import Path
import numpy as np,pandas as pd,h5py
from scipy.io import mmwrite
from takeover_audit import frame,sparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
SEED=270231
PANELS={
'T_NK':['CD2','CD3D','CD3E','CD247','BCL11B','THEMIS','NKG7','KLRD1','GNLY'],
'B_Plasma':['CD79A','MS4A1','CD79B','MZB1','POU2AF1','JCHAIN'],
'Myeloid':['LYZ','LST1','MS4A7','CD163','CSF1R','MERTK','MSR1','SPI1','TYROBP'],
'DC':['CLEC4C','LILRA4','FLT3','FCER1A','IRF8'],
'Endo':['PECAM1','VWF','CLDN5','EMCN','KDR','TEK','EGFL7','CDH5','PTPRB'],
'Epi':['EPCAM','KRT8','KRT18','KRT19','SFTPC','SFTPB','AGER','ABCA3'],
'Fib':['DCN','LUM','PDGFRA','DPT','COL3A1','COL6A1','COL14A1','PI16','TCF21'],
'Osteo':['RUNX2','SP7','ALPL','IBSP','BGLAP','DMP1','SATB2','PHEX','PANX3'],
'Mural':['RGS5','CSPG4','NOTCH3','MCAM','MYH11','PDGFRB','ABCC9','KCNJ8','COX4I2'],
'OC':['ACP5','CTSK','ATP6V0D2','OCSTAMP','DCSTAMP','TNFRSF11A','CALCR','RUFY4'],
'Cycle':['MKI67','TOP2A','AURKB','ASPM','CDK1','CENPF','UBE2C'],
'CTA':['MAGEA3','MAGEA10','MAGEA12','SSX1','CTAG2','PAGE2B','PAGE5','CSAG1','FMR1NB'],
'ECM':['COL1A1','COL1A2','COL3A1','COL11A1','FAP','POSTN','THBS2','MMP2','SPARC'],
'Contract':['ACTA2','TAGLN','MYL9','TPM2','CNN1','MYH11'],
'Mineral':['MGP','ENPP1','STC1','ALPL','IBSP','BGLAP','DMP1','PHEX','SPP1'],
}
COLORS={'T_NK':'#165DB1','B_Plasma':'#8545AB','Myeloid':'#C15B00','OC_like':'#964329','Endo':'#007D78','Epi':'#C8283B','Mural':'#C6288B','FibLike':'#238329','OsteoLike':'#8D7500','MesAmb':'#495B9C','Mixed':'#555555','Unresolved':'#999999','Mast':'#C45C7E','DC':'#498500','Cycling':'#7A3567'}
# Fixed review proposals; require at least two named markers in the existing top-50 list.
REVIEW={
('S0217','1'):('T_NK',['THEMIS','CD2','SKAP1'],'T-lineage cluster previously missed'),
('S0217','8'):('Myeloid',['MSR1','MS4A7','DOCK2'],'Myeloid cluster previously missed'),
('S0058','3'):('Mural',['MYH11','ACTG2','AVPR1A'],'Smooth muscle lineage'),
('S0059','18'):('Mural',['COX4I2','HIGD1B','LMOD1'],'Mural lineage'),
('S0217','9'):('Mural',['EBF2','ABCC9','CDH6'],'Mural candidate; retain CNV target'),
('SC069','11'):('Mural',['RGS5','PDGFRB','MCAM'],'Pericyte lineage'),
('SC072','9'):('Mural',['PDGFRB','TRPC6','HEYL'],'Mural lineage'),
('SC073','10'):('Mural',['RGS5','PDGFRB','TRPC6'],'Mural lineage'),
('S0058','17'):('MesAmb',['CTAG2','PAGE2B','PAGE5'],'FibLike with cancer-testis markers; malignancy unresolved'),
('SC069','5'):('Cycling',['ASPM','MKI67','TOP2A'],'Cycling mesenchymal cells; fibroblast identity unresolved'),
('S0059','21'):('Cycling',['NEK2','KIF20A','PLK1'],'Cycling mesenchymal cells'),
('S0217','4'):('Cycling',['PBK','AURKB','MELK'],'Cycling mesenchymal cells'),
('S0218','6'):('Cycling',['TOP2A','ASPM','MKI67'],'Cycling mesenchymal cells'),
('SC073','8'):('Cycling',['CENPF','TOP2A','HIST1H2BG'],'Cycling mesenchymal cells'),
('S0218','0'):('Mixed',['ABCA3','SFTPB','LRP2'],'Alveolar program conflicts with OsteoLike; retain target'),
}
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    for d in ['inputs','review','figures']: (a.out/d).mkdir(parents=True,exist_ok=True)
    summary=[];audit=[];clusters=[];evidence=[];coverage=[]
    samples=pd.read_csv(a.root/'results/sample_summary.tsv',sep='\t')
    for r in samples.itertuples():
        ip=a.root/'results'/r.sample
        with h5py.File(ip/'annotated.h5ad') as h:
            o=frame(h['obs']);v=frame(h['var']);x=sparse(h['layers/counts']);logx=sparse(h['X']);umap=h['obsm/X_umap'][()]
        o['cluster']=o.cluster.astype(str);o['phase2_label']=o.lineage_reviewed;o['lineage_phase3']=o.lineage_reviewed
        mk=pd.read_csv(ip/'cluster_markers.tsv',sep='\t',dtype={'group':str})
        for (short,c),(label,genes,reason) in REVIEW.items():
            if short!=r.short:continue
            hits=set(genes)&set(mk.loc[mk.group==c,'symbol']);assert len(hits)>=2,(short,c,hits)
            o.loc[o.cluster==c,'lineage_phase3']=label
            audit.append({'sample':r.sample,'short':short,'cluster':c,'phase2_label':o.loc[o.cluster==c,'phase2_label'].iloc[0],'phase3_label':label,'marker_hits':';'.join(sorted(hits)),'reason':reason})
        for k,genes in PANELS.items():
            ix=np.flatnonzero(v.symbol.isin(genes));mat=x[:,ix]
            o['det3_'+k]=np.asarray((mat>0).sum(1)).ravel()
            o['pctUMI_'+k]=np.asarray(mat.sum(1)).ravel()/o.n_counts*100
            o['meanlog_'+k]=np.asarray(logx[:,ix].sum(1)).ravel()/max(len(ix),1)
            coverage.append({'sample':r.sample,'module':k,'requested':len(genes),'covered':len(set(v.symbol.iloc[ix])),'missing':';'.join(sorted(set(genes)-set(v.symbol)))})
        # Follow the original conservative mesenchymal exclusion for reference cells.
        # More permissive low-depth cells remain targets/candidates, never forced into normal references.
        mes=(o.det_FibLike>=2)|(o.det_OsteoLike>=2)|(o.det_Mural>=2)
        immune_good=np.zeros(len(o),dtype=bool)
        for k in ['T_NK','B_Plasma','Myeloid','DC']:immune_good |= ((o.lineage_phase3==k)&(o['det3_'+k]>=2)).values
        ref_qc=(o.n_genes>=500)&(o.n_counts>=750)&(o.doublet_score<=o.doublet_score.quantile(.9))
        clean=~mes&(o.det3_OC<2)&(o.det3_Epi<2)
        o['ref_eligible_immune']=immune_good&clean&ref_qc&(o.det3_Endo<2)
        o['ref_eligible_endo']=(o.lineage_phase3=='Endo')&(o.det3_Endo>=3)&~mes&(o.det3_OC<2)&(o.det3_Epi<2)&ref_qc
        o['cnv_target_expanded']=o.cnv_target_reviewed|o.lineage_phase3.isin(['FibLike','MesAmb','Mural','OsteoLike','Mixed','Unresolved','Cycling','OC_like'])
        # Flag combined lineage evidence; no fixed-fraction doublet deletion.
        immune_strong=(o.det3_T_NK>=3)|(o.det3_B_Plasma>=3)|(o.det3_Myeloid>=4)
        mes_strong=(o.det3_Fib>=4)|(o.det3_Osteo>=3)|(o.det3_Mural>=4)
        o['mixed_lineage_flag']=(immune_strong&mes_strong)|((o.det3_Endo>=4)&(o.det3_Osteo>=3))
        o['doublet_review_flag']=o.mixed_lineage_flag & ((o.doublet_score>=o.doublet_score.quantile(.9))|o.high_complexity_flag)
        o['cnv_role']='not_selected';o.loc[o.cnv_target_expanded,'cnv_role']='target'
        rng=np.random.default_rng(SEED)
        for typ,cap in [('immune',400),('endo',160)]:
            ids=o.index[o['ref_eligible_'+typ]].to_numpy();rng.shuffle(ids)
            ids=ids[:cap]
            n_hold=max(5,int(len(ids)*.2)) if len(ids)>=30 else 0
            o.loc[ids[:n_hold],'cnv_role']='holdout_'+typ
            o.loc[ids[n_hold:],'cnv_role']='ref_'+typ
        refs=o.cnv_role.str.startswith('ref_');holds=o.cnv_role.str.startswith('holdout_')
        assert refs.sum()>=30,(r.short,'insufficient references',refs.sum())
        sel=(o.cnv_role!='not_selected')
        out=a.out/'inputs'/r.short;out.mkdir(exist_ok=True)
        # Tool-safe IDs avoid R name mangling; retain exact full sample barcode mapping.
        o['tool_id']=[f'{r.short}_c{i:05d}' for i in range(len(o))]
        o['cell_id']=o.index;o['short']=r.short
        o['patient_mapping_status']=np.where(o.patient_id=='UNCONFIRMED','unresolved','known_within_dataset')
        o['phase3_note']='CNV pending; labels are lineage candidates'
        o.to_csv(a.out/'review'/f'{r.short}_cell_evidence.tsv.gz',sep='\t',index=False)
        o.loc[sel].to_csv(out/'cells.tsv',sep='\t',index=False)
        genes=v[['symbol']].copy();genes.insert(0,'ensembl_id',v.index);genes.to_csv(out/'genes.tsv',sep='\t',index=False)
        with gzip.open(out/'counts.mtx.gz','wb') as f:mmwrite(f,x[sel].T,field='integer')
        ann=pd.DataFrame({'tool_id':o.loc[sel,'tool_id'],'annotation':np.where(refs[sel],o.loc[sel,'cnv_role'],'C'+o.loc[sel,'cluster']+'_'+o.loc[sel,'lineage_phase3'])})
        ann.to_csv(out/'annotations.tsv',sep='\t',index=False,header=False)
        summary.append({'sample':r.sample,'short':r.short,'retained':len(o),'phase2_FibLike':int((o.phase2_label=='FibLike').sum()),'phase3_FibLike':int((o.lineage_phase3=='FibLike').sum()),'target_expanded_total':int(o.cnv_target_expanded.sum()),'target_run':int((o.cnv_role=='target').sum()),'ref_immune':int((o.cnv_role=='ref_immune').sum()),'ref_endo':int((o.cnv_role=='ref_endo').sum()),'holdout':int(holds.sum()),'run_cells':int(sel.sum()),'doublet_review_flag':int(o.doublet_review_flag.sum()),'ref_strategy':'same_sample; heldout_validation'})
        for c,cc in o.groupby('cluster'):
            ci=(o.cluster==c).values
            cr={'sample':r.sample,'short':r.short,'cluster':c,'cells':len(cc),'phase2_label':cc.phase2_label.iloc[0],'lineage_phase3':cc.lineage_phase3.iloc[0],'ref_eligible':int((cc.ref_eligible_immune|cc.ref_eligible_endo).sum()),'doublet_review_flag':int(cc.doublet_review_flag.sum())}
            for k in PANELS:cr[k+'_pct2']=float((cc['det3_'+k]>=2).mean()*100);cr[k+'_pctUMI']=float(cc['pctUMI_'+k].median())
            clusters.append(cr)
            for g in sorted(set(sum(PANELS.values(),[]))):
                ix=np.flatnonzero(v.symbol==g);vals=np.asarray(x[:,ix].sum(1)).ravel();expr=np.asarray(logx[:,ix].sum(1)).ravel()
                evidence.append({'sample':r.sample,'short':r.short,'cluster':c,'gene':g,'pct':float((vals[ci]>0).mean()*100),'meanlog':float(expr[ci].mean())})
        fig,axs=plt.subplots(1,3,figsize=(15,4.5))
        for label,color in COLORS.items():
            z=o.lineage_phase3==label
            if z.any():axs[0].scatter(umap[z,0],umap[z,1],s=3,c=color,label=label,rasterized=True)
        axs[0].legend(fontsize=7,markerscale=3,loc='best');axs[0].set_title(r.short+' lineage')
        rc={'target':'#555555','not_selected':'#DDDDDD','ref_immune':'#165DB1','ref_endo':'#007D78','holdout_immune':'#D69C00','holdout_endo':'#AA4499'}
        for role,color in rc.items():
            z=o.cnv_role==role
            if z.any():axs[1].scatter(umap[z,0],umap[z,1],s=3,c=color,label=role.replace('not_selected','Other').replace('holdout_','HO_'))
        axs[1].legend(fontsize=7,markerscale=3,loc='best');axs[1].set_title(r.short+' CNV input')
        im=axs[2].scatter(umap[:,0],umap[:,1],s=3,c=o.doublet_score,cmap='magma');fig.colorbar(im,ax=axs[2],label='DB score');axs[2].set_title(r.short+' DB score')
        for ax in axs:ax.set_xlabel('UMAP1');ax.set_ylabel('UMAP2')
        fig.tight_layout();fig.savefig(a.out/'figures'/f'{r.short}_review.png',dpi=170);plt.close(fig)
        print('INPUT_REVIEW_DONE',r.short,json.dumps(summary[-1]),flush=True)
    pd.DataFrame(summary).to_csv(a.out/'review/input_summary.tsv',sep='\t',index=False)
    pd.DataFrame(audit).to_csv(a.out/'review/annotation_changes.tsv',sep='\t',index=False)
    pd.DataFrame(clusters).to_csv(a.out/'review/cluster_evidence.tsv',sep='\t',index=False)
    pd.DataFrame(evidence).to_csv(a.out/'review/marker_evidence.tsv.gz',sep='\t',index=False)
    pd.DataFrame(coverage).to_csv(a.out/'review/panel_coverage.tsv',sep='\t',index=False)
    (a.out/'review/parameters.json').write_text(json.dumps({'seed':SEED,'panels':PANELS,'reference_min_genes':500,'reference_min_umi':750,'reference_max_doublet_quantile':.9,'reference_caps':{'immune':400,'endo':160},'holdout_fraction':.2,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2),encoding='utf-8')
    print('CNV_INPUT_REVIEW_COMPLETE')
if __name__=='__main__':main()
