"""SYNTAX_CHECK: python -m py_compile plot_cnv_summary.py
RUN_COMMAND: python plot_cnv_summary.py <phase3>
EXPECTED_SUCCESS_MARKER: CNV_SUMMARY_PLOTS_COMPLETE
"""
import sys,json
from pathlib import Path
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from matplotlib.patches import Patch
from review_inputs import COLORS
p=Path(sys.argv[1]);summary=[]
for f in sorted((p/'identity/plot_data').glob('*_profiles.tsv.gz')):
    short,method=f.name.removesuffix('_profiles.tsv.gz').split('_',1)
    z=pd.read_csv(f,sep='\t');chrom=z['chr'].astype(str)
    ids=z.columns[2:];o=pd.read_csv(p/method/short/'cell_results.tsv',sep='\t').set_index('tool_id').loc[ids].copy()
    o['role_order']=o.cnv_role.map(lambda x:0 if x.startswith('ref_') else (1 if x.startswith('holdout_') else (2 if x.startswith('check_') else 3)))
    o=o.sort_values(['role_order','lineage_phase3','cluster']);v=z[o.index].to_numpy().T
    fig=plt.figure(figsize=(13,7.5));gs=fig.add_gridspec(1,3,width_ratios=[.14,.14,10],wspace=.02)
    ax1=fig.add_subplot(gs[0,0]);ax2=fig.add_subplot(gs[0,1]);ax=fig.add_subplot(gs[0,2])
    rolecolors={0:'#165DB1',1:'#D69C00',2:'#92278F',3:'#555555'}
    ax1.imshow(np.array([to_rgb(rolecolors[x]) for x in o.role_order])[:,None,:],aspect='auto',interpolation='nearest')
    ax2.imshow(np.array([to_rgb(COLORS[x]) for x in o.lineage_phase3])[:,None,:],aspect='auto',interpolation='nearest')
    for aa in [ax1,ax2]:aa.set_xticks([]);aa.set_yticks([])
    ax1.set_ylabel('Cells')
    ax1.text(.5,1.012,'Role',rotation=90,ha='center',va='bottom',fontsize=8,transform=ax1.transAxes)
    ax2.text(.5,1.012,'Type',rotation=90,ha='center',va='bottom',fontsize=8,transform=ax2.transAxes)
    limit=.15 if method.startswith('copykat') else .25
    im=ax.imshow(v,aspect='auto',interpolation='nearest',cmap='RdBu_r',vmin=-limit,vmax=limit,rasterized=True)
    positions=[];labs=[]
    for c in chrom.drop_duplicates():
        ix=np.flatnonzero(chrom==c);positions.append((ix[0]+ix[-1])/2);labs.append(c)
        if ix[0]>0:ax.axvline(ix[0]-.5,color='#333333',linewidth=.4)
    ax.set_xticks(positions,labs,rotation=90,fontsize=8);ax.set_yticks([]);ax.set_xlabel('Chr');ax.set_title(short+' '+{'infercnv':'inferCNV','copykat':'CopyKAT','copykat_balanced':'CK_bal','scevan':'SCEVAN'}[method])
    boundaries=np.flatnonzero(np.diff(o.role_order.to_numpy()))+.5
    for b in boundaries:ax.axhline(b,color='black',linewidth=.6)
    fig.colorbar(im,ax=ax,shrink=.6,pad=.015,label='Log2 dev' if method=='infercnv' else 'CNV dev')
    labels=list(dict.fromkeys(o.lineage_phase3));handles=[Patch(color=rolecolors[k],label=lab) for k,lab in [(0,'Ref'),(1,'Hold'),(2,'Check'),(3,'Target')] if k in o.role_order.values]+[Patch(color=COLORS[k],label=k) for k in labels]
    fig.legend(handles=handles,loc='lower center',ncol=8,fontsize=8,frameon=False)
    fig.subplots_adjust(bottom=.14,top=.91,left=.05,right=.98)
    dst=p/'identity/figures'/f'{short}_{method}_CNV.png';fig.savefig(dst,dpi=180);plt.close(fig)
    summary.append({'short':short,'method':method,'displayed_cells':len(o),'displayed_bins':len(z),'selection':'seed270231; at most25 cells per actual tool role/cluster; tool-evaluable cells only','aggregation':'inferCNV/SCEVAN 100-gene windows; CopyKAT median of20 adjacent genomic bins','symmetric_color_limit':limit})
pd.DataFrame(summary).to_csv(p/'identity/CNV_plot_manifest.tsv',sep='\t',index=False)
print('CNV_SUMMARY_PLOTS_COMPLETE')
