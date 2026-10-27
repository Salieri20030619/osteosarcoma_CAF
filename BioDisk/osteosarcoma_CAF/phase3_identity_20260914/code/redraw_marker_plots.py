"""Redraw saved marker summaries with cell counts and a detection-size legend.
SYNTAX_CHECK: python -m py_compile redraw_marker_plots.py
RUN_COMMAND: python redraw_marker_plots.py <phase3>
EXPECTED_SUCCESS_MARKER: MARKER_PLOTS_VERIFIED_LAYOUT
"""
import sys
from pathlib import Path
import pandas as pd,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=Path(sys.argv[1]);data=pd.read_csv(p/'identity/identity_marker_expression.tsv',sep='\t')
genes=['DCN','LUM','PDGFRA','FAP','COL1A1','COL11A1','ACTA2','TAGLN','MYH11','RGS5','RUNX2','ALPL','IBSP','MGP','ENPP1','SPP1','CD36','PECAM1','PTPRC']
for short,d in data.groupby('short',sort=False):
    groups=[x for x in ['FibCand','MalOS','Mural','Unk','Other'] if x in d.identity.values]
    rows=[];labels=[]
    for i,group in enumerate(groups):
        sub=d[d.identity==group].set_index('gene').loc[genes]
        labels.append(f'{group} ({int(sub.cells.iloc[0])})')
        for j,r in enumerate(sub.itertuples()):rows.append((j,i,r.pct_expressed,r.mean_log1p))
    v=np.array(rows)
    fig=plt.figure(figsize=(13.5,3.8));gs=fig.add_gridspec(1,3,width_ratios=[20,.4,1.8],wspace=.5)
    ax=fig.add_subplot(gs[0,0]);cx=fig.add_subplot(gs[0,1]);lx=fig.add_subplot(gs[0,2]);lx.axis('off')
    im=ax.scatter(v[:,0],v[:,1],s=v[:,2],c=v[:,3],cmap='Reds',vmin=0)
    ax.set_xticks(range(len(genes)),genes,rotation=90);ax.set_yticks(range(len(groups)),labels);ax.invert_yaxis();ax.set_xlabel('Gene');ax.set_ylabel('Identity');ax.set_title(short+' markers')
    fig.colorbar(im,cax=cx,label='Mean log1p')
    handles=[lx.scatter([],[],s=n,c='#8A2222',label=str(n)) for n in [25,50,100]]
    lx.legend(handles=handles,title='Expr %',loc='center',frameon=False,labelspacing=1.1)
    fig.subplots_adjust(left=.1,right=.97,bottom=.3,top=.87)
    fig.savefig(p/'identity/figures'/f'{short}_identity_markers.png',dpi=180);plt.close(fig)
print('MARKER_PLOTS_VERIFIED_LAYOUT')
