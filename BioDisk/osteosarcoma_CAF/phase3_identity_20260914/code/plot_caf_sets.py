"""SYNTAX_CHECK: python -m py_compile plot_caf_sets.py
RUN_COMMAND: python plot_caf_sets.py <phase2> <phase3>
EXPECTED_SUCCESS_MARKER: CAF_SET_PLOTS_COMPLETE
"""
import sys
from pathlib import Path
import pandas as pd,h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(sys.argv[1]);p=Path(sys.argv[2])
for f in sorted((p/'identity').glob('*_identity.tsv.gz')):
    o=pd.read_csv(f,sep='\t');short=o.short.iloc[0];sample=o['sample'].iloc[0]
    with h5py.File(root/'results'/sample/'annotated.h5ad') as h:umap=h['obsm/X_umap'][()]
    fig,axs=plt.subplots(1,2,figsize=(10,4.5))
    for ax,field,label in zip(axs,['caf_strict','caf_extended'],['Strict','Expanded']):
        z=o[field].to_numpy();ax.scatter(umap[:,0],umap[:,1],s=2,c='#D5D5D5',rasterized=True)
        ax.scatter(umap[z,0],umap[z,1],s=8,c='#196B31',rasterized=True)
        ax.set_title(f'{short} {label} n={z.sum()}');ax.set_xlabel('UMAP1');ax.set_ylabel('UMAP2')
    fig.tight_layout();fig.savefig(p/'identity/figures'/f'{short}_CAF_sets.png',dpi=180);plt.close(fig)
print('CAF_SET_PLOTS_COMPLETE')
