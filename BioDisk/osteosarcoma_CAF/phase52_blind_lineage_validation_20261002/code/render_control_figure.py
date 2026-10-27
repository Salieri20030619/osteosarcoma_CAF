from pathlib import Path
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
R=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase52_blind_lineage_validation_20261002')
b=pd.read_csv(R/'tables/control_agreement_summary.tsv',sep='\t').sort_values('annotation')
fig,ax=plt.subplots(figsize=(10,6.5));y=np.arange(len(b))
ax.barh(y-.18,b.mean_full_agreement*100,height=.34,label='Full-panel agreement',color='#467994')
ax.barh(y+.18,b.mean_correct_consensus_fraction*100,height=.34,label='All controls: agreeing consensus',color='#6eae99')
ax.set_yticks(y,[f'{r.annotation} (N={r.N_patients})' for r in b.itertuples()]);ax.invert_yaxis();ax.set_xlim(0,110)
for j,r in enumerate(b.itertuples()):
 ax.text(r.mean_full_agreement*100+.8,j-.18,f'{100*r.mean_full_agreement:.1f}%',va='center',fontsize=8)
 ax.text(r.mean_correct_consensus_fraction*100+.8,j+.18,f'{100*r.mean_correct_consensus_fraction:.1f}%',va='center',fontsize=8)
ax.set_xlabel('Patient-mean agreement with existing annotations (%)');ax.set_title('Target-gene-heldout panel benchmark');ax.legend(loc='lower right',fontsize=9,frameon=False)
fig.tight_layout(rect=[0,.055,1,1]);fig.text(.07,.015,'Annotation agreement is not external accuracy. Osteogenic controls: none.',fontsize=9)
fig.savefig(R/'figures/control_annotation_agreement.png',dpi=170);fig.savefig(R/'figures/control_annotation_agreement.pdf');plt.close(fig)
