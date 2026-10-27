from pathlib import Path
import copy
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
B=Path(__file__).parent;O=B/'report_update';O.mkdir(exist_ok=True)
g=pd.read_csv(B/'inputs/reused_full_transcriptome_withinF_GSEA.tsv.gz',sep='\t')
s=pd.read_csv(B/'inputs/subgroup_summary.tsv',sep='\t')
source=(B/'update_report.py').read_text(encoding='utf-8')
start=source.index("terms=['REACTOME_")
stop=source.index('nom.to_csv',start)
exec(source[start:stop])
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
start=source.index('order=s.subgroup_id.tolist()')
stop=source.index("assert (O/'GSEA_33_all.png')",start)
exec(source[start:stop])
heatmap(order,'all',O/'GSEA_33_all.png',(8.3,9.5))
heatmap(focus,'Fibro',O/'GSEA_Fibro_focus.png',(8.3,4.0))
print('Two figures saved')
