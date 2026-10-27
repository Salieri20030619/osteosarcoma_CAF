import pathlib,json,pandas as pd,numpy as np,h5py
from scipy.io import mmread
P=pathlib.Path(__file__).parent;O=P/'review_output'
core={
'T':'CD3D CD3E CD3G TRAC TRBC1 TRBC2', 'NK':'GNLY KLRD1 KLRF1 NKG7 PRF1', 'B':'MS4A1 CD79A CD79B CD22 PAX5', 'Plasma':'JCHAIN MZB1 DERL3 TNFRSF17 IGKC', 'Myeloid':'LST1 TYROBP FCER1G CSF1R AIF1 LYZ C1QA C1QB CD163', 'DC':'CD1C FCER1A CLEC10A CLEC9A XCR1 LAMP3', 'pDC':'LILRA4 CLEC4C IL3RA GZMB PLD4 TCF4', 'Mast':'TPSAB1 TPSB2 KIT CPA3 MS4A2 HDC', 'Osteoclast':'ACP5 CTSK ATP6V0D2 MMP9 DCSTAMP TCIRG1', 'Endothelial':'PECAM1 CDH5 VWF KDR EMCN CLDN5', 'Lymphatic':'PROX1 FLT4 CCL21 LYVE1', 'Fibroblast':'DCN LUM DPT PDGFRA COL3A1 COL14A1 FBLN1 PI16', 'Mural':'RGS5 NOTCH3 PDGFRB ABCC9 KCNJ8 MCAM', 'SmoothMuscle':'MYH11 CNN1 ACTA2 TAGLN LMOD1', 'Osteogenic':'SP7 RUNX2 SATB2 ALPL IBSP IFITM5 PTH1R DMP1 MEPE BGLAP', 'Chondrogenic':'SOX9 COL2A1 ACAN COL9A1 COL9A2 COL9A3 COL10A1', 'Epithelial':'EPCAM KRT8 KRT18 KRT19 CDH1 KRT7', 'AT2':'SFTPC SFTPB SFTPA1 SFTPA2 ABCA3 NAPSA SLC34A2', 'AT1':'AGER CAV1 CLIC5 EMP2 AQP5 PDPN', 'Ciliated':'FOXJ1 PIFO TPPP3 DNAH5 DNAH12 RSPH1 CFAP157', 'Renal':'PAX8 LRP2 CUBN SLC34A1 SLC17A1 AQP2 FXYD4 SLC4A9 SLC26A4', 'Glial':'SOX10 S100B PLP1 MPZ MBP GFAP', 'Myogenic':'MYOG MYOD1 MYF5 PAX7 ACTA1 TNNT3 MYH3 MYH7', 'Erythroid':'HBB HBA1 HBA2 ALAS2 GYPA AHSP', 'Platelet':'PF4 PPBP GP9 ITGA2B', 'Cycling':'MKI67 TOP2A UBE2C CENPF PCLAF TK1', 'Stress':'FOS JUN JUNB ATF3 HSPA1A DDIT3', 'IFN':'ISG15 IFIT1 IFIT2 IFIT3 MX1 OAS1'}
(O/'core_panels.json').write_text(json.dumps(core,indent=2),encoding='utf-8')
rows=[]; summaries=[]; qc=[]
for d in sorted(P.glob('D0?')):
 top=pd.read_csv(d/'top_descriptive_markers.tsv',sep='\t'); c=pd.read_csv(d/'clusters.tsv',sep='\t');
 with h5py.File(d/'all_cluster_marker_evidence.h5') as h:
  sy=h['gene_symbol'].asstr()[:];fi=h['feature_id'].asstr()[:];ids=h['cluster_id'].asstr()[:];N=h['N_patients_matched'][:]; a={k:h[k][:] for k in ['mean_log_expression','rest_mean_log_expression','detection_fraction','patient_mean_log2_difference','positive_patients','negative_patients']};
  wanted=set(' '.join(core.values()).split())|set(' '.join(json.loads((O/'marker_panels.json').read_text()).values()).split()); gmap={g:np.where(sy==g)[0] for g in wanted}
  for j,cl in enumerate(ids):
   wantedj=wanted|set(top.loc[top.cluster_id==cl,'symbol'].dropna())
   for i,g in enumerate(sy):
    if g in wantedj: rows.append(dict(source_id=d.name,cluster_id=cl,feature_id=fi[i],gene=g,N_patients_matched=int(N[j]),**{k:v[i,j] for k,v in a.items()}))
   panel=[]
   for name,gs in core.items():
    ii=np.concatenate([gmap[g] for g in gs.split()]);
    if not len(ii):continue
    ii=sorted(ii,key=lambda i:a['detection_fraction'][i,j]*max(a['patient_mean_log2_difference'][i,j],0),reverse=True)
    score=np.mean([a['detection_fraction'][i,j]*max(a['patient_mean_log2_difference'][i,j],0) for i in ii]); active=[i for i in ii if a['patient_mean_log2_difference'][i,j]>.2 and a['detection_fraction'][i,j]>.03]
    panel.append((score,name,active,ii))
   panel=sorted(panel,reverse=True)
   ss=[]
   for score,name,active,ii in panel[:5]:ss.append(name+':'+','.join(sy[i]+'('+format(a['detection_fraction'][i,j],'.2f')+'/'+format(a['patient_mean_log2_difference'][i,j],'.1f')+')' for i in ii[:3]))
   summaries.append(cl+' | '+' | '.join(ss))
  m=mmread(d/'primary_cluster_counts.mtx.gz').tocsr(); total=np.asarray(m.sum(axis=0)).ravel(); mt=np.char.startswith(sy.astype(str),'MT-');ribo=np.array([g.startswith(('RPL','RPS')) for g in sy]);
  for j,cl in enumerate(ids): qc.append({'source_id':d.name,'cluster_id':cl,'matrix_library':total[j],'reported_full_library':c.full_library.iloc[j],'full_minus_matrix':c.full_library.iloc[j]-total[j],'matrix_mt_fraction':np.asarray(m[mt].sum(axis=0)).ravel()[j]/total[j],'matrix_ribosomal_prefix_fraction':np.asarray(m[ribo].sum(axis=0)).ravel()[j]/total[j],'counts_per_cell':total[j]/c.cells.iloc[j]})
pd.DataFrame(rows).to_csv(O/'selected_gene_evidence.tsv.gz',sep='\t',index=False)
pd.DataFrame(qc).to_csv(O/'cluster_expression_qc.tsv',sep='\t',index=False)
(O/'panel_review.txt').write_text('\n'.join(sorted(summaries)),encoding='utf-8')
print('COMPLETE',len(rows))
print(pd.DataFrame(qc).groupby('source_id')[['full_minus_matrix','matrix_mt_fraction']].agg(['min','median','max']).to_string())
print((P/'D08/technical_summary.json').read_text())
