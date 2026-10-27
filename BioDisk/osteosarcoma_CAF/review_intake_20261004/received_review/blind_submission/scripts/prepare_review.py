import pathlib,hashlib,json
import pandas as pd,numpy as np,h5py
from scipy.io import mmread
P=pathlib.Path(__file__).parent; O=P/'review_output'; O.mkdir(exist_ok=True)
panels={
'T':'CD3D CD3E CD3G TRAC TRBC1 TRBC2 CD247 LCK CD2 IL7R CCR7 FOXP3 CD4 CD8A CD8B',
'NK':'NKG7 GNLY KLRD1 KLRF1 NCAM1 PRF1 GZMB GZMH FCGR3A XCL1 XCL2',
'B':'MS4A1 CD79A CD79B CD19 CD22 CD37 CD74 BANK1 PAX5',
'Plasma':'JCHAIN MZB1 SDC1 XBP1 DERL3 TNFRSF17 IGKC IGHG1 IGHA1',
'Myeloid':'LST1 TYROBP FCER1G AIF1 CSF1R LYZ CTSS SPI1 CD68 C1QA C1QB C1QC APOE CD163 MRC1 SPP1 TREM2 FCN1 S100A8 S100A9 VCAN',
'DC':'CD1C FCER1A CLEC10A CLEC9A XCR1 BATF3 LAMP3 FSCN1 CCR7 GZMB CLEC4C IL3RA TCF4',
'Neutrophil':'CSF3R FCGR3B CXCR2 S100A8 S100A9 MPO ELANE PRTN3 CTSG',
'Mast':'KIT TPSAB1 TPSB2 CPA3 MS4A2 HDC',
'Osteoclast':'CTSK ACP5 ATP6V0D2 DCSTAMP OCSTAMP CALCR MMP9 NFATC1',
'Endothelial':'PECAM1 VWF CDH5 KDR EMCN ESAM RAMP2 ERG CLDN5 PLVAP ACKR1 CA4',
'Lymphatic':'PROX1 PDPN FLT4 LYVE1 CCL21',
'Fibroblast':'DCN LUM COL1A1 COL1A2 COL3A1 DPT PDGFRA COL6A1 COL6A2 COL6A3 C7 PI16 CFD APOD FBLN1 MFAP4 COL14A1 FAP POSTN CTHRC1 CXCL12 LOX MMP14',
'Mural':'RGS5 CSPG4 MCAM NOTCH3 PDGFRB KCNJ8 ABCC9 RBP1 ACTA2 TAGLN MYH11 CNN1 DES MYL9',
'Osteogenic':'RUNX2 SP7 ALPL BGLAP IBSP SPP1 SATB2 IFITM5 PTH1R DMP1 MEPE COL1A1',
'Chondrogenic':'SOX9 COL2A1 ACAN COMP COL9A1 COL9A2 COL9A3 COL10A1 MATN3 CHAD',
'Epithelial':'EPCAM KRT8 KRT18 KRT19 KRT7 KRT5 KRT14 KRT17 TP63 CDH1 MUC1 SFTPC SFTPB SCGB1A1 FOXJ1',
'Glial':'SOX10 S100B PLP1 MPZ MBP GFAP SLC1A3 NGFR',
'Myogenic':'MYOD1 MYOG MYF5 TNNT1 TNNT3 ACTA1 MYH1 MYH2 MYH3 MYH7',
'Adipocyte':'ADIPOQ PLIN1 FABP4 LPL PPARG LEP',
'Erythroid':'HBB HBA1 HBA2 ALAS2 GYPA AHSP KLF1',
'Megakaryocyte':'PPBP PF4 ITGA2B GP9 GP1BA NRGN',
'Cycling':'MKI67 TOP2A UBE2C CENPF BIRC5 PCLAF STMN1 TYMS TK1',
'Stress':'FOS JUN JUNB ATF3 HSPA1A HSPA1B DDIT3 HSPB1 DNAJB1',
'Interferon':'ISG15 IFIT1 IFIT2 IFIT3 MX1 OAS1 STAT1',
'Hypoxia':'CA9 VEGFA SLC2A1 BNIP3 NDRG1 PDK1 LDHA'
}
(O/'marker_panels.json').write_text(json.dumps(panels,ensure_ascii=False,indent=2),encoding='utf-8')
manifest=pd.read_csv(P/'PACKAGE_FILE_MANIFEST.tsv',sep='\t'); check=[]
for r in manifest.itertuples(index=False):
 f=P/r.relative_path; ok=f.is_file(); b=f.stat().st_size if ok else None; sha=hashlib.file_digest(open(f,'rb'),'sha256').hexdigest() if ok else None
 check.append({'path':r.relative_path,'exists':ok,'bytes_match':b==r.bytes,'sha256_match':sha==r.SHA256})
pd.DataFrame(check).to_csv(O/'input_integrity.tsv',sep='\t',index=False)
print('MANIFEST',len(check),pd.DataFrame(check).iloc[:,1:].all().to_dict(),flush=True)
markers=[]; cluster_rows=[]; lines=[]; structs=[]; validation=[]
for d in sorted(P.glob('D0?')):
 c=pd.read_csv(d/'clusters.tsv',sep='\t'); s=pd.read_csv(d/'cell_structure.tsv.gz',sep='\t'); structs.append(s)
 top=pd.read_csv(d/'top_descriptive_markers.tsv',sep='\t'); genes=pd.read_csv(d/'genes.tsv',sep='\t')
 with h5py.File(d/'all_cluster_marker_evidence.h5') as h:
  sy=h['gene_symbol'].asstr()[:]; fi=h['feature_id'].asstr()[:]; ids=h['cluster_id'].asstr()[:]; n=h['N_patients_matched'][:]
  arr={k:h[k][:] for k in ['mean_log_expression','rest_mean_log_expression','detection_fraction','patient_mean_log2_difference','positive_patients','negative_patients']}
  wanted=set(' '.join(panels.values()).split()); ix=np.array([i for i,g in enumerate(sy) if g in wanted])
  for j,cl in enumerate(ids):
   for i in ix: markers.append(dict(source_id=d.name,cluster_id=cl,feature_id=fi[i],gene=sy[i],N_patients_matched=n[j],**{k:v[i,j] for k,v in arr.items()}))
   cr=c.set_index('cluster_id').loc[cl].to_dict(); ss=s[s.r10==cl]; cr.update(source_id=d.name,cluster_id=cl,N_patients_matched=int(n[j]))
   for res in ['r05','r20']:
    vc=ss[res].value_counts(); major=vc.index[0]; cr[res+'_major']=major; cr[res+'_major_fraction']=vc.iloc[0]/len(ss); cr[res+'_n_ge10pct']=int((vc/len(ss)>=.1).sum()); cr[res+'_major_reciprocal']=vc.iloc[0]/(s[res]==major).sum()
   cr['largest_patient']=ss.patient_id.value_counts().index[0]; cr['largest_sample']=ss.sample_id.value_counts().index[0]; cr['tissues']=';'.join(ss.tissue.value_counts().index.astype(str)); cluster_rows.append(cr)
   tt=top[top.cluster_id==cl].head(22)
   lines.append(cl+' n='+str(cr['cells'])+' p='+str(cr['patients'])+' maxP='+format(cr['max_patient_fraction'],'.2f')+' | '+','.join(tt.symbol.fillna(tt.feature).astype(str)))
  m=mmread(d/'primary_cluster_counts.mtx.gz'); sums=np.asarray(m.sum(axis=0)).ravel()
  validation.append({'source':d.name,'cells':len(s),'clusters':len(c),'samples':s.sample_id.nunique(),'patients':s.patient_id.nunique(),'h5_order_matches':list(ids)==list(c.cluster_id),'matrix_rows':m.shape[0],'genes_rows':len(genes),'full_library_matches':bool(np.allclose(sums,c.full_library,rtol=0,atol=0)),'cluster_counts_match':bool(s.r10.value_counts().reindex(c.cluster_id).equals(pd.Series(c.cells.values,index=c.cluster_id)))})
 print(d.name,'done',flush=True)
pd.DataFrame(markers).to_csv(O/'marker_panel_evidence.tsv.gz',sep='\t',index=False)
pd.DataFrame(cluster_rows).sort_values('cluster_id').to_csv(O/'cluster_structure_review.tsv',sep='\t',index=False)
pd.concat(structs).to_pickle(O/'cell_structure.pkl')
pd.DataFrame(validation).to_csv(O/'data_reconciliation.tsv',sep='\t',index=False)
(O/'all_top_markers.txt').write_text('\n'.join(sorted(lines)),encoding='utf-8')
print(pd.DataFrame(validation).to_string(index=False),flush=True)
