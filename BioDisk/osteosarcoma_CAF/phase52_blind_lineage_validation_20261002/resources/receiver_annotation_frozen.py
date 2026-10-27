provider=meta.get('provider_consensus_celltype_annotation',pd.Series('',index=meta.index)).astype(str).str.lower();hist=meta.get('historical_broad_lineage',pd.Series('',index=meta.index)).astype(str).str.lower()
known={
'T_NK':provider.str.contains(r'\bt cell\b|natural killer|innate lymphoid',regex=True).to_numpy()|hist.isin(['t_nk','t','nk','tnk']).to_numpy(),
'Myeloid':provider.str.contains(r'macrophage|monocyte|myeloid|neutrophil|dendritic|phagocyte',regex=True).to_numpy()|hist.isin(['myeloid','dc']).to_numpy(),
'Endothelial':provider.str.contains('endothelial').to_numpy()|hist.isin(['endo','endothelial']).to_numpy(),
'Mural':provider.str.contains(r'pericyte|smooth muscle|perivascular',regex=True).to_numpy()|hist.isin(['mural','pericyte','smooth_muscle']).to_numpy(),
'Fibroblast_other':provider.str.contains(r'fibroblast|stromal cell',regex=True).to_numpy()|hist.isin(['fiblike','fibroblast','fibroblast_like']).to_numpy()}
bits=np.column_stack(list(known.values()));ambiguous=bits.sum(1)>1;known_union=bits.any(1);mask={k:v&~ambiguous&nonfrozen for k,v in known.items()}
evidence=pd.read_csv(P/'phase26_client_reanalysis_20260922/tables'/f'{row.sample}_identity_evidence.tsv.gz',sep='\t',keep_default_na=False).set_index('cell_id');assert evidence.index.is_unique and set(evidence.index)==set(meta.index)
alteration=evidence.CNV_evidence.reindex(meta.index).astype(str).eq('concordant_alteration').to_numpy();reference=flag(cells.loc[meta.index,'reference']).to_numpy();mask['CNV_altered_unresolved']=alteration&~reference&~known_union&nonfrozen;used=np.column_stack(list(mask.values())).any(1);mask['Other_nonfrozen_unresolved']=nonfrozen&~used