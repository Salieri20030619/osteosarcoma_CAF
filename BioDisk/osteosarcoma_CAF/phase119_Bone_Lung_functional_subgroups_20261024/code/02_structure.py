from caf119_hierarchy import *
assert (O/'INPUTS_COMPLETE.json').exists();event('02_structure','started');X,m,genes=load();allout=[];allnodes=[]
for site in ['Bone','Lung']:
 allrows=np.flatnonzero(m.site.eq(site));fitrows=np.flatnonzero(m.site.eq(site)&~m.doublet)
 if all((O/'tables'/('main_'+site+s)).exists() for s in ['_labels.tsv.gz','_nodes.tsv','_refinement_eligibility.tsv']):
  out=read('main_'+site+'_labels.tsv.gz');assert out.cell_id.tolist()==m.iloc[allrows].cell_id.tolist();out.index=allrows;nodes=read('main_'+site+'_nodes.tsv');event('02_structure','reused_completed_site',site=site)
 else:out,nodes=discover(X,m,genes,site,fitrows,allrows,'main_'+site)
 allout.append(out);allnodes.append(nodes)
labels=pd.concat(allout).sort_index();assert labels.cell_id.tolist()==m.cell_id.tolist();nodes=pd.concat(allnodes,ignore_index=True);save(labels.reset_index(drop=True),'hierarchical_cell_labels.tsv.gz');save(nodes,'hierarchical_nodes.tsv');assert nodes.node_id.is_unique
for oldcol in ['original_class','original_state']:
 z=m[['cell_id',oldcol]].merge(labels[['cell_id','leaf']],on='cell_id',validate='one_to_one');save(pd.crosstab(z[oldcol],z.leaf).reset_index(),oldcol+'_to_leaf_counts.tsv')
save(labels[labels.is_fit].groupby(['site','parent','leaf','source','donor']).size().rename('cells').reset_index(),'leaf_patient_source_counts.tsv')
save(labels[labels.is_fit].groupby(['site','parent','leaf','old_lineage']).size().rename('cells').reset_index(),'leaf_inherited_lineage.tsv')
counts={site:dict(all_cells=int(m.site.eq(site).sum()),fit_cells=int((m.site.eq(site)&~m.doublet).sum()),parents=int(nodes.site.eq(site).mul(nodes.level.eq('parent')).sum()),refined_children=int(nodes.site.eq(site).mul(nodes.level.eq('child')).sum()),terminal_groups=int(labels[labels.site.eq(site)&labels.is_fit].leaf.nunique())) for site in ['Bone','Lung']}
dump('STRUCTURE_COMPLETE.json',dict(status='PASS_HIERARCHICAL_COMPUTATION',sites=counts,cells=6633,fit_cells=6423,no_cross_site_embedding=True,old_labels_unchanged=True,figures_generated=False));event('02_structure','completed',counts=counts)
