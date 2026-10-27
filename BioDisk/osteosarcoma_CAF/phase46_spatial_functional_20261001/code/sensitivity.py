from common import *
z=pd.read_csv(O/'tables/all_spot_RNA_and_components.tsv.gz',sep='\t');pr=json.loads((O/'resources/FUNCTIONAL_SIGNATURES_FROZEN.json').read_text());names=list(pr['programs']);design={'library_xy':['library_UMI','x','y'],'cancer_library_xy':['cancer','library_UMI','x','y'],'fibroblast_library_xy':['fibroblast','library_UMI','x','y']}
dump(O/'resources/SPATIAL_MODEL_SENSITIVITY_FROZEN.json',{'at':datetime.datetime.now().astimezone().isoformat(),'reason':'primary broad null and selective3gene negative only after full adjustment; assess dependence on cancer/fibroblast covariates after seeing main results','models':design,'all15programs_included':True,'family':'separate by model; unchanged primary results','post_primary_sensitivity':True})
rows=[]
for sm,t in z.groupby('sample'):
 for model,cols in design.items():
  c=t[cols].to_numpy(float)
  for j,col in enumerate(cols):
   if col=='library_UMI':c[:,j]=np.log1p(c[:,j])
  for name in names:rows.append(dict(sample=sm,patient_id=t.patient_id.iloc[0],site_group=t.site_group.iloc[0],event_month=t.event_month.iloc[0],program=name,model=model,effect=corr(t[name],t.CD8_RNA,c)))
f=pd.DataFrame(rows);f.to_csv(O/'tables/spatial_sensitivity_sample_effects.tsv',sep='\t',index=False);ev=f.groupby(['patient_id','site_group','event_month','program','model'],dropna=False).effect.mean().reset_index();si=ev.groupby(['patient_id','site_group','program','model']).effect.mean().reset_index();pa=si.groupby(['patient_id','program','model']).effect.mean().reset_index();pa.to_csv(O/'tables/spatial_sensitivity_patient_effects.tsv',sep='\t',index=False);s=test_table(pa,['program','model'],'effect',lambda r:'spatial_postprimary_model_'+r.model);s.to_csv(O/'tables/spatial_model_sensitivity_statistics.tsv',sep='\t',index=False);(O/'SENSITIVITY_COMPLETE').write_text('SUCCESS');print(s[s.program.isin(['F3_SELECTED_REMODELING','NABA_CORE_MATRISOME'])].to_string(index=False))
