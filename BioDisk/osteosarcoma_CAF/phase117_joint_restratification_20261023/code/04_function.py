from caf117_function_utils import *
assert (O/'CORRESPONDENCE_COMPLETE.json').exists();event('04_function','started');X,m,genes=load();labs=read('joint_cell_labels.tsv.gz').set_index('cell_id');m['joint_cluster']=m.cell_id.map(labs.joint_cluster);v,names,family=endpoints(m,X,genes)
df=evaluate(m,v,names,family,'joint_cluster',['all','singlet','singlet_min5'],'joint_function')
summary=df[df.arm.eq('singlet')&df.family.eq('functional_program')];key=['fibrillar_I_III_V','collagen_biosynthesis','collagen_crosslinking','myCAF_customer','F2C3_collagen_support','Osteogenic_customer','iCAF_customer']
save(summary[summary.endpoint.isin(key)],'joint_key_function_statistics.tsv');save(df[df.nominal&df.mean_effect.gt(0)&df.arm.eq('singlet')],'joint_nominal_positive_results.tsv')
# Descriptive source-stratified functional profiles, rather than cross-study pooled inference.
profile=pd.concat([m[['source','donor','joint_cluster','doublet']],pd.DataFrame(v,columns=names)],axis=1);save(profile[~profile.doublet].groupby(['source','donor','joint_cluster'])[names].mean().reset_index(),'joint_patient_function_profiles.tsv')
dump('FUNCTION_COMPLETE.json',dict(status='PASS',endpoints=len(names),families=sorted(set(family.values())),statistic_rows=len(df),valid_tests=int(df.P.notna().sum()),nominal_positive=int((df.arm.eq('singlet')&df.nominal&df.mean_effect.gt(0)).sum()),source_stratified=True,expression_not_integrated=True,scope='Conditional functional evidence; not independent subtype validation'))
event('04_function','completed',valid_tests=int(df.P.notna().sum()))
