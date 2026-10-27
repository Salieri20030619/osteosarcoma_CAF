from pathlib import Path
import pandas as pd,json,hashlib,datetime
P=Path('C:/codex_test');O=P/'CAF_external_delivery_20260921'
g=pd.read_csv(P/'CAF_evidence_20260921/tables/gene_candidate_evidence.tsv',sep='\t')
signatures={}
for cl in [7,9]:
 d=g[(g.layer=='supported')&(g.cluster==cl)&(g.evidence_level=='G2_exploratory_recurrent')].sort_values(['median_effect','symbol'],ascending=[False,True])
 signatures['C'+str(cl)]=d.symbol.head(15).tolist()
assert all(len(v)>=5 for v in signatures.values())
protocol={'frozen_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'candidate_source':'Phase18 supported-layer C7/C9 recurrent positive marker genes; up to top15 by within-donor median effect; no outcome-based selection',
 'signatures':signatures,'primary_external_cohort':'GSE21257','primary_endpoint':'original any-metastasis within5y group, including baseline cases as in source; not lung-specific',
 'score':'mean of per-gene standardized author-normalized expression across cohort; standardize composite to1SD; >80% signature coverage required; probe median collapse',
 'primary_statistics':'unadjusted logistic OR per1SD with95%WaldCI and two-sidedP; BH across2 candidate signatures for GSE21257 primary endpoint',
 'sensitivities':['GSE21257 remove baseline metastasis cases','GSE21257 adjust age and sex','rank-based within-sample percentile score replacing cohort-gene z scores'],
 'new_replication_cohort':'GSE39055; use official clinical recurrence-free and overall survival fields if available; do not relabel recurrence as metastasis',
 'replication_statistics':'Cox PH continuous score HR per1SD, Efron ties; separate BH family2 signatures per endpoint; survival time<=0 excluded with count; no optimized cutpoints',
 'scope':'external tissue-level outcome association; does not validate CAF identity, single-cell subtype transfer or protein sorting purity',
 'prior_use':'GSE21257 previously used for other frozen generic programs in Phase8, so this is new-candidate retrospective external cohort testing, not a prospectively untouched cohort. GSE39055 not used in discovery.',
 'source_checksum':hashlib.sha256((P/'CAF_evidence_20260921/tables/gene_candidate_evidence.tsv').read_bytes()).hexdigest()}
f=O/'FROZEN_PROTOCOL.json'
if f.exists():assert json.loads(f.read_text())['signatures']==signatures
else:f.write_text(json.dumps(protocol,ensure_ascii=False,indent=2),encoding='utf-8')
pd.DataFrame([{'signature':k,'gene':x,'rank':i+1} for k,v in signatures.items() for i,x in enumerate(v)]).to_csv(O/'tables/frozen_signatures.tsv',sep='\t',index=False)
print(json.dumps(signatures),flush=True)
