from caf117_common import *
import re,shutil,zipfile
event('07_report_review','started')
path=O/'tables/joint_cluster_origins_and_stability.tsv';before=pd.read_csv(path,sep='\t');after=before.copy();mask=after.old_route.eq('dominant_old_state')&after.leading_old_fraction_among_old.lt(.7)
after.loc[mask,'old_route']='mixed_old_states'
assert before.drop(columns='old_route').equals(after.drop(columns='old_route'))
backup=O/'resources/report_review_before';backup.mkdir(exist_ok=True)
for p in [path,O/'delivery/CAF_joint_restratification_report.html',O/'delivery/CAF_joint_full_tables.html',O/'FINAL_STATUS.json']:
 assert not (backup/p.name).exists();shutil.copy2(p,backup/p.name)
save(after,'joint_cluster_origins_and_stability.tsv')
def tab(d):return '<div class="scroll">'+d.to_html(index=False,border=0,escape=True,float_format=lambda v:f'{v:.6g}',na_rep='不可评估')+'</div>'
routes=read('old33_joint_correspondence.tsv');n=int((~routes.count_patient_route_agree).sum());c3=routes[routes.old_state.eq('F2_C3')].iloc[0];cl=c3.leading_joint_cluster
focus=after[after.joint_cluster.eq(cl)].iloc[0]
p=O/'delivery/CAF_joint_restratification_report.html';s=p.read_text(encoding='utf-8');assert tab(before) in s;s=s.replace(tab(before),tab(after));a=before[before.joint_cluster.eq(cl)];b=after[after.joint_cluster.eq(cl)];assert tab(a) in s;s=s.replace(tab(a),tab(b))
anchor='<h2>3. 原33亚群的真实细胞去向</h2>';assert anchor in s;s=s.replace(anchor,anchor+f'<p><b>患者加权的重要差异：</b>{n}/33个原亚群按细胞数与患者等权得到的首要去向不同。下表保留两种口径，不能只凭细胞多数即宣称跨患者稳定复现。</p>')
anchor='<h2>4. F2-C3及基质相关亚群的重点对应</h2>';assert anchor in s;s=s.replace(anchor,anchor+f'<p><b>{cl}是包含F2-C3的较宽共同群。</b>该群有{int(focus.cells)}个主拟合细胞，其中旧细胞{int(focus.original_cells)}个、NCC新增{int(focus.NCC_cells)}个、NCI新增{int(focus.NCI_cells)}个。F2-C3虽为其最大的单个旧成分，也仅占旧细胞的{focus.leading_old_fraction_among_old:.1%}；其余来自多个原亚群。因此不能将整个{cl}重命名为F2-C3。未达到70%旧标签主导且不满足两旧群各20%汇合规则者，表中明确标为mixed_old_states。</p>')
anchor='<h2>5. 留出患者的发现与验证</h2>';assert anchor in s;s=s.replace(anchor,anchor+'<p><b>本轮留出结果仅提供有限方向复核：</b>NCC留出1位患者，不作推断检验；NCI留出3位患者，J07交联和收缩面板均3/3同向，纤维胶原及原F2-C3签名2/3同向，未达名义显著。不能把全数据中的显著性写成已经通过独立患者验证。</p>')
p.write_text(s,encoding='utf-8')
p=O/'delivery/CAF_joint_full_tables.html';s=p.read_text(encoding='utf-8');match=re.search(r'const datasets=(.*?);const sel=',s,re.S);assert match;data=json.loads(match.group(1));data['共同新群来源和稳定性']=json.loads(after.to_json(orient='records',double_precision=15));s=s[:match.start(1)]+json.dumps(data,ensure_ascii=False,allow_nan=False)+s[match.end(1):];p.write_text(s,encoding='utf-8')
amend=dict(id='P117-01',time=now(),scope='Interpretation field and report clarity only',evidence='Some new groups had no old state >=20pct, yet fallback enum said dominant_old_state; J07 largest old state is only18pct',correction='Preserve existing >=2 old states each20pct convergence rule; fallback is mixed_old_states unless leading old state>=70pct',changed_groups=after.loc[mask,'joint_cluster'].tolist(),numeric_statistics_changed=False,analysis_clusters_changed=False,original_plan_unchanged=True,report_notes=['11 count/patient primary route disagreements made prominent','J07 mixed old composition explicitly stated','Limited heldoutN1/N3 direction evidence explicitly stated'])
dump('resources/AMENDMENT_01_REPORT_INTERPRETATION.json',amend)
archive=O/'delivery/CAF_joint_restratification_Phase117.zip'
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
 for f in [O/'delivery/CAF_joint_restratification_report.html',O/'delivery/CAF_joint_full_tables.html']:z.write(f,f.name)
with zipfile.ZipFile(archive) as z:assert z.testzip() is None
status=json.loads((O/'FINAL_STATUS.json').read_text());status.update(archive_sha256=sha(archive),archive_bytes=archive.stat().st_size,report_revision='interpretation_review_v2',reviewed_at=now());dump('FINAL_STATUS.json',status);dump('receipts/07_report_review.json',dict(exit_code=0,status='PASS',code_sha256=sha(Path(__file__)),amendment=amend,archive_sha256=status['archive_sha256']))
event('07_report_review','completed',archive_sha256=status['archive_sha256']);print(json.dumps(status))
