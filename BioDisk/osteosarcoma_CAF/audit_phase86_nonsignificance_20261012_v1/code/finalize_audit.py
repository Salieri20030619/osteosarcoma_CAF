
from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,zipfile,datetime
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'audit_phase86_nonsignificance_20261012_v1';B=P/'phase86_functional_bone_lung_20261012'
def r(n):return pd.read_csv(O/'tables'/n,sep='\t')
summary=json.loads((O/'DIAGNOSTIC_SUMMARY.json').read_text())
for n in ['exit_code','sample_weight.exit_code','sensitivity_QA.exit_code']:assert (O/'logs'/n).read_text().strip()=='0'
hashes=json.loads((O/'INPUT_SHA256.json').read_text());assert all(hashlib.sha256(Path(k).read_bytes()).hexdigest()==v for k,v in hashes.items())
c=pd.read_csv(B/'tables/cell_scores_and_original_labels.tsv.gz',sep='\t');pat=pd.read_csv(B/'tables/patient_site_programs.tsv',sep='\t');th=pd.read_csv(B/'tables/state_thresholds.tsv',sep='\t').set_index(['scope','dataset','matrix','program']).threshold.to_dict();errors=[]
for v in pat.itertuples():
 t=c[(c.dataset==v.dataset)&(c.donor_record==v.donor_record)&(c.strict_site==v.site)]
 if v.scope=='Fibro_reviewed':t=t[t.review_lineage_code=='Fibro']
 score=t[v.matrix+'__'+v.program];valid=score.notna();frac=(score[valid]>th[v.scope,v.dataset,v.matrix,v.program]).mean();errors.append(abs(frac-v.high_fraction));assert errors[-1]<1e-10
main=r('03_primary_effects_and_resolution.tsv');s17=main[(main.scope=='Fibro_reviewed')&(main.dataset=='SCPCP000017')];assert s17.LOPO_direction_retained.eq(True).sum()==5
text='''CAF 骨—肺功能分析未达名义显著的原因调查
审计对象：Phase86（2026-10-12）。本次仅作原因调查和诊断性敏感性分析，未修改原主分析、细胞标签及项目状态。

结论
原分析500项有效功能检验均无P<0.05，复算无误。其原因不是FDR过严：一半有效检验受小样本精确检验分辨率限制；另一半还受到肺端患者少、患者波动及宽功能集合内相反变化的影响。细胞身份限定和临床取材差异削弱跨来源验证。
新发现：患者内多样本的权重会影响结果。样本等权后收缩程序出现4条相互重叠的名义显著结果，均非FDR显著；但该变化主要由一个仅1细胞的样本增权驱动。应保留为权重敏感的探索性线索，不能替换原结果成为稳定阳性结论。

1. 500项有效检验中，250项数学上无法达到P<0.05
800项登记中，300项无足够生物学重复、没有P值；500项有P，其中250项的最小可达P不低于0.05，另250项有达到P<0.05的数学可能。原始最小P=0.088911。
S17 Fibro主非配对为骨11、肺3名独立患者：364种分配，现有绝对均值差双侧检验的设计P下界为1/364=0.002747。
S17 Fibro主配对只有2对：4种符号组合，双侧最小P=2/4=0.5。
S23 Fibro主非配对只有骨2、肺2人：6种分配，双侧最小P=2/6=0.333333。
S23完整候选池骨4、肺3人时P下界1/35=0.028571；每患者部位至少5细胞时变为骨2、肺3人，P下界1/10=0.1。
这些下界针对本项目现有双侧定义，不把不等组设计一律写为2/组合数。所有原始/去污染、身份层、覆盖阈值和高评分比例比较互相重叠，500条不是500次独立研究，不能按“应有25个阳性”推断算法异常。

2. S17胶原/基质的净效应弱，患者一致性不足
S17主Fibro原始平均评分14项中，只有5项在逐患者剔除后方向始终保留，另9项存在方向翻转。
胶原合成肺减骨=-0.041807，描述性标准化效应=-0.188，P=0.774725，FDR=0.901099。
胶原装配=-0.029473，标准化效应=-0.138，P=0.829670，FDR=0.901099。
ECM组织=+0.030537，标准化效应=+0.303，P=0.637363，FDR=0.901099。
效应单位是平均log1p(CP10k)评分差，不是log2倍数变化；标准化效应为组间差除以合并患者标准差，仅用于描述。
IL6-JAK-STAT3和炎症方向更稳定、标准化效应分别1.137和1.056，但P分别0.107143和0.129121，FDR均0.901099，仍未达到筛选标准。P>0.05不得称为名义显著，也不等于没有变化。
S23的14项主Fibro方向只有5项与S17同向。S23细胞/患者极少，方向比较仅作探索性描述，不是可靠的阴性验证。

3. 同一功能集合中相反基因变化被均值抵消
在原公式下逐基因分解S17主Fibro效应，胶原装配61基因中27项正贡献、32项负贡献，其余为零，算术抵消比例83.4%；胶原合成76.3%、交联74.6%、降解87.0%、胶原组成88.9%。
抵消比例=1-|贡献总和|/绝对贡献总和，是描述性拆解，不表示单基因已统计显著。这支持“整套程序未一致上升”，不能扩展为“所有局部功能均未变化”；不应事后删除反向基因来制造显著评分。
S17主Fibro的91个程序配对中，20对患者部位评分的Spearman绝对相关达到0.8，亦说明各功能结果不独立。

4. S23的Fibro限定仅保留1.42%候选细胞
S17完整候选池3732细胞，Fibro审核层1368细胞。
S23完整候选池1834细胞（骨311、肺1523），Fibro层仅26细胞（骨7、肺19），约1.42%，骨肺各2人。每组各有一个样本仅2个Fibro细胞。
其余多为Mesench、Osteo等历史审核标签。本审计不重注释、不据标签直接宣判其真实谱系；该差异说明跨来源验证对象并不充分可比。
两个来源本轮候选池没有再因双联体阳性排除损失细胞，阴性不能归咎于本轮双联体剔除。
完整候选池的功能敏感性同样没有P<0.05，所以也不能把所有阴性仅归因于Fibro限定。

5. 样本权重可产生收缩名义信号，但受单细胞样本驱动
原流程在同一患者同一部位内合并细胞取均值：患者之间等权，患者内部不同样本按细胞数加权，并未把每个样本当成独立患者。
本次对全部800项登记统一检查样本等权：保留患者名单、细胞、程序、原高评分阈值和检验，只改为患者内部先按样本求均值，再样本等权。500项可检验记录出现4条名义显著，均属于S17 Fibro收缩程序，均未FDR显著。
原始评分：效应+0.074930、P=0.269231，变为+0.172546、P=0.043956、FDR=0.615385。
去污染评分：效应+0.084517、P=0.206044，变为+0.199093、P=0.024725、FDR=0.271978。
去污染和FDR多重校正分别记录。这4条是重叠分析设置，不能称为4次独立验证。
关键肺患者CCG1106_007有两份复发样本，分别1个和14个Fibro细胞，原始收缩评分0.916074与0.242230。样本等权使单细胞样本从1/15增至该患者50%权重。
仅恢复该患者原权重，其余患者仍样本等权，原始P回到0.269231，去污染P=0.211538；剔除该患者的影响诊断也未显著。此影响诊断用于解释敏感性，不作为额外筛选家族或新的主检验。
另有4项S17胶原相关方向随样本权重从负转正。故样本聚合规则影响解释，但“等权更显著”不是其更正确的依据；需先按取材性质、时点和最低细胞可靠性确定规则。

6. 临床和技术混杂有证据，但不能直接认定其导致全部阴性
S17主肺组仅3人，其中2人为转移切除、1人为复发；骨组混有初诊、局部控制手术及治疗后切除。两个骨患者还包含跨病程取材。S17治疗字段缺失，无法据此充分调整治疗影响。
S23主Fibro骨肺的病程与治疗分布也不同。临床患者可跨病程出现，各类别患者计数不能直接相加。详见实际主分析名单表21—22。
S17主Fibro评分与检测基因数的描述性Spearman相关：TGF-beta 0.767、矿化0.593、炎症0.538；胶原装配0.044、合成-0.036。技术复杂度对不同程序影响不同，不能统一解释全部阴性。这些是患者部位层描述，未作独立患者推断。

7. 未发现P值、细胞平均分或高评分比例的实现错误
500个原始有效P全部独立枚举复现，最大误差1.11e-16。患者评分从细胞数据重建最大误差2.22e-16；原高评分比例逐条重建一致。未发现同一患者重复进入独立两组、均值计算或细胞索引错误。
两个来源14程序的样本特征测量覆盖为100%，不意味着每细胞均表达；基因缺测补零不是本轮阴性的解释。
S17主Fibro11项有原始和去污染版本的程序效应方向全部一致、均未P<0.05，因此不是去污染把原先明确功能阳性全部抹掉。

8. 不能用本轮阴性否定此前14候选，整个项目并非“都不显著”
此前14候选来自完整候选池单基因骨肺比较，两来源均名义显著、均未FDR显著。本轮检验的是预设功能评分或高评分比例，主细胞层也变为Fibro审核层。
14候选中只有BEST1直接被这14程序之一收录（炎症）。这不证明其他候选与这些功能无关，只说明本轮没有逐一验证14基因所代表的功能状态。
Phase86患者内标志—功能关联另有39149条有效检验，其中2854条名义显著，包含重复设置和关联条目。因此“全部不显著”只能限定为原骨肺功能比较。

后续优先顺序
保留14候选及其原始证据；先明确完整候选池单基因与Fibro内功能状态两个目标。
先复核S23上游身份标签和成纤维细胞覆盖，避免将26个细胞承担强重复验证；本次未更改标签或扩大细胞池。
核定同患者多样本的生物学性质、时点及细胞量要求，再预设聚合口径。收缩保留为权重敏感探索性线索，不以单细胞增权追求阳性。
优先增加肺端独立患者、配对样本及临床信息；增加同患者细胞数不能替代生物学样本量。具体样本量另按目标效应设计，不从本次观察P值倒推。
不以更换阈值、只保留同向基因、合并细胞作为独立重复的方式获取显著。

方法、材料与复核位置
审计输入：Phase86原统计、患者部位评分、细胞评分与标签、基因覆盖；Phase63冻结机制矩阵；Phase82原14候选表。全部原输入保持不变。
新增诊断：精确检验分辨率、患者逐一影响、原评分逐基因贡献、同患者样本等权敏感性。FDR分别在来源/细胞层/矩阵/覆盖阈值/设计/端点家族内计算；不与原检验混池。
新增4条名义敏感性结果用SciPy独立精确置换复核；绝对均值差作为统计量，使用greater尾概率复现原双侧定义。官方方法依据：https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.permutation_test.html
表01—04保留全部原800项P/FDR和检验分辨率；表05—16为纳入、覆盖、评分重建、技术协变和基因抵消；表17—24为完整敏感性统计、临床取材与单患者影响。脚本、日志和输入SHA256均归档。
'''
text+='\n附表A：原始矩阵、主覆盖阈值、非配对平均评分完整结果\n'+main[['scope','dataset','program','bone_N','lung_N','mean_lung_minus_bone','P','FDR','LOPO_direction_retained']].to_string(index=False)+'\n'
w=r('19_all800_equal_sample_weight_sensitivity.tsv');text+='\n附表B：新增等样本权重敏感性全部名义显著条目\n'+w[w.nominal_sensitivity][['scope','dataset','matrix','program','minimum_cells','design','endpoint','equal_sample_effect','equal_sample_P','equal_sample_FDR']].to_string(index=False)+'\n'
(O/'AUDIT_REPORT_ZH.txt').write_text(text,encoding='utf-8-sig')
qa={'status':'PASS','completed_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'original_P_reproduced':500,'original_patient_score_max_error':summary['aggregation_max_error'],'original_high_fraction_rows_checked':len(errors),'high_fraction_max_error':max(errors),'all_inputs_unchanged':True,'new_nominal_sensitivity_rows':4,'sensitivity_QA':json.loads((O/'SENSITIVITY_QA.json').read_text()),'scope':'Independent diagnosis, no primary pipeline replacement'}
(O/'FINAL_ACCEPTANCE.json').write_text(json.dumps(qa,indent=2))
files=[f for f in sorted(O.rglob('*')) if f.is_file() and f.suffix!='.zip' and f.name!='FILE_MANIFEST.json']
manifest=[dict(path=str(f.relative_to(O)),bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest()) for f in files];(O/'FILE_MANIFEST.json').write_text(json.dumps(manifest,indent=2))
archive=O/'CAF_phase86_nonsignificance_audit_20261012.zip'
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
 for f in files+[O/'FILE_MANIFEST.json']:z.write(f,str(f.relative_to(O)))
with zipfile.ZipFile(archive) as z:assert z.testzip() is None
print(json.dumps({'archive':str(archive),'bytes':archive.stat().st_size,'sha256':hashlib.sha256(archive.read_bytes()).hexdigest(),'files':len(files)+1,'QA':qa},ensure_ascii=False))
