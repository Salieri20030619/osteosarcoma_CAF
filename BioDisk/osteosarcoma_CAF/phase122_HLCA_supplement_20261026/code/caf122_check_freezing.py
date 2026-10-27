from pathlib import Path
import torch,inspect,json
from scvi.model.base import ArchesMixin,BaseModelClass
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase122_HLCA_supplement_20261026');a=torch.load(O/'objects/query_model/model.pt',map_location='cpu',weights_only=False);b=torch.load(O/'objects/query_model_extended/model.pt',map_location='cpu',weights_only=False);a=a['model_state_dict'];b=b['model_state_dict'];changes=[]
for k in a:
 if torch.is_tensor(a[k]) and a[k].shape==b[k].shape and a[k].numel() and a[k].dtype.is_floating_point:
  changes.append(dict(parameter=k,max_difference=float((a[k]-b[k]).abs().max())))
src=inspect.getsource(BaseModelClass.load);arch=inspect.getsource(ArchesMixin.load_query_data);rec=dict(load_mentions_online_update='_set_params_online_update' in src,load_query_mentions_online_update='_set_params_online_update' in arch,changed_parameters=[v for v in changes if v['max_difference']>1e-6]);(O/'resources/EXTENDED_MODEL_FREEZE_AUDIT.json').write_text(json.dumps(rec,indent=2),encoding='utf-8');print(json.dumps(rec)[:6000])
