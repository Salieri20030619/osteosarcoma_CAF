from pathlib import Path
import os,json,threading,time,datetime,math
def available_gib():
 return int(next(x for x in Path('/proc/meminfo').read_text().splitlines() if x.startswith('MemAvailable:')).split()[1])/2**20
def tree_pss(pid):
 seen=set();queue=[pid];total=0
 while queue:
  p=queue.pop()
  if p in seen:continue
  seen.add(p)
  try:
   lines=(Path('/proc')/str(p)/'smaps_rollup').read_text().splitlines()
   total+=sum(int(x.split()[1]) for x in lines if x.startswith('Pss:'))/2**20
   queue.extend(map(int,(Path('/proc')/str(p)/'task'/str(p)/'children').read_text().split()))
  except (OSError,ValueError):pass
 return total
def estimate_gib(qa):
 n=int(qa['cells']);g=int(qa['genes'])
 # Admission estimate only, never a biological inclusion rule or an absolute cap.
 # Dense matrix copies + pairwise workspace + per-cell processing + interpreter overhead.
 return math.ceil(max(12,8+4*n*g*8/2**30+3*n*n*8/2**30+.001*n))
class MemoryGate:
 def __init__(self,out,reserve=100,available_fn=available_gib,pss_fn=tree_pss):
  self.out=Path(out);self.reserve=reserve;self.available=available_fn;self.pss=pss_fn
  self.cond=threading.Condition();self.tokens={};self.waiting={}
 def can_admit(self,available,uncommitted,budget):return available-uncommitted-budget>=self.reserve
 def status(self,key=None,budget=None):
  available=self.available()
  future=sum(max(0,v['budget_GiB']-(self.pss(v['pid']) if v['pid'] else 0)) for v in self.tokens.values())
  x={'at':datetime.datetime.now().astimezone().isoformat(),'reserve_available_GiB':self.reserve,'available_GiB':available,'uncommitted_reserved_GiB':future,'allocations':self.tokens,'waiting':self.waiting,'memory_pressure':available<self.reserve,'policy':'Do not launch if estimated additional allocations would reduce available memory below reserve. Existing tasks and other users can change availability. Estimates are not hard cgroup isolation.'}
  tmp=self.out.with_suffix('.tmp');tmp.write_text(json.dumps(x,indent=2));tmp.replace(self.out)
  return available,future
 def admit(self,key,budget):
  with self.cond:
   self.waiting[key]=budget
   while True:
    available,future=self.status()
    if self.can_admit(available,future,budget):
     self.tokens[key]={'budget_GiB':budget,'pid':None};self.waiting.pop(key);self.status();return
    self.cond.wait(30)
 def attach(self,key,pid):
  with self.cond:self.tokens[key]['pid']=pid;self.status()
 def release(self,key):
  with self.cond:self.tokens.pop(key,None);self.status();self.cond.notify_all()
