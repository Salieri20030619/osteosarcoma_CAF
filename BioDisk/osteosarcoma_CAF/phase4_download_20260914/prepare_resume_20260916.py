from pathlib import Path
p=Path('download_validate_http1.py')
s=p.read_text()
s=s.replace('import numpy as np','import numpy as np\nfrom threadpoolctl import threadpool_limits\nthreadpool_limits(limits=2)')
a=s.index("  part=p.with_name(p.name+'.geo.part')")
b=s.index('  assert part.stat().st_size>0',a)
s=s[:a]+'''  part=p.with_name(p.name+'.resume20260916.part')
  transport=url
  print('DOWNLOAD_START',p.name,'resumable_official_https',flush=True)
  for attempt in range(1,81):
   before=part.stat().st_size if part.exists() else 0
   result=subprocess.run(['curl','--http1.1','--continue-at','-','--max-time','240','--fail','--location','--connect-timeout','30','--speed-limit','1024','--speed-time','90','--output',str(part),'--silent','--show-error',transport])
   after=part.stat().st_size if part.exists() else 0
   print('TRANSFER_ATTEMPT',p.name,attempt,'exit',result.returncode,'bytes',after,'added',after-before,flush=True)
   if result.returncode==0:break
   if result.returncode==33:raise RuntimeError('RANGE_NOT_SUPPORTED '+p.name)
   time.sleep(5 if after>before else 20)
  else:raise RuntimeError('RETRY_LIMIT '+p.name)
'''+s[b:]
s=s.replace('download_validate_http1.py','download_validate_resume_20260916.py')
Path('download_validate_resume_20260916.py').write_text(s)
