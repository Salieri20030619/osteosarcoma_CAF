from pathlib import Path
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF')
for name in ['phase26_client_reanalysis_20260922/code/native_cnv.R','phase26_client_reanalysis_20260922/code/run_background.py','phase26_client_reanalysis_20260922/code/run_background_v6.py']:
    f=P/name;lines=f.read_text().splitlines();print('\nCODE',name)
    hits=[i for i,s in enumerate(lines) if any(k in s.lower() for k in ['pred_cnv','infercnv_call','state','bayes','hmm','norm.cell','annotation','ngene.chr','alteration_called','no_alteration_called'])]
    idx=sorted(set(j for i in hits for j in range(max(0,i-2),min(len(lines),i+5))))
    print('\n'.join(f'{i+1}: {lines[i]}' for i in idx))
print('\nCNV DIRECTORY STRUCTURE')
for d in [P/'phase26_client_reanalysis_20260922/cnv',P/'phase26_client_reanalysis_20260922/retry_once_20260924']:
    print(str(d),[f.name for f in d.iterdir()][:30])
    for child in d.iterdir():
        if child.is_dir() and child.name in ['infercnv','copykat','cnv']:
            print(child.name,[f.name for f in child.iterdir()][:15])
