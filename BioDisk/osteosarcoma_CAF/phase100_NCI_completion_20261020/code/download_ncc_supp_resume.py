import pathlib,subprocess,concurrent.futures,hashlib,json
p=pathlib.Path(__file__).resolve().parents[1]
src=p/'resources/ncc_supp1_springer.pdf';prefix=src.read_bytes();size=21165445
assert prefix[:4]==b'%PDF' and len(prefix)<size
parts=p/'resources/ncc_supp1_parts';parts.mkdir(exist_ok=True)
url='https://static-content.springer.com/esm/art%3A10.1186%2Fs13045-024-01598-7/MediaObjects/13045_2024_1598_MOESM1_ESM.pdf'
jobs=[(s,min(s+1048575,size-1)) for s in range(len(prefix),size,1048576)]
def fetch(item):
 s,e=item;f=parts/(str(s)+'.bin')
 r=subprocess.run(['curl','-fLsS','--connect-timeout','10','--max-time','120','--retry','2','--retry-delay','3','-r',str(s)+'-'+str(e),'-o',str(f),url],capture_output=True,text=True)
 assert r.returncode==0 and f.stat().st_size==e-s+1,(s,r.returncode,r.stderr)
 return f
with concurrent.futures.ThreadPoolExecutor(6) as ex:files=list(ex.map(fetch,jobs))
out=p/'resources/ncc_supp1_verified.pdf'
with out.open('wb') as w:
 w.write(prefix)
 for f in files:w.write(f.read_bytes())
assert out.stat().st_size==size and hashlib.md5(out.read_bytes()).hexdigest()=='734d6b81da4646c2ecd11ced4071f554'
r=subprocess.run(['pdftotext','-layout',str(out),str(p/'resources/ncc_supp1_text.txt')],capture_output=True,text=True)
assert r.returncode==0,r.stderr
(p/'resources/NCC_SUPP_PDF_COMPLETE.json').write_text(json.dumps(dict(url=url,bytes=size,md5='734d6b81da4646c2ecd11ced4071f554',matches_publisher_etag=True,pdf_parse_exit=0)),encoding='utf-8')
print('NCC_SUPP_VERIFIED',flush=True)
print((p/'resources/ncc_supp1_text.txt').read_text(encoding='utf-8')[:16000],flush=True)
