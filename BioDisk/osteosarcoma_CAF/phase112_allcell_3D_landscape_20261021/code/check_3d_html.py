from pathlib import Path
import sys,json
sys.path.insert(0,'C:/Users/31937/CAF_paper_figures_20261021/browserqa_deps')
from playwright.sync_api import sync_playwright
P=Path('C:/Users/31937/CAF_paper_figures_20261021/allcell_3d')
O=Path('C:/Users/31937/Downloads/CAF_allcell_3D_Fig1b_1d_20261021')
expected={'SCPCP000017':{'All cells':167180,'Bone':104767,'Lung':41767,'Unknown / mixed site':20646},'SCPCP000023':{'All cells':178039,'Bone':90351,'Lung':82549,'Unknown / mixed site':5139}}
checks=[]
with sync_playwright() as pw:
 browser=pw.chromium.launch(channel='msedge',headless=True,args=['--use-gl=angle','--use-angle=swiftshader','--enable-unsafe-swiftshader'])
 for src,counts in expected.items():
  page=browser.new_page(viewport={'width':1500,'height':1050},device_scale_factor=1);errors=[];page.on('pageerror',lambda e:errors.append(str(e)))
  page.goto((O/f'Fig1b_{src}_interactive_3D.html').as_uri(),wait_until='load',timeout=60000)
  page.wait_for_function("document.querySelector('.js-plotly-plot')?._fullData?.length > 0",timeout=60000)
  def count():return page.evaluate("() => document.querySelector('.js-plotly-plot')._fullData.filter(t=>t.visible!==false&&t.visible!=='legendonly').reduce((n,t)=>n+t.x.length,0)")
  assert count()==counts['All cells'],(src,count())
  page.screenshot(path=str(P/'previews'/f'{src}_interactive_all.png'))
  for label in ['Bone','Lung','Unknown / mixed site','All cells']:
   page.get_by_text(label,exact=True).click()
   page.wait_for_function('(n) => document.querySelector(\'.js-plotly-plot\')._fullData.filter(t=>t.visible!==false&&t.visible!==\'legendonly\').reduce((v,t)=>v+t.x.length,0)===n',arg=counts[label],timeout=60000)
   checks.append(dict(source=src,filter=label,visible_cells=count(),expected=counts[label],status='PASS'))
  page.evaluate("async () => {let g=document.querySelector('.js-plotly-plot');await Plotly.relayout(g,{'scene.camera.eye':{x:-1.4,y:1.6,z:.9}})}")
  page.screenshot(path=str(P/'previews'/f'{src}_interactive_rotated.png'))
  assert not errors,errors
  checks.append(dict(source=src,rotation='PASS',page_errors=errors));page.close()
 browser.close()
(P/'resources/INTERACTIVE_BROWSER_QA.json').write_text(json.dumps(dict(status='PASS',browser='headless Microsoft Edge',checks=checks),indent=2),encoding='utf-8')
print(json.dumps(checks))
