from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote
class Parser(HTMLParser):
 def __init__(self): super().__init__();self.images=[];self.links=[];self.text=[]
 def handle_starttag(self,tag,attrs):
  a=dict(attrs)
  if tag=='img' and 'src'in a:self.images.append(a['src'])
  if tag=='a' and 'href'in a:self.links.append(a['href'])
 def handle_data(self,data):self.text.append(data)
p=Parser();p.feed(Path('index.html').read_text());text=''.join(p.text)
assert 'Output created: index.html' in Path('revision_checks/equalto_2026-10-01/render_final.log').read_text()
assert 'glmmTMB_1.1.15.2' in text
assert 'The five equalto() models' in text
assert text.count('vi_named <- setNames')==5
assert 'Users should not rely on the modelling software' in text
assert 'revision_checks/equalto_2026-10-01/REPORT.md' in p.links
for src in p.images:
 if not src.startswith(('http:','https:','data:')): assert Path(unquote(src)).is_file(),src
print('FINAL_HTML_VERSION_ALIGNMENT_AND_LOCAL_IMAGE_REFERENCES_PASSED')
