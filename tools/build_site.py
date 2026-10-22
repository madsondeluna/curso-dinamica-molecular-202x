"""Monta site_src/ (entrada do MkDocs) a partir do README, dos guias e das imagens."""
import re, shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / 'site_src'
REPO = 'https://github.com/madsondeluna/curso-md'

# links relativos para arquivos que nao entram no site viram links do GitHub
def repo_links(text, prefix=''):
    text = re.sub(rf'\]\({prefix}(mdp_martini|legado)(/[^)]*)?\)',
                  lambda m: f']({REPO}/tree/main/{m.group(1)}{m.group(2) or ""})', text)
    return re.sub(rf'\]\({prefix}LICENSE\)', f']({REPO}/blob/main/LICENSE)', text)

if OUT.exists():
    shutil.rmtree(OUT, ignore_errors=True)
(OUT / 'guias').mkdir(parents=True)
shutil.copytree(ROOT / 'imagens', OUT / 'imagens', ignore=shutil.ignore_patterns('._*'))
shutil.copytree(ROOT / 'site_theme', OUT, dirs_exist_ok=True, ignore=shutil.ignore_patterns('._*'))

(OUT / 'index.md').write_text(repo_links((ROOT / 'README.md').read_text()))
for guia in sorted((ROOT / 'guias').glob('[!._]*.md')):
    text = guia.read_text().replace('[Voltar ao README](../README.md)', '[Voltar à aula](../index.md)')
    (OUT / 'guias' / guia.name).write_text(repo_links(text, r'\.\./'))
