"""Monta site_src/ (entrada do MkDocs) a partir do README, dos guias, das imagens e dos slides.

O README vira uma pagina por secao `##`; a capa recebe o titulo, a introducao,
os cartoes de navegacao e o roteiro. Gera tambem mkdocs.gen.yml com o nav.
"""
import re, shutil, unicodedata
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / 'site_src'
REPO = 'https://github.com/madsondeluna/curso-md'
COLAB = 'https://colab.research.google.com/github/madsondeluna/curso-md/blob/main/MD_Martini.ipynb'

# secoes que ficam na capa ou somem do site (o menu lateral substitui o indice)
HOME = {'Roteiro da aula'}
DROP = {'Índice'}
# secoes que se juntam a anterior em vez de virar pagina propria
MERGE = {'Licença'}
# grupos do menu lateral, pelo titulo da secao
REFERENCE = {'Material legado: curso atomístico', 'Referência técnica', 'Referências'}

ICONS = ['material/atom', 'material/layers-triple', 'material/terrain', 'material/scale-balance',
         'material/molecule', 'material/water', 'material/google-circles-communities', 'material/dice-multiple']


def slug(text):
    text = unicodedata.normalize('NFKD', text).encode('ascii', 'ignore').decode()
    return re.sub(r'[^a-z0-9]+', '-', text.lower()).strip('-')[:48].rstrip('-')


def repo_links(text, prefix=''):
    """Links relativos para arquivos que nao entram no site viram links do GitHub."""
    text = re.sub(rf'\]\({prefix}(mdp_martini|legado)(/[^)]*)?\)',
                  lambda m: f']({REPO}/tree/main/{m.group(1)}{m.group(2) or ""})', text)
    return re.sub(rf'\]\({prefix}LICENSE\)', f']({REPO}/blob/main/LICENSE)', text)


def promote(body):
    """Sobe cada titulo um nivel (## vira #), sem tocar em blocos de codigo."""
    out, fence = [], False
    for line in body.split('\n'):
        if line.startswith('```'):
            fence = not fence
        if not fence and re.match(r'#{2,6} ', line):
            line = line[1:]
        out.append(line)
    return '\n'.join(out)


def label(title):
    """O numero da secao ("1. ") fica no titulo da pagina, nao no menu nem na URL."""
    return re.sub(r'^\d+\.\s*', '', title)


def split_readme(text):
    """Devolve (preambulo, [(titulo, corpo com o titulo)]) cortando em `## `."""
    parts = re.split(r'^(?=## )', text, flags=re.M)
    sections = []
    for part in parts[1:]:
        title = part.split('\n', 1)[0][3:].strip()
        sections.append((title, part.rstrip() + '\n'))
    return parts[0], sections


def home_page(preamble, roteiro, pages):
    lines = preamble.strip().split('\n')
    title = lines[0].lstrip('# ').strip()
    rest = '\n'.join(lines[1:]).strip()
    intro, _, rest = rest.partition('\n\n')
    # a tabela do notebook e substituida pelos botoes; a de autoria fica
    autoria = rest[rest.index('| Autoria'):].strip()
    cards = '\n\n'.join(
        f'-   :{ICONS[i % len(ICONS)].replace("/", "-")}:{{ .lg .middle }} [{t}]({f})'
        for i, (t, f) in enumerate(pages))
    return f'''---
hide:
  - navigation
  - toc
---

<div class="hero" markdown>

# {title}

{intro}

[Abrir no Colab]({COLAB}){{ .md-button .md-button--primary target="_blank" }}
[Começar a aula]({pages[0][1]}){{ .md-button }}
[Slides da aula (PDF)](slides/aula-dinamica-molecular.pdf){{ .md-button target="_blank" }}

</div>

## Conteúdo

<div class="grid cards" markdown>

{cards}

-   :material-file-document-outline:{{ .lg .middle }} [Passo a passo: proteína](guias/proteina.md)

-   :material-file-document-outline:{{ .lg .middle }} [Passo a passo: membrana](guias/membrana.md)

</div>

{promote(roteiro).replace('# Roteiro da aula', '## Roteiro da aula', 1)}

## Autoria

{autoria}
'''


def main():
    shutil.rmtree(OUT, ignore_errors=True)
    (OUT / 'guias').mkdir(parents=True)
    shutil.copytree(ROOT / 'imagens', OUT / 'imagens', ignore=shutil.ignore_patterns('._*'))
    shutil.copytree(ROOT / 'slides', OUT / 'slides', ignore=shutil.ignore_patterns('._*'))
    shutil.copytree(ROOT / 'site_theme', OUT, dirs_exist_ok=True, ignore=shutil.ignore_patterns('._*'))

    preamble, sections = split_readme(repo_links((ROOT / 'README.md').read_text()))

    # junta as secoes MERGE a anterior
    merged = []
    for title, body in sections:
        if title in MERGE and merged:
            merged[-1] = (merged[-1][0], merged[-1][1] + '\n' + body)
        else:
            merged.append((title, body))

    roteiro, aula, refs = '', [], []
    for title, body in merged:
        if title in DROP:
            continue
        if title in HOME:
            roteiro = body
            continue
        name = f'{slug(label(title))}.md'
        (OUT / name).write_text(promote(body))
        (refs if title in REFERENCE else aula).append((title, name))

    (OUT / 'index.md').write_text(home_page(preamble, roteiro, aula))

    for guia in sorted((ROOT / 'guias').glob('[!._]*.md')):
        text = guia.read_text().replace('[Voltar ao README](../README.md)', '[Voltar à aula](../index.md)')
        (OUT / 'guias' / guia.name).write_text(repo_links(text, r'\.\./'))

    def entries(items, indent):
        return '\n'.join(f'{indent}- "{label(t)}": {f}' for t, f in items)

    nav = f'''INHERIT: mkdocs.yml
nav:
  - Início: index.md
  - Aula:
{entries(aula, '      ')}
  - Passo a passo:
      - "Proteína": guias/proteina.md
      - "Membrana": guias/membrana.md
  - Referência:
{entries(refs, '      ')}
'''
    (ROOT / 'mkdocs.gen.yml').write_text(nav)


if __name__ == '__main__':
    main()
