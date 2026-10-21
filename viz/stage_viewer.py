"""Cartoes de etapa em Pure Design para o notebook Martini.

Cada saida do Colab e um iframe isolado, entao cada cartao carrega o seu
proprio CSS (tokens, patterns, light, stage) e o light.js inline.

Uso:
    from stage_viewer import stage, compare_resolutions, BEADS
    stage(part='A', step=3, title='Solvatacao', structure='ubq_solv.gro',
          text=['...'], facts={'beads': 2304}, compare={...})
"""
from __future__ import annotations

import html
from pathlib import Path

import py3Dmol
from IPython.display import HTML, display

from pure import palette

HERE = Path(__file__).resolve().parent
PURE = HERE / 'pure'

# Series de dados: as oito cores de --chart-1..8, sempre na mesma ordem
CHART = palette.CATEGORICAL
INK = palette.tokens('light')['secondary']

FONTS = ('<link rel="preconnect" href="https://fonts.googleapis.com">'
         '<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>'
         '<link href="https://fonts.googleapis.com/css2?family=Archivo:wdth,wght@125,300..400'
         '&family=Public+Sans:wght@300;400;500;600&family=Spline+Sans+Mono:wght@400;500'
         '&display=swap" rel="stylesheet">')

# Modo segue o sistema: escuro vira .dark, claro fica em :root
MODE_JS = ("<script>if (window.matchMedia('(prefers-color-scheme: dark)').matches)"
           "{document.documentElement.classList.add('dark');}</script>")

_ASSETS = None


def _assets():
    global _ASSETS
    if _ASSETS is None:
        css = '\n'.join((PURE / f).read_text() for f in ('tokens.css', 'patterns.css', 'light.css'))
        css += '\n' + (HERE / 'stage.css').read_text()
        # um comentario do light.js cita a tag de fechamento de script e
        # encerraria o bloco inline antes da hora
        js = (PURE / 'light.js').read_text().replace('</script', '<\\/script')
        _ASSETS = (css, js)
    return _ASSETS


# Beads Martini: nome na legenda, slot da serie, selecao py3Dmol, raio (A)
_TAILS = [f'C{i}{s}' for i in range(1, 6) for s in 'AB'] + [f'D{i}{s}' for i in range(1, 6) for s in 'AB']
BEADS = {
    'BB':  ('Backbone (BB)',          1, {'atom': 'BB'}, 2.0),
    'SC':  ('Cadeia lateral (SC)',    2, {'atom': ['SC1', 'SC2', 'SC3', 'SC4', 'SC5']}, 1.7),
    'NC3': ('Colina (NC3)',           3, {'atom': 'NC3'}, 2.1),
    'PO4': ('Fosfato (PO4)',          4, {'atom': 'PO4'}, 2.1),
    'GL':  ('Glicerol (GL1, GL2)',    5, {'atom': ['GL1', 'GL2']}, 2.0),
    'C':   ('Caudas (C1A a C4B)',     6, {'atom': _TAILS}, 2.0),
    'W':   ('Água W (4 H₂O)',         7, {'resn': 'W'}, 1.6),
    'ION': ('Íons Na⁺ e Cl⁻',         8, {'resn': ['NA', 'CL', 'ION']}, 1.9),
}

STEPS = {
    'A': ['Estrutura', 'Mapeamento', 'Solvatação', 'Íons', 'Minimização', 'NPT', 'Produção', 'Análise'],
    'B': ['Lipídio', 'Inserção', 'Solvatação', 'Minimização', 'Automontagem', 'Equilibração', 'Análise'],
}
PART_NAME = {'A': 'Parte A, ubiquitina', 'B': 'Parte B, bicamada de DSPC'}


def style_beads(view, model=None, water=True, water_opacity=0.25, viewer=None):
    """Uma esfera por bead, cor pelo slot da serie."""
    def sel(d):
        d = dict(d)
        if model is not None:
            d['model'] = model
        return d
    kw = {'viewer': viewer} if viewer is not None else {}
    view.setStyle(sel({}), {'sphere': {'color': INK, 'radius': 1.8}}, **kw)
    for key, (_, slot, q, r) in BEADS.items():
        if key == 'W':
            continue
        view.setStyle(sel(q), {'sphere': {'color': CHART[slot - 1], 'radius': r}}, **kw)
    w = BEADS['W']
    if water:
        view.setStyle(sel(w[2]), {'sphere': {'color': CHART[w[1] - 1], 'radius': w[3], 'opacity': water_opacity}}, **kw)
    else:
        view.setStyle(sel(w[2]), {}, **kw)


def draw_box(view, gro_text, viewer=None):
    """Arestas da caixa periodica a partir da ultima linha de um .gro."""
    try:
        a, b, c = [float(x) * 10 for x in gro_text.strip().splitlines()[-1].split()[:3]]
    except (ValueError, IndexError):
        return
    P = [(0, 0, 0), (a, 0, 0), (0, b, 0), (0, 0, c), (a, b, 0), (a, 0, c), (0, b, c), (a, b, c)]
    E = [(0, 1), (0, 2), (0, 3), (1, 4), (1, 5), (2, 4), (2, 6), (3, 5), (3, 6), (4, 7), (5, 7), (6, 7)]
    kw = {'viewer': viewer} if viewer is not None else {}
    for i, j in E:
        view.addLine({'start': dict(zip('xyz', P[i])), 'end': dict(zip('xyz', P[j])),
                      'color': INK, 'linewidth': 1.5}, **kw)


def _viewer(structure, fmt=None, frames=False, water=True, box=True, height=440, style=None):
    text = Path(structure).read_text()
    fmt = fmt or structure.rsplit('.', 1)[-1]
    view = py3Dmol.view(width='100%', height=height)
    if frames:
        view.addModelsAsFrames(text, fmt)
    else:
        view.addModel(text, fmt)
    if style is not None:
        style(view)
    else:
        style_beads(view, water=water)
    if box and fmt == 'gro':
        draw_box(view, text)
    if frames:
        view.animate({'loop': 'forward', 'reps': 0, 'step': 1})
    view.setBackgroundColor('white', 0)
    view.zoomTo()
    return view._make_html()


def _legend(keys):
    items = []
    for k in keys:
        if k in BEADS:
            label, slot = BEADS[k][0], BEADS[k][1]
            items.append(f'<li><span class="swatch swatch-{slot}" aria-hidden="true"></span>{html.escape(label)}</li>')
        else:
            label, slot = k
            cls = 'swatch-ink' if slot == 'ink' else f'swatch-{slot}'
            items.append(f'<li><span class="swatch {cls}" aria-hidden="true"></span>{html.escape(label)}</li>')
    return '<ul class="stage-legend">' + ''.join(items) + '</ul>'


def _facts(facts):
    rows = ''.join(f'<div><dt>{html.escape(str(k))}</dt><dd>{html.escape(str(v))}</dd></div>'
                   for k, v in facts.items())
    return f'<dl class="stage-facts">{rows}</dl>'


def _compare(compare):
    """compare = {'aa': (forcefield, valor, texto), 'ua': ..., 'cg': ...}"""
    names = {'aa': 'All-atom', 'ua': 'United-atom', 'cg': 'Coarse-grained, esta aula'}
    cols = []
    for key in ('aa', 'ua', 'cg'):
        ff, value, text = compare[key]
        case = ' is-case' if key == 'cg' else ''
        cols.append(
            f'<div class="compare-col surface{case}">'
            f'<p class="compare-name">{names[key]}</p>'
            f'<p class="compare-ff">{html.escape(ff)}</p>'
            + (f'<p class="compare-value">{html.escape(str(value))}</p>' if value not in (None, '') else '')
            + f'<p class="compare-text">{html.escape(text)}</p></div>')
    return ('<section class="stage-compare" aria-label="Comparação entre resoluções">'
            '<p class="eyebrow stage-compare-head">Mesma etapa em três resoluções</p>'
            f'<div class="compare-grid">{"".join(cols)}</div></section>')


def _page(inner):
    css, js = _assets()
    return (f'{FONTS}<style>{css}</style>{MODE_JS}'
            f'<div class="lit-field stage-root" data-lit>{inner}</div><script>{js}</script>')


def _head(part, step, title):
    steps = STEPS[part]
    current = ' aria-current="step"'
    pills = ''.join(
        f'<li><span class="pill pill-sm glass lit lit-swell" data-lit'
        f'{current if i + 1 == step else ""}>{i + 1}. {s}</span></li>'
        for i, s in enumerate(steps))
    return (f'<header class="stage-head">'
            f'<p class="eyebrow">{PART_NAME[part]} · etapa {step} de {len(steps)}</p>'
            f'<h2 class="display stage-title">{html.escape(title)}</h2>'
            f'<ol class="stage-steps" aria-label="Etapas">{pills}</ol></header>')


def stage(part, step, title, structure=None, text=(), legend=(), facts=None, compare=None,
          frames=False, water=True, box=True, style=None, height=440, figure_html=None, wide=False):
    """Mostra um cartao de etapa: visualizacao 3D, explicacao, numeros e comparacao."""
    fig = figure_html or (_viewer(structure, frames=frames, water=water, box=box, style=style, height=height)
                          if structure else '')
    paras = ''.join(f'<p class="prose">{t}</p>' for t in text)
    body = (f'<div class="stage-body">'
            f'<figure class="stage-figure{" stage-figure-wide" if wide else ""} glass glass-frost lit lit-edge" data-lit>{fig}</figure>'
            f'<div class="stage-text{" stage-text-wide" if wide else ""}">{paras}'
            + (_legend(legend) if legend else '')
            + (_facts(facts) if facts else '')
            + '</div></div>')
    inner = (f'<article class="stage card-glass lit" data-lit>{_head(part, step, title)}{body}'
             + (_compare(compare) if compare else '') + '</article>')
    display(HTML(_page(inner)))
