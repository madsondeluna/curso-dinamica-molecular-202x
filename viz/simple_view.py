"""Visualizacao simples das etapas: estrutura em py3Dmol e explicacao em Markdown.

Sem cartoes nem CSS: um visualizador com fundo branco, a agua num modelo
separado e translucido (no 3Dmol a opacidade vale por modelo), a caixa
periodica em linha fina e, abaixo, legenda, numeros e explicacao formatados.
"""
from __future__ import annotations

import html
from pathlib import Path

import py3Dmol
from IPython.display import HTML, Markdown, display

from pure import palette
from replica_views import SEM

CHART = palette.CATEGORICAL
INK = '#4b5563'

_TAILS = [f'C{i}{s}' for i in range(1, 6) for s in 'AB'] + [f'D{i}{s}' for i in range(1, 6) for s in 'AB']
# nome na legenda, cor quimica, selecao py3Dmol, raio (A)
BEADS = {
    'BB':  ('Backbone (BB)',       SEM['protein'],    {'atom': 'BB'}, 1.9),
    'SC':  ('Cadeia lateral (SC)', SEM['protein_sc'], {'atom': ['SC1', 'SC2', 'SC3', 'SC4', 'SC5']}, 1.6),
    'NC3': ('Colina (NC3)',        SEM['head_cat'],   {'atom': 'NC3'}, 2.0),
    'PO4': ('Fosfato (PO4)',       SEM['po4'],        {'atom': 'PO4'}, 2.0),
    'GL':  ('Glicerol (GL1, GL2)', SEM['glyc'],       {'atom': ['GL1', 'GL2']}, 1.8),
    'C':   ('Caudas (C1A a C4B)',  SEM['tail'],       {'atom': _TAILS}, 1.6),
    'W':   ('Água W (4 H₂O)',      SEM['water'],      {'resn': 'W'}, 1.6),
    'NA':  ('Na⁺',                 SEM['na'],         {'atom': 'NA'}, 1.6),
    'CL':  ('Cl⁻',                 SEM['cl'],         {'atom': 'CL'}, 1.6),
}


def style_beads(view, model=None, water=True, water_opacity=0.15, viewer=None):
    """Uma esfera por bead, cor quimica do bead (agua tratada a parte)."""
    def sel(d):
        d = dict(d)
        if model is not None:
            d['model'] = model
        return d
    kw = {'viewer': viewer} if viewer is not None else {}
    view.setStyle(sel({}), {'sphere': {'color': INK, 'radius': 1.6}}, **kw)
    for key, (_, color, q, r) in BEADS.items():
        if key != 'W':
            view.setStyle(sel(q), {'sphere': {'color': color, 'radius': r}}, **kw)


def _split_gro(text):
    """Separa um .gro em (sistema sem agua, so agua), cada um um .gro valido."""
    lines = text.splitlines()
    title, atoms, box = lines[0], lines[2:-1], lines[-1]
    solute = [l for l in atoms if l[5:10].strip() != 'W']
    water = [l for l in atoms if l[5:10].strip() == 'W']
    mk = lambda rows: '\n'.join([title, f'{len(rows):5d}'] + rows + [box]) + '\n'
    return mk(solute), (mk(water) if water else None)


def _box(view, gro_text):
    try:
        a, b, c = [float(x) * 10 for x in gro_text.strip().splitlines()[-1].split()[:3]]
    except (ValueError, IndexError):
        return
    P = [(0, 0, 0), (a, 0, 0), (0, b, 0), (0, 0, c), (a, b, 0), (a, 0, c), (0, b, c), (a, b, c)]
    for i, j in [(0, 1), (0, 2), (0, 3), (1, 4), (1, 5), (2, 4), (2, 6), (3, 5), (3, 6), (4, 7), (5, 7), (6, 7)]:
        view.addLine({'start': dict(zip('xyz', P[i])), 'end': dict(zip('xyz', P[j])), 'color': '#9aa0a6'})


def show(structure, frames=False, water=True, box=True, style=None, height=440, width=820):
    """Mostra uma estrutura: agua translucida em modelo proprio, caixa em linha fina."""
    text = Path(structure).read_text()
    fmt = structure.rsplit('.', 1)[-1]
    view = py3Dmol.view(width=width, height=height)
    water_text = None
    if fmt == 'gro' and not frames:
        text, water_text = _split_gro(text)
    if frames:
        view.addModelsAsFrames(text, fmt)
    else:
        view.addModel(text, fmt)
    if style is not None:
        style(view)
    else:
        style_beads(view)
    if water and water_text:
        view.addModel(water_text, 'gro')
        # agua como lamina continua (efeito QuickSurf): esferas grandes, sobrepostas e translucidas
        view.setStyle({'resn': 'W'}, {'sphere': {'color': '#5b9bd5', 'radius': 3.0, 'opacity': 0.3}})
    if box and fmt == 'gro':
        _box(view, text)
    if frames:
        view.animate({'loop': 'forward', 'reps': 0, 'interval': 150})
    view.setBackgroundColor('white')
    view.zoomTo()
    view.show()


def _legend(keys):
    items = []
    for k in keys:
        if isinstance(k, str) and k in BEADS:
            label, color = BEADS[k][0], BEADS[k][1]
        else:
            label, color = k[0], (k[1] if isinstance(k[1], str) and k[1].startswith('#') else INK)
        items.append(f'<span style="color:{color};font-size:1.2em">&#9679;</span>&nbsp;{html.escape(label)}')
    display(HTML('<p>' + ' &nbsp;&nbsp; '.join(items) + '</p>'))


def _explain(title=None, text=(), facts=None, compare=None):
    md = []
    if facts:
        md.append('| Medida | Valor |\n|---|---|\n' + '\n'.join(f'| {k} | {v} |' for k, v in facts.items()))
    md += list(text)
    if compare:
        names = {'aa': 'All-atom', 'ua': 'United-atom', 'cg': '**Coarse-grained (esta aula)**'}
        rows = ['| | ' + ' | '.join(names[k] for k in ('aa', 'ua', 'cg')) + ' |', '|---|---|---|---|',
                '| Modelo | ' + ' | '.join(compare[k][0] for k in ('aa', 'ua', 'cg')) + ' |']
        if any(compare[k][1] for k in compare):
            rows.append('| Valor | ' + ' | '.join(str(compare[k][1]) for k in ('aa', 'ua', 'cg')) + ' |')
        rows.append('| O que muda | ' + ' | '.join(compare[k][2] for k in ('aa', 'ua', 'cg')) + ' |')
        md.append('**Mesma etapa em três resoluções**\n\n' + '\n'.join(rows))
    display(Markdown('\n\n'.join(md)))


def stage(part, step, title, structure=None, text=(), legend=(), facts=None, compare=None,
          frames=False, water=True, box=True, style=None, height=440, figure_html=None, wide=False):
    """Estrutura da etapa seguida da explicacao formatada."""
    display(Markdown(f'#### {title}'))
    if callable(figure_html):
        figure_html()                      # funcao que desenha o visualizador
    elif figure_html:
        display(HTML(figure_html))
    elif structure:
        show(structure, frames=frames, water=water, box=box, style=style, height=height)
    if legend:
        _legend(legend)
    _explain(text=text, facts=facts, compare=compare)


def card(eyebrow, title, figure_html, text=(), legend=(), facts=None):
    """Visualizadores das replicas (HTML do py3Dmol) seguidos da explicacao."""
    display(Markdown(f'#### {title}'))
    display(HTML(figure_html))
    if legend:
        _legend(legend)
    _explain(text=text, facts=facts)
