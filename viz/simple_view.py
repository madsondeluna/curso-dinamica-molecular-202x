"""Visualizacao simples das etapas: estrutura em py3Dmol e explicacao em Markdown.

Sem CSS: um visualizador com fundo branco, a agua num modelo
separado e translucido (no 3Dmol a opacidade vale por modelo), a caixa
periodica em linha fina e, abaixo, legenda, numeros e explicacao formatados.
"""
from __future__ import annotations

import html
from pathlib import Path

import numpy as np
import py3Dmol
from IPython.display import HTML, Markdown, display

from pure import palette
from replica_views import SEM, show_white, white_box

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


_IONS = {'NA', 'CL', 'NA+', 'CL-'}
# cores do cartoon pela estrutura secundaria do DSSP
SS_COLORS = {'h': ('Hélice', '#9e527f'), 's': ('Folha β', '#9f8322'), 'c': ('Alça', '#9aa0a6')}


def _gro_rows(text):
    lines = text.splitlines()
    return lines[0], lines[2:2 + int(lines[1])], lines[2 + int(lines[1])]


def _xyz(row):
    return np.array([float(row[20:28]), float(row[28:36]), float(row[36:44])])


def _set_xyz(row, x):
    return f'{row[:20]}{x[0]:8.3f}{x[1]:8.3f}{x[2]:8.3f}{row[44:]}'


def _center_protein(rows, box):
    """Proteina inteira e no centro da caixa; agua e ions reembrulhados dentro dela."""
    L = np.array([float(v) for v in box.split()[:3]])
    prot = [i for i, r in enumerate(rows) if r[5:10].strip() != 'W' and r[10:15].strip() not in _IONS]
    X = np.array([_xyz(r) for r in rows])
    for k in range(1, len(prot)):            # desembrulha bead a bead ao longo da cadeia
        d = X[prot[k]] - X[prot[k - 1]]
        X[prot[k]] = X[prot[k - 1]] + d - L * np.round(d / L)
    X += L / 2 - X[prot].mean(axis=0)
    rest = np.setdiff1d(np.arange(len(rows)), prot)
    X[rest] -= L * np.floor(X[rest] / L)
    return [_set_xyz(r, x) for r, x in zip(rows, X)]


def _orient_membrane(rows, box):
    """Bicamada inteira no centro da caixa, normal na vertical da tela (eixo y), como nos filmes.
    Sem bicamada (caudas espalhadas) devolve as linhas como estao."""
    L = np.array([float(v) for v in box.split()[:3]])
    X = np.array([_xyz(r) for r in rows])
    tails = np.array([r[10:15].strip()[:1] in 'CD' and r[10:15].strip()[1:2].isdigit() for r in rows])
    ang = 2 * np.pi * X[tails] / L
    z = np.exp(1j * ang).mean(axis=0)                 # concentracao circular das caudas por eixo
    a = int(np.argmax(np.abs(z)))
    if np.abs(z[a]) < 0.5:
        return rows, box
    X[:, a] += L[a] / 2 - np.angle(z[a]) / (2 * np.pi) * L[a]
    # cada molecula inteira e com o centro dentro da caixa
    keys = [r[:10] for r in rows]
    start = 0
    for i in range(1, len(rows) + 1):
        if i == len(rows) or keys[i] != keys[start]:
            idx = np.arange(start, i)
            for k in idx[1:]:
                d = X[k] - X[start]
                X[k] = X[start] + d - L * np.round(d / L)
            X[idx] -= L * np.floor(X[idx].mean(axis=0) / L)
            start = i
    # rotacao propria levando a normal para y; a caixa acompanha
    if a == 2:
        X = np.stack([X[:, 0], X[:, 2], L[1] - X[:, 1]], axis=1); L = L[[0, 2, 1]]
    elif a == 0:
        X = np.stack([L[1] - X[:, 1], X[:, 0], X[:, 2]], axis=1); L = L[[1, 0, 2]]
    box = f'{L[0]:10.5f}{L[1]:10.5f}{L[2]:10.5f}'
    return [_set_xyz(r, x) for r, x in zip(rows, X)], box


def _split_gro(text):
    """Separa um .gro em (soluto, ions, agua), cada um um .gro valido; centra a proteina se houver."""
    title, atoms, box = _gro_rows(text)
    if any(r[10:15].strip() == 'BB' for r in atoms):
        atoms = _center_protein(atoms, box)
    elif any(r[10:15].strip() == 'PO4' for r in atoms):
        atoms, box = _orient_membrane(atoms, box)
    water = [l for l in atoms if l[5:10].strip() == 'W']
    ions = [l for l in atoms if l[5:10].strip() != 'W' and l[10:15].strip() in _IONS]
    solute = [l for l in atoms if l[5:10].strip() != 'W' and l[10:15].strip() not in _IONS]
    mk = lambda rows: ('\n'.join([title, f'{len(rows):5d}'] + rows + [box]) + '\n') if rows else None
    return mk(solute), mk(ions), mk(water)


def _box(view, gro_text):
    try:
        a, b, c = [float(x) * 10 for x in gro_text.strip().splitlines()[-1].split()[:3]]
    except (ValueError, IndexError):
        return
    P = [(0, 0, 0), (a, 0, 0), (0, b, 0), (0, 0, c), (a, b, 0), (a, 0, c), (0, b, c), (a, b, c)]
    for i, j in [(0, 1), (0, 2), (0, 3), (1, 4), (1, 5), (2, 4), (2, 6), (3, 5), (3, 6), (4, 7), (5, 7), (6, 7)]:
        view.addLine({'start': dict(zip('xyz', P[i])), 'end': dict(zip('xyz', P[j])), 'color': '#9aa0a6'})


def show(structure, frames=False, water=True, box=True, style=None, height=440, width=820):
    """Mostra uma estrutura: proteina ou bicamada centrada, agua translucida em modelo proprio, caixa."""
    text = Path(structure).read_text()
    fmt = structure.rsplit('.', 1)[-1]
    view = py3Dmol.view(width=width, height=height)
    ions_text = water_text = None
    if fmt == 'gro' and not frames:
        text, ions_text, water_text = _split_gro(text)
    if frames:
        view.addModelsAsFrames(text, fmt)
    else:
        view.addModel(text, fmt)
    if ions_text:
        view.addModel(ions_text, 'gro')
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
    view.zoomTo()
    show_white(view, Path(structure).stem)


def _legend(keys):
    items = []
    # 'SS' vira as tres entradas do cartoon (helice, folha, alca)
    keys = [e for k in keys for e in (list(SS_COLORS.values()) if k == 'SS' else [k])]
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
        display(HTML(white_box(figure_html)))
    elif structure:
        show(structure, frames=frames, water=water, box=box, style=style, height=height)
    if legend:
        _legend(legend)
    _explain(text=text, facts=facts, compare=compare)




# ---------------------------------------------------------------- compatibilidade
# Notebooks abertos antes de uma atualizacao baixam este arquivo novo e ainda importam
# estes nomes. Mantidos para que uma copia antiga do notebook nao quebre na importacao.
SS = None


def card(eyebrow, title, figure_html, text=(), legend=(), facts=None):
    """Visualizador (HTML do py3Dmol) seguido da explicacao."""
    display(Markdown(f'#### {title}'))
    display(HTML(white_box(figure_html)))
    if legend:
        _legend(legend)
    _explain(text=text, facts=facts)


def ss_cartoon(view, model, ss):
    """Cartoon colorido pela string do DSSP (H/G/I helice, E folha) no modelo indicado."""
    view.setStyle({'model': model}, {'cartoon': {'color': SS_COLORS['c'][1]}})
    helix = [i + 1 for i, s in enumerate(ss) if s in 'HGI']
    sheet = [i + 1 for i, s in enumerate(ss) if s == 'E']
    for key, res in (('h', helix), ('s', sheet)):
        if res:
            view.setStyle({'model': model, 'resi': res}, {'cartoon': {'color': SS_COLORS[key][1]}})
