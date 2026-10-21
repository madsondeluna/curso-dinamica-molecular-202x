"""Graficos de replicas no estilo de figura de artigo.

Moldura completa com ticks para dentro nos quatro lados, titulo centrado,
uma cor por sistema e uma linha por replica (continua, tracejada,
pontilhada). A figura de media mostra a media com faixa de +- 1 DP na
mesma cor e legenda fora do eixo, a direita.
"""
from __future__ import annotations

import json
import os
import re

import matplotlib.pyplot as plt
import numpy as np

from pure import palette

INK = palette.tokens('light')['text']
BAND = palette.tokens('light')['dim']
COLORS = {'protein': '#de8552', 'membrane': '#4c72b0'}   # cores quimicas do vmd_view.sh (proteina, colina)
NAMES = {'protein': 'Ubiquitina (Martini 3)', 'membrane': 'Bicamada de DSPC (Martini 3)'}
LINESTYLES = ['-', (0, (6, 3)), (0, (1, 2))]
MARKERS = ['o', 's', '^']

STYLE = {
    'axes.spines.top': True, 'axes.spines.right': True,
    'axes.spines.left': True, 'axes.spines.bottom': True,
    'axes.edgecolor': INK, 'axes.linewidth': 0.8, 'axes.grid': False,
    'axes.titlelocation': 'center', 'axes.titlesize': 11, 'axes.titleweight': 'normal',
    'axes.labelcolor': INK, 'axes.labelsize': 10,
    'xtick.direction': 'in', 'ytick.direction': 'in',
    'xtick.top': True, 'ytick.right': True,
    'xtick.color': INK, 'ytick.color': INK,
    'xtick.labelcolor': INK, 'ytick.labelcolor': INK, 'text.color': INK,
    'figure.facecolor': 'white', 'axes.facecolor': 'white', 'savefig.facecolor': 'white',
    'xtick.major.size': 4, 'ytick.major.size': 4,
    'figure.constrained_layout.use': False,
    'legend.frameon': False, 'legend.handlelength': 3.2, 'legend.fontsize': 9,
}


def load_summaries(base_dir, repo=None, replicas=(1, 2, 3)):
    """summary.json de cada replica: local se existir, senao baixado do repositorio."""
    sums = {}
    cache = os.path.join(base_dir, 'summaries')
    os.makedirs(cache, exist_ok=True)
    for r in replicas:
        own = os.path.join(base_dir, f'rep{r}', 'summary.json')
        remote = os.path.join(cache, f'rep{r}.json')
        if os.path.exists(own):
            sums[r] = json.load(open(own))
            continue
        if repo and (not os.path.exists(remote) or os.path.getsize(remote) == 0):
            os.system(f'wget -q {repo}/results/summary_rep{r}.json -O {remote}')
        if os.path.exists(remote) and os.path.getsize(remote):
            sums[r] = json.load(open(remote))
    return sums


def helix_bands(ss, min_len=4):
    """Intervalos de residuos (base 1) em alfa-helice segundo a string DSSP."""
    return [(m.start() + 1, m.end()) for m in re.finditer('H{%d,}' % min_len, ss or '')]


def _runs(sums, key, xk, yk):
    runs = [(r, np.array(s[key][xk]), np.array(s[key][yk])) for r, s in sorted(sums.items()) if key in s]
    if not runs:
        return [], None, None
    n = min(len(y) for _, _, y in runs)
    return runs, runs[0][1][:n], np.array([y[:n] for _, _, y in runs])


def _energy_ylim(ax, runs):
    """Ignora o primeiro quadro (energia da estrutura minimizada) no limite do eixo."""
    vals = np.concatenate([y[x > 0.2] for _, x, y in runs])
    pad = (vals.max() - vals.min()) * 0.15
    ax.set_ylim(vals.min() - pad, vals.max() + pad)


def white(fig):
    """Fundo branco na tela e no arquivo, independente do tema global."""
    fig.patch.set_facecolor('white')
    for ax in fig.axes:
        ax.set_facecolor('white')
    return fig


def _bands(ax, bands):
    for i, (a, b) in enumerate(bands):
        ax.axvspan(a - 0.5, b + 0.5, color=BAND, lw=0, zorder=0,
                   label='α-hélice (DSSP)' if i == 0 else None)


PANELS = [
    # sistema, x, y, titulo, rotulo x, rotulo y, escala x, posicao da legenda
    ('protein', 'rmsf_res', 'rmsf', 'RMSF por resíduo', 'Resíduo (1UBQ)', 'BB RMSF (nm)', 1.0, 'upper left'),
    ('protein', 'rmsd_t', 'rmsd', 'RMSD', 'Tempo (ns)', 'BB RMSD (nm)', 1e-3, 'lower right'),
    ('membrane', 't_ns', 'S', 'Parâmetro de ordem', 'Tempo (ns)', 'S', 1.0, 'upper left'),
    ('membrane', 'pot_t', 'pot', 'Energia potencial', 'Tempo (ns)', 'Energia potencial (kJ/mol)', 1.0, 'upper right'),
]


def plot_replicas(sums, ss=None):
    """Uma linha por replica, mesma cor por sistema."""
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(2, 2, figsize=(12, 8.5))
        for ax, (key, xk, yk, title, xl, yl, xs, loc) in zip(axes.ravel(), PANELS):
            runs, _, _ = _runs(sums, key, xk, yk)
            if key == 'protein' and xk == 'rmsf_res':
                _bands(ax, helix_bands(ss))
            for (r, x, y), ls in zip(runs, LINESTYLES):
                ax.plot(x * xs, y, color=COLORS[key], ls=ls, lw=1.1, label=f'Réplica {r}')
            if yk == 'S':
                ax.set_ylim(0, 1)
            if yk == 'pot' and runs:
                _energy_ylim(ax, runs)
            ax.set_title(f'{NAMES[key]}: {title}')
            ax.set_xlabel(xl); ax.set_ylabel(yl)
            ax.legend(loc=loc)
        fig.tight_layout()
    return white(fig)


def plot_means(sums, ss=None):
    """Media +- DP sobre as replicas, legenda fora do eixo."""
    rows = [PANELS[1], PANELS[0], PANELS[2], PANELS[3]]
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(len(rows), 1, figsize=(11, 13))
        for ax, (key, xk, yk, title, xl, yl, xs, _) in zip(axes, rows):
            runs, x, Y = _runs(sums, key, xk, yk)
            if not runs:
                continue
            if key == 'protein' and xk == 'rmsf_res':
                _bands(ax, helix_bands(ss))
            m = Y.mean(0)
            sd = Y.std(0, ddof=1) if len(Y) > 1 else np.zeros_like(m)
            ax.fill_between(x * xs, m - sd, m + sd, color=COLORS[key], alpha=0.18, lw=0)
            ax.plot(x * xs, m, color=COLORS[key], lw=1.2, label=NAMES[key])
            if yk == 'S':
                ax.set_ylim(0, 1)
            if yk == 'pot':
                _energy_ylim(ax, runs)
            ax.set_title(title)
            ax.set_xlabel(xl); ax.set_ylabel(yl)
            ax.legend(loc='upper left', bbox_to_anchor=(1.02, 1.0), borderaxespad=0)
        fig.suptitle(f'Média ± DP sobre {len(sums)} réplicas', y=0.995, fontsize=11)
        fig.tight_layout()
    return white(fig)


def plot_scalars(sums):
    """Valores finais da membrana: um marcador por replica, media e DP em tinta."""
    metrics = [('t_form', 'Formação da bicamada', 'Tempo (ns)'),
               ('d_pp', 'Espessura D_PP', 'D_PP (nm)'),
               ('apl', 'Área por lipídio', 'Área (nm²)')]
    stats = {}
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(1, 3, figsize=(12, 3.4))
        for ax, (k, title, xl) in zip(axes, metrics):
            vals = [(r, s['membrane'][k]) for r, s in sorted(sums.items())
                    if 'membrane' in s and s['membrane'].get(k) is not None]
            for (r, v), mk in zip(vals, MARKERS):
                ax.plot([v], [1], ls='none', marker=mk, ms=8, mfc='none', mec=COLORS['membrane'],
                        mew=1.4, label=f'Réplica {r}')
            if vals:
                v = np.array([x for _, x in vals])
                sd = v.std(ddof=1) if len(v) > 1 else 0.0
                ax.errorbar([v.mean()], [0], xerr=[sd], fmt='o', color=INK, capsize=4, ms=6, label='Média ± DP')
                stats[k] = (v.mean(), sd, len(v))
                pad = max(sd, abs(v.mean()) * 0.02, 1e-3) * 3
                ax.set_xlim(max(v.min() - pad, 0) if k == 't_form' else v.min() - pad, v.max() + pad)
            missing = [r for r, s in sorted(sums.items()) if 'membrane' in s and s['membrane'].get(k) is None]
            if missing:
                ax.text(0.02, 0.04, 'sem bicamada: réplica ' + ', '.join(map(str, missing)),
                        transform=ax.transAxes, fontsize=8, color=INK, va='bottom')
            ax.set_ylim(-0.8, 1.8); ax.set_yticks([])
            ax.set_title(title); ax.set_xlabel(xl)
            ax.legend(loc='upper right', fontsize=8, handlelength=1.5)
        fig.tight_layout()
    return white(fig), stats


def use_paper_style(dpi=130):
    """Aplica o estilo de figura de artigo ao matplotlib inteiro (graficos do notebook)."""
    plt.rcParams.update(STYLE)
    plt.rcParams.update({'figure.dpi': dpi, 'savefig.dpi': 200, 'font.size': 10,
                         'lines.linewidth': 1.2, 'legend.fontsize': 9})
    try:   # telas de alta densidade no Jupyter e no Colab
        from matplotlib_inline.backend_inline import set_matplotlib_formats
        set_matplotlib_formats('retina')
    except Exception:
        pass


def finish(ax, title=None, subtitle=None, xlabel=None, ylabel=None, legend=None):
    """Acabamento no estilo de artigo: titulo centrado, unidade no eixo y, legenda sem caixa."""
    if title:
        ax.set_title(title)
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel or subtitle:
        ax.set_ylabel(f'{ylabel} ({subtitle})' if ylabel and subtitle else (ylabel or subtitle))
    handles, labels = ax.get_legend_handles_labels()
    if (legend is None and len(labels) >= 2) or legend:
        ax.legend(frameon=False)


# ---------------------------------------------------------------- um grafico por metrica
METRICS = {
    'rmsd': PANELS[1], 'rmsf': PANELS[0], 'S': PANELS[2], 'pot': PANELS[3],
}
SCALARS = {
    't_form': ('Formação da bicamada', 'Tempo (ns)'),
    'd_pp': ('Espessura D_PP', 'D_PP (nm)'),
    'apl': ('Área por lipídio', 'Área (nm²)'),
}


def plot_metric(sums, metric, mode='replicas', ss=None, figsize=(8, 4)):
    """Um painel: 'replicas' (uma linha por replica) ou 'mean' (media com faixa de +- 1 DP)."""
    key, xk, yk, title, xl, yl, xs, loc = METRICS[metric]
    with plt.rc_context(STYLE):
        fig, ax = plt.subplots(figsize=figsize)
        runs, x, Y = _runs(sums, key, xk, yk)
        if key == 'protein' and xk == 'rmsf_res':
            _bands(ax, helix_bands(ss))
        if mode == 'replicas':
            for (r, rx, ry), ls in zip(runs, LINESTYLES):
                ax.plot(rx * xs, ry, color=COLORS[key], ls=ls, lw=1.1, label=f'Réplica {r}')
            ax.legend(loc=loc)
            ax.set_title(f'{NAMES[key]}: {title}')
        else:
            m = Y.mean(0)
            sd = Y.std(0, ddof=1) if len(Y) > 1 else np.zeros_like(m)
            ax.fill_between(x * xs, m - sd, m + sd, color=COLORS[key], alpha=0.18, lw=0, label='± 1 DP')
            ax.plot(x * xs, m, color=COLORS[key], lw=1.3, label=f'Média de {len(Y)} réplicas')
            ax.legend(loc=loc)
            ax.set_title(f'{NAMES[key]}: {title} (média ± DP)')
        if yk == 'S':
            ax.set_ylim(0, 1)
        if yk == 'pot' and runs:
            _energy_ylim(ax, runs)
        ax.set_xlabel(xl); ax.set_ylabel(yl)
        fig.tight_layout()
    return white(fig)


def plot_scalar(sums, k, figsize=(7, 2.8)):
    """Um valor final da membrana: um marcador por replica, media e DP em tinta."""
    title, xl = SCALARS[k]
    with plt.rc_context(STYLE):
        fig, ax = plt.subplots(figsize=figsize)
        vals = [(r, s['membrane'][k]) for r, s in sorted(sums.items())
                if 'membrane' in s and s['membrane'].get(k) is not None]
        for (r, v), mk in zip(vals, MARKERS):
            ax.plot([v], [1], ls='none', marker=mk, ms=9, mfc='none', mec=COLORS['membrane'], mew=1.4,
                    label=f'Réplica {r}')
        stat = None
        if vals:
            v = np.array([x for _, x in vals])
            sd = v.std(ddof=1) if len(v) > 1 else 0.0
            ax.errorbar([v.mean()], [0], xerr=[sd], fmt='o', color=INK, capsize=4, ms=6, label='Média ± DP')
            stat = (v.mean(), sd, len(v))
            pad = max(sd, abs(v.mean()) * 0.02, 1e-3) * 3
            ax.set_xlim(max(v.min() - pad, 0) if k == 't_form' else v.min() - pad, v.max() + pad)
        missing = [r for r, s in sorted(sums.items()) if 'membrane' in s and s['membrane'].get(k) is None]
        if missing:
            ax.text(0.01, 0.05, 'sem bicamada: réplica ' + ', '.join(map(str, missing)),
                    transform=ax.transAxes, fontsize=8, color=INK)
        ax.set_ylim(-0.8, 1.8); ax.set_yticks([])
        ax.set_title(f'{NAMES["membrane"]}: {title}'); ax.set_xlabel(xl)
        ax.legend(loc='center left', bbox_to_anchor=(1.01, 0.5), fontsize=8, handlelength=1.5)
        fig.tight_layout()
    return white(fig), stat

