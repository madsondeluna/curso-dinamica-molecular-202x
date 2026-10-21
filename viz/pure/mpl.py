"""Tema matplotlib da linguagem pure.

Uso:

    from pure import mpl
    mpl.use("light")          # ou "paper-like", "deep-blue", "dark"

    fig, ax = plt.subplots()
    ax.plot(x, y)             # ja sai com a cor do slot 1
    mpl.finish(ax, title="Cobertura por amostra", ylabel="reads")

`finish` faz o que o rcParams sozinho nao alcanca: tira a grade
vertical, recolhe as bordas, empurra o titulo para a esquerda e
deixa os rotulos no tom de tinta, nunca na cor da serie.
"""

from __future__ import annotations

from typing import Iterable

import matplotlib as _matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from cycler import cycler

from . import palette as p

_CMAPS_REGISTERED = False


def _register_cmaps() -> None:
    global _CMAPS_REGISTERED
    if _CMAPS_REGISTERED:
        return
    pairs = [
        ("pure", p.SEQUENTIAL),
        ("pure_r", list(reversed(p.SEQUENTIAL))),
        ("pure_div", p.DIVERGING),
        ("pure_div_r", list(reversed(p.DIVERGING))),
    ]
    for name, colors in pairs:
        cmap = LinearSegmentedColormap.from_list(name, colors, N=256)
        try:
            _matplotlib.colormaps.register(cmap, name=name)
        except ValueError:
            pass  # ja registrado numa sessao anterior
    _CMAPS_REGISTERED = True


def rc(mode: str = "light") -> dict[str, object]:
    """rcParams do modo pedido, sem aplicar."""
    t = p.tokens(mode)
    ink, muted = t["text"], t["muted"]
    bg, grid = t["bg"], t["border"]

    return {
        "figure.facecolor": bg,
        "figure.edgecolor": bg,
        "figure.dpi": 120,
        "figure.constrained_layout.use": True,
        "savefig.facecolor": bg,
        "savefig.edgecolor": bg,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.transparent": False,

        "axes.facecolor": bg,
        "axes.edgecolor": grid,
        "axes.linewidth": 1.0,
        "axes.labelcolor": muted,
        "axes.titlecolor": ink,
        "axes.grid": True,
        "axes.grid.axis": "y",
        "axes.axisbelow": True,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.spines.left": False,
        "axes.spines.bottom": True,
        "axes.prop_cycle": cycler(color=p.CATEGORICAL),
        "axes.titlelocation": "left",
        "axes.titlesize": 13,
        "axes.titlepad": 14,
        "axes.labelsize": 10,
        "axes.labelpad": 8,

        "grid.color": grid,
        "grid.linewidth": 0.8,
        "grid.alpha": 0.7,

        "xtick.color": muted,
        "ytick.color": muted,
        "xtick.labelcolor": muted,
        "ytick.labelcolor": muted,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.size": 4,
        "ytick.major.size": 0,
        "xtick.major.width": 1.0,
        "xtick.minor.visible": False,
        "ytick.minor.visible": False,

        "lines.linewidth": 2.0,
        "lines.markersize": 5.0,
        "lines.solid_capstyle": "round",
        "lines.antialiased": True,

        "patch.linewidth": 0.0,
        "patch.edgecolor": bg,
        "patch.force_edgecolor": False,

        "scatter.marker": "o",

        "legend.frameon": False,
        "legend.fontsize": 9,
        "legend.labelcolor": muted,
        "legend.handlelength": 1.2,
        "legend.handleheight": 1.2,
        "legend.borderpad": 0.0,
        "legend.columnspacing": 1.4,
        "legend.labelspacing": 0.5,

        "font.family": "sans-serif",
        "font.sans-serif": p.FONT_SANS,
        "font.size": 10,

        "text.color": ink,

        "image.cmap": "pure",

        "boxplot.boxprops.color": muted,
        "boxplot.whiskerprops.color": muted,
        "boxplot.capprops.color": muted,
        "boxplot.medianprops.color": ink,
        "boxplot.flierprops.markeredgecolor": muted,
    }


def use(mode: str = "light") -> None:
    """Aplica o tema globalmente e registra os colormaps."""
    _register_cmaps()
    plt.rcParams.update(rc(mode))


def context(mode: str = "light"):
    """Context manager: aplica o tema so dentro do bloco.

        with mpl.context("dark"):
            fig, ax = plt.subplots()
    """
    _register_cmaps()
    return plt.rc_context(rc(mode))


def finish(
    ax,
    title: str | None = None,
    subtitle: str | None = None,
    xlabel: str | None = None,
    ylabel: str | None = None,
    mode: str = "light",
    legend: bool | None = None,
) -> None:
    """Acabamento do eixo: hierarquia de titulo, rotulos em tinta, legenda sem caixa.

    `legend=None` deixa a legenda automatica: aparece quando ha duas ou
    mais series rotuladas, some quando ha uma so (o titulo ja a nomeia).
    """
    t = p.tokens(mode)
    ink, muted = t["text"], t["muted"]

    # heatmap e mapa de contato nao levam grade: a linha atravessaria as celulas
    if ax.get_images() or any(
        type(c).__name__ in ("QuadMesh", "PolyQuadMesh") for c in ax.collections
    ):
        ax.grid(False)

    if title:
        ax.set_title(title, color=ink, fontsize=13, loc="left", pad=18 if subtitle else 14)
    if subtitle:
        ax.text(
            0.0, 1.02, subtitle,
            transform=ax.transAxes, ha="left", va="bottom",
            color=muted, fontsize=10,
        )
    if xlabel:
        ax.set_xlabel(xlabel, color=muted)
    if ylabel:
        ax.set_ylabel(ylabel, color=muted)

    handles, labels = ax.get_legend_handles_labels()
    show = (len(labels) >= 2) if legend is None else legend
    if show and labels:
        ax.legend(handles, labels, frameon=False, loc="best", labelcolor=muted)
    elif ax.get_legend() is not None:
        ax.get_legend().remove()


def bar_gap(ax, gap: float = 0.02, mode: str = "light") -> None:
    """Insere o vao de superficie entre barras encostadas de uma pilha.

    O vao e desenhado como borda na cor do fundo, entao segmentos
    vizinhos ficam separados por uma linha do proprio papel em vez de
    encostarem cor com cor.
    """
    bg = p.tokens(mode)["bg"]
    for container in ax.containers:
        for patch in container:
            patch.set_edgecolor(bg)
            patch.set_linewidth(max(1.0, gap * 100))
