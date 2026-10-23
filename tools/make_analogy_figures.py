"""Gera as ilustracoes de analogias e exemplos de cada secao da aula.

Esquemas ilustrativos, no mesmo estilo de tools/make_figures.py; a unica
figura calculada e a das replicas (Langevin em poco duplo, modelo de brinquedo).

Uso: uv run --no-project --with matplotlib --with numpy --with mdtraj python tools/make_analogy_figures.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Ellipse, FancyArrowPatch, FancyBboxPatch, Polygon, Rectangle

sys.path.insert(0, str(Path(__file__).resolve().parent))
from make_figures import C, INK, BAND, OUT, save, frame_off, _lip, LIPID  # noqa: E402
from replica_plots import LINESTYLES  # noqa: E402

WATER = '#a6cfe3'
GREY = '#9a9a9a'


def arrow(ax, p0, p1, color=INK, lw=1.4, ms=14, z=4, style='-|>'):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle=style, mutation_scale=ms, color=color, lw=lw, zorder=z))


# ---------------------------------------------------------------- 1. fotografia e filme
def _blob(ax, cx, cy, width, depth=0.55, scale=1.0, ligand=False):
    """Proteina como contorno fechado com um sulco em theta = 0 de largura `width`."""
    th = np.linspace(-np.pi, np.pi, 400)
    r = 1 + 0.07 * np.sin(5 * th + 0.6) + 0.04 * np.cos(9 * th) - depth * np.exp(-(th / width) ** 2)
    x, y = cx + scale * r * np.cos(th), cy + scale * r * np.sin(th)
    ax.fill(x, y, color=C[0], alpha=0.22, lw=0, zorder=1)
    ax.plot(x, y, color=C[0], lw=1.6, zorder=2)
    if ligand:
        ax.add_patch(Circle((cx + scale * 0.72, cy), 0.17 * scale, color=C[5], ec=INK, lw=0.6, zorder=3))


def fig_photo_film():
    fig, ax = plt.subplots(figsize=(13, 3.6))
    frame_off(ax); ax.set_aspect('equal')
    ax.set_xlim(-1.5, 15.6); ax.set_ylim(-2.0, 1.75)
    _blob(ax, 0, 0, 0.05)
    ax.text(0, -1.55, 'Fotografia: estrutura do PDB\n(sulco fechado)', ha='center', va='top', fontsize=9)
    ax.plot([1.75, 1.75], [-1.8, 1.5], color=BAND, lw=1)
    widths = [0.05, 0.12, 0.3, 0.42, 0.1]
    for k, w in enumerate(widths):
        cx = 3.4 + 2.75 * k
        _blob(ax, cx, 0, w, ligand=(k == 3))
        ax.text(cx, 1.3, f'quadro {k + 1}', ha='center', va='bottom', fontsize=8.5)
    ax.text(3.4 + 2.75 * 2, -1.55, 'Filme: dinâmica molecular. O sulco abre e fecha; no quadro 4 um ligante cabe nele',
            ha='center', va='top', fontsize=9)
    arrow(ax, (2.6, -1.25), (14.8, -1.25), color=GREY, lw=1.0, ms=10)
    fig.suptitle('A estrutura cristalográfica é uma fotografia; a MD é o filme (esquema ilustrativo)', fontsize=11, y=1.0)
    save(fig, 'fund_foto_filme.png')


# ---------------------------------------------------------------- 2. anemia falciforme
def _hb(ax, x, y, knob=None, notch=True, w=1.5, h=1.1):
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle='round,pad=0,rounding_size=0.45',
                                fc=C[7], alpha=0.25, ec=C[7], lw=1.2, zorder=1))
    if notch:   # bolsao hidrofobico do lado esquerdo
        ax.add_patch(Ellipse((x - w / 2 + 0.12, y - 0.1), 0.26, 0.36, fc='#e6d9a8', ec=INK, lw=0.5, zorder=2))
    if knob == 'glu':
        ax.add_patch(Circle((x + w / 2 + 0.06, y + 0.15), 0.16, fc='#c44e52', ec=INK, lw=0.5, zorder=3))
        ax.text(x + w / 2 + 0.06, y + 0.15, '−', ha='center', va='center', color='white', fontsize=9, fontweight='bold', zorder=4)
    elif knob == 'val':
        ax.add_patch(Circle((x + w / 2 + 0.06, y - 0.1), 0.16, fc='#e6d9a8', ec=INK, lw=0.5, zorder=3))


def _sickle(ax, x, y, r=0.95):
    t = np.linspace(np.pi + 0.3, 2 * np.pi - 0.3, 120)
    outer = np.c_[x + r * np.cos(t), y + 0.7 * r * np.sin(t)]
    xl, xr = outer[0, 0], outer[-1, 0]
    ye = outer[0, 1]
    xs = np.linspace(xr, xl, 120)
    inner = np.c_[xs, ye - 0.27 * r * (1 - ((xs - x) / (xr - x)) ** 2)]
    ax.add_patch(Polygon(np.r_[outer, inner], closed=True, fc='#c44e52', alpha=0.75, ec='#8f2f33', lw=1, zorder=2))


def fig_sickle():
    fig, axs = plt.subplots(1, 3, figsize=(14, 3.8), gridspec_kw={'width_ratios': [1, 1.6, 1]})
    for ax in axs:
        frame_off(ax); ax.set_aspect('equal')
    a, b, c = axs
    a.set_xlim(-0.5, 4.6); a.set_ylim(-1.6, 1.6)
    _hb(a, 1.0, 0.4, knob='glu'); _hb(a, 3.2, -0.5, knob='glu')
    a.set_title('HbA: Glu6 com carga, na água', fontsize=10)
    a.text(2.05, -1.45, 'moléculas livres em solução', ha='center', fontsize=8.5)
    b.set_xlim(-0.6, 7.4); b.set_ylim(-1.6, 1.6)
    for k in range(4):
        _hb(b, 0.5 + 1.82 * k, 0.0, knob='val' if k < 3 else None)
    b.set_title('HbS desoxigenada: Val6 se encaixa no bolsão hidrofóbico da vizinha', fontsize=10)
    b.text(3.2, -1.45, 'as moléculas se encadeiam em fibras longas', ha='center', fontsize=8.5)
    c.set_xlim(-1.4, 4.0); c.set_ylim(-1.6, 1.6)
    c.add_patch(Ellipse((0, 0.1), 1.9, 1.25, fc='#c44e52', alpha=0.75, ec='#8f2f33', lw=1))
    c.add_patch(Ellipse((0, 0.1), 0.8, 0.45, fc='#e8a3a5', lw=0))
    _sickle(c, 2.5, 0.35)
    c.text(0, -1.05, 'normal', ha='center', fontsize=8.5); c.text(2.5, -1.05, 'em foice', ha='center', fontsize=8.5)
    c.set_title('Hemácias', fontsize=10)
    from matplotlib.lines import Line2D
    hs = [Line2D([], [], ls='none', marker='o', ms=9, mfc='#c44e52', mec=INK, mew=0.5),
          Line2D([], [], ls='none', marker='o', ms=9, mfc='#e6d9a8', mec=INK, mew=0.5)]
    fig.legend(hs, ['Glu (ácido glutâmico, carga −)', 'Val e bolsão hidrofóbico'], loc='lower center', ncol=2,
               frameon=False, fontsize=8.5, bbox_to_anchor=(0.5, -0.04))
    fig.suptitle('Anemia falciforme: uma troca de resíduo (Glu6Val na cadeia β) muda o comportamento coletivo (esquema)',
                 fontsize=11, y=1.02)
    save(fig, 'fund_falciforme.png')


# ---------------------------------------------------------------- 3. paisagem e doencas de enovelamento
def _well(x, x0, depth, width):
    return -depth * np.exp(-((x - x0) / width) ** 2)


def fig_misfolding():
    x = np.linspace(0, 10, 800)
    base = 0.35 * np.exp(-((x - 0.6) / 2.2) ** 2) * 2.4 + 0.05 * np.sin(6 * x) * np.exp(-((x - 1.2) / 1.5) ** 2)
    wild = base + _well(x, 3.2, 3.0, 0.55) + _well(x, 5.6, 1.2, 0.4) + _well(x, 8.4, 3.6, 0.5) + 0.6
    mut = base + _well(x, 3.2, 1.6, 0.55) + _well(x, 5.6, 1.6, 0.4) + _well(x, 8.4, 3.6, 0.5) + 0.6
    fig, ax = plt.subplots(figsize=(11, 4.4))
    ax.plot(x, wild, color=C[0], lw=2, label='proteína selvagem')
    ax.plot(x, mut, color=C[7], lw=2, ls=(0, (6, 3)), label='mutante (ex.: CFTR ΔF508): vale nativo mais raso')
    i = np.argmin(np.abs(x - 3.2))
    ax.plot(3.2, wild[i] + 0.2, 'o', ms=14, color=C[0], mec=INK, mew=0.6, zorder=5)
    for xx, lab in [(0.9, 'desenovelada'), (3.2, 'nativa\n(funcional)'), (5.6, 'mal enovelada\n(armadilha)'),
                    (8.4, 'agregado amiloide\n(Alzheimer, príons)')]:
        ax.text(xx, 1.55, lab, ha='center', va='bottom', fontsize=9)
    ax.set_ylim(-3.4, 2.4); ax.set_xlim(0, 10)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel('Conformação (coordenada ilustrativa)'); ax.set_ylabel('Energia livre')
    ax.legend(frameon=False, loc='lower left', fontsize=9)
    ax.set_title('A proteína como uma bolinha no relevo: uma mutação muda a profundidade dos vales (esquema)')
    save(fig, 'fund_paisagem_doenca.png')


# ---------------------------------------------------------------- 4. marcos de tempo simulado
def fig_milestones():
    fig, ax = plt.subplots(figsize=(12, 4.2))
    ax.set_xscale('log'); ax.set_xlim(1e-15, 10); ax.set_ylim(-2.9, 3.0)
    ax.tick_params(axis='x', which='both', top=False)
    ax.tick_params(axis='x', which='minor', bottom=False)
    ax.set_yticks([])
    for sp in ('left', 'right', 'top'):
        ax.spines[sp].set_visible(False)
    ax.set_xticks([1e-15, 1e-12, 1e-9, 1e-6, 1e-3, 1])
    ax.set_xticklabels(['1 fs', '1 ps', '1 ns', '1 µs', '1 ms', '1 s'])
    ax.set_xlabel('Tempo simulado (escala logarítmica)')
    events = [(5e-15, 1e-13, 'vibração de ligação'), (1e-12, 1e-9, 'rotação de cadeias laterais'),
              (1e-9, 1e-6, 'movimento de alças'), (1e-6, 1e-3, 'enovelamento de proteínas pequenas'),
              (1e-3, 1, 'abertura de canais')]
    for k, (a, b, lab) in enumerate(events):
        y = -0.8 - 1.0 * (k % 2)
        ax.plot([a, b], [y, y], color=BAND, lw=7, solid_capstyle='butt')
        ax.text(np.sqrt(a * b), y - 0.3, lab, ha='center', va='top', fontsize=8.2)
    marks = [(2e-15, 'passo de tempo\nall-atom: 2 fs', C[2], 0.5),
             (9e-12, '1977: primeira MD de\nproteína (BPTI), < 10 ps', C[1], 1.5),
             (15e-9, 'esta aula: 15 ns de Martini\npor réplica, num notebook', C[3], 0.5),
             (1e-3, '2010: Anton, BPTI,\n1 ms all-atom', C[0], 1.5),
             (0.1, '2020: Folding@home,\nSARS-CoV-2, 0,1 s somados', C[7], 0.5)]
    for t, lab, col, h in marks:
        ax.plot([t, t], [0, h], color=col, lw=1.3)
        ax.plot(t, 0, 'o', color=col, ms=7, zorder=5)
        ax.text(t, h + 0.08, lab, ha='center', va='bottom', fontsize=8.5, color=INK)
    ax.axhline(0, color=INK, lw=0.8)
    ax.spines['bottom'].set_position(('data', -2.9))
    ax.set_title('Quanto tempo a MD alcança, e o que acontece em cada escala')
    save(fig, 'fund_marcos_md.png')


# ---------------------------------------------------------------- 5. ubiquitina como etiqueta
def _ub_chain(ax, x, y, n=4, r=0.22):
    for k in range(n):
        ax.add_patch(Circle((x + k * 0.36, y + 0.12 * (k % 2)), r, fc=C[3], ec=INK, lw=0.6, zorder=4))


def _substrate(ax, x, y, w=1.2, h=0.8, alpha=0.35):
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle='round,pad=0,rounding_size=0.3',
                                fc=C[4], alpha=alpha, ec=C[4], lw=1.2, zorder=2))


def _proteasome(ax, x, y):
    for k, col in enumerate([C[6], C[1], C[1], C[6]]):
        ax.add_patch(FancyBboxPatch((x - 0.55, y - 1.0 + k * 0.5), 1.1, 0.46, boxstyle='round,pad=0,rounding_size=0.12',
                                    fc=col, alpha=0.45, ec=INK, lw=0.5, zorder=2))


def fig_ubiquitin_tag():
    fig, ax = plt.subplots(figsize=(13, 5.2))
    frame_off(ax); ax.set_aspect('equal')
    ax.set_xlim(-1, 15.5); ax.set_ylim(-4.6, 2.4)
    y = 0.6
    _substrate(ax, 0.5, y); ax.text(0.5, y - 0.75, 'proteína velha\nou danificada', ha='center', va='top', fontsize=8.5)
    ax.add_patch(Circle((3.3, y), 0.55, fc=C[5], alpha=0.5, ec=INK, lw=0.6))
    ax.text(3.3, y, 'E3', ha='center', va='center', fontsize=9, fontweight='bold')
    ax.text(3.3, y - 0.75, 'ligase E3\n(ex.: parkina)', ha='center', va='top', fontsize=8.5)
    arrow(ax, (1.25, y), (2.65, y))
    arrow(ax, (3.95, y), (5.35, y))
    _substrate(ax, 6.2, y); _ub_chain(ax, 6.55, y + 0.62)
    ax.text(6.2, y - 0.75, 'etiqueta: cadeia de\n≥ 4 ubiquitinas (Lys48)', ha='center', va='top', fontsize=8.5)
    arrow(ax, (7.4, y), (8.8, y))
    _proteasome(ax, 9.8, y)
    ax.text(9.8, y - 1.25, 'proteassomo', ha='center', va='top', fontsize=8.5)
    arrow(ax, (10.5, y), (11.9, y))
    for k in range(5):
        ax.add_patch(FancyBboxPatch((12.2 + 0.5 * k, y - 0.4 + 0.18 * (k % 3)), 0.32, 0.18,
                                    boxstyle='round,pad=0,rounding_size=0.08', fc=C[4], alpha=0.5, lw=0))
    ax.text(13.2, y - 0.75, 'peptídeos\nreciclados', ha='center', va='top', fontsize=8.5)
    ax.text(-0.9, 1.9, 'Célula saudável', fontsize=10, fontweight='bold')
    y2 = -2.9
    ax.text(-0.9, -1.75, 'Parkina mutada (parkinsonismo juvenil)', fontsize=10, fontweight='bold')
    _substrate(ax, 0.5, y2)
    ax.add_patch(Circle((3.3, y2), 0.55, fc=C[5], alpha=0.25, ec=INK, lw=0.6, ls='--'))
    ax.text(3.3, y2, 'E3', ha='center', va='center', fontsize=9, color=GREY)
    ax.plot([2.9, 3.7], [y2 - 0.4, y2 + 0.4], color='#c44e52', lw=2)
    arrow(ax, (1.25, y2), (2.65, y2))
    arrow(ax, (3.95, y2), (5.35, y2), color=GREY)
    for k, (dx, dy) in enumerate([(0, 0), (0.9, 0.35), (1.7, -0.2), (0.6, -0.65), (1.5, 0.75)]):
        _substrate(ax, 6.2 + dx, y2 + dy, w=0.8, h=0.5, alpha=0.45)
    ax.text(7.0, y2 - 1.0, 'sem etiqueta, a proteína não é degradada e se acumula', ha='center', va='top', fontsize=8.5)
    fig.suptitle('A ubiquitina como etiqueta de endereço para o proteassomo (esquema)', fontsize=11, y=0.98)
    save(fig, 'fund_ubiquitina_etiqueta.png')


# ---------------------------------------------------------------- 6. oleo, bicamada e surfactante
def _alveolus(ax, cx, cy, R, n_lip):
    """Alveolo em corte: parede, camada de agua e surfactante com as caudas voltadas para o ar."""
    ax.add_patch(Circle((cx, cy), R + 0.28, fc='#e8b4b0', ec='#b5706b', lw=1, zorder=0))
    ax.add_patch(Circle((cx, cy), R, fc=WATER, ec='none', alpha=0.9, zorder=1))
    ax.add_patch(Circle((cx, cy), R - 0.22, fc='white', ec='none', zorder=1))
    ax.text(cx, cy, 'ar', ha='center', va='center', fontsize=9, color=INK)
    for k in range(n_lip):
        t = 2 * np.pi * k / n_lip + 0.3
        u = np.array([np.cos(t), np.sin(t)])
        head = np.array([cx, cy]) + (R - 0.2) * u
        for off in (-0.06, 0.06):
            perp = np.array([-u[1], u[0]]) * off
            tail = head - 0.42 * u + perp
            ax.plot([head[0] + perp[0], tail[0]], [head[1] + perp[1], tail[1]], color=GREY, lw=1.0, zorder=2)
        ax.add_patch(Circle(head, 0.1, fc=LIPID['PC'][0], ec=INK, lw=0.4, zorder=3))


def fig_oil_surfactant():
    fig, axs = plt.subplots(1, 4, figsize=(15, 3.9))
    rng = np.random.default_rng(3)
    for ax in axs:
        frame_off(ax); ax.set_aspect('equal'); ax.set_xlim(0, 6); ax.set_ylim(0, 4.6)
    a, b, c, d = axs
    a.add_patch(Rectangle((0, 0), 6, 4.2, fc=WATER, alpha=0.35, lw=0))
    for _ in range(9):
        a.add_patch(Circle((rng.uniform(0.5, 2.5), rng.uniform(0.5, 3.7)), rng.uniform(0.12, 0.3), fc='#e9c46a', ec='#b8902f', lw=0.6))
    arrow(a, (2.8, 2.1), (3.6, 2.1))
    a.add_patch(Circle((4.7, 2.1), 0.85, fc='#e9c46a', ec='#b8902f', lw=0.8))
    a.set_title('Óleo na água: as gotas se juntam', fontsize=10)
    b.add_patch(Rectangle((0, 0), 6, 4.2, fc=WATER, alpha=0.35, lw=0))
    for x in np.arange(0.5, 6, 0.62):
        _lip(b, x, 2.95, True, LIPID['PC'][0], tlen=0.75, r=0.2, lw=1.1)
        _lip(b, x, 1.25, False, LIPID['PC'][0], tlen=0.75, r=0.2, lw=1.1)
    b.set_title('Lipídios: as caudas se escondem entre\nduas camadas de cabeças (bicamada)', fontsize=10)
    for ax, R, n, title in [(c, 2.0, 22, 'Com surfactante: o alvéolo\nfica aberto'),
                            (d, 1.05, 4, 'Prematuro, sem surfactante:\no alvéolo se fecha')]:
        _alveolus(ax, 3.0, 2.2, R, n)
        ax.set_ylim(-0.7, 4.9)
        ax.set_title(title, fontsize=10)
    for ang in np.radians([90, 210, 330]):
        d.annotate('', xy=(3.0 + 0.75 * np.cos(ang), 2.2 + 0.75 * np.sin(ang)),
                   xytext=(3.0 + 1.75 * np.cos(ang), 2.2 + 1.75 * np.sin(ang)),
                   arrowprops=dict(arrowstyle='-|>', color='#c44e52', lw=1.6))
    d.text(3.0, -0.45, 'a tensão da água puxa a parede para dentro', ha='center', fontsize=8.5, color='#c44e52')
    fig.suptitle('Efeito hidrofóbico: da gota de óleo ao surfactante pulmonar (esquema)', fontsize=11, y=1.04)
    save(fig, 'fund_oleo_surfactante.png')


# ---------------------------------------------------------------- 7. replicas divergem
def fig_replica_divergence(n_steps=60000, dt=1e-3, kT=0.25, seeds=(1, 4, 7), n_hist=300):
    """Langevin superamortecido em U(x) = (x^2 - 1)^2: mesmo ponto de partida, sementes diferentes."""
    def run(seed, n):
        rng = np.random.default_rng(seed)
        x = np.empty(n); x[0] = -1.0
        noise = rng.normal(0, np.sqrt(2 * kT * dt), n)
        for i in range(1, n):
            xi = x[i - 1]
            x[i] = xi - 4 * xi * (xi * xi - 1) * dt + noise[i]
        return x
    fig, (a, b) = plt.subplots(1, 2, figsize=(13, 3.8), gridspec_kw={'width_ratios': [2.2, 1]})
    t = np.arange(n_steps) * dt
    for k, s in enumerate(seeds):
        x = np.convolve(run(s, n_steps), np.ones(400) / 400, mode='same')[::50]   # media movel, 1 ponto a cada 50
        a.plot(t[::50], x, color=C[0], lw=1.4, ls=LINESTYLES[k], label=f'réplica {k + 1}')
    a.axhline(-1, color=BAND, lw=1); a.axhline(1, color=BAND, lw=1)
    a.set_xlabel('Tempo (unidades reduzidas)'); a.set_ylabel('Coordenada')
    a.set_yticks([-1, 0, 1]); a.set_yticklabels(['vale A', 'barreira', 'vale B'])
    a.legend(frameon=False, loc='lower right', bbox_to_anchor=(1.0, 1.0), ncol=3, fontsize=8.5)
    a.set_title('Mesmo ponto de partida, sementes diferentes', loc='left')
    first = []
    for s in range(100, 100 + n_hist):
        x = run(s, n_steps)
        hit = np.argmax(x > 0.9)
        if x[hit] > 0.9:
            first.append(hit * dt)
    b.hist(first, bins=20, color=C[0], alpha=0.6, ec='white')
    b.set_xlabel('Tempo da primeira passagem para o vale B'); b.set_ylabel(f'Réplicas (de {n_hist})')
    b.set_title('A passagem é um evento aleatório')
    fig.suptitle('Por que réplicas: a mesma física, trajetórias diferentes (modelo de brinquedo, Langevin em poço duplo)',
                 fontsize=11, y=1.03)
    fig.tight_layout()
    save(fig, 'fund_replicas_divergem.png')


# ---------------------------------------------------------------- p53: estrutura real em esferas
PDB_CACHE = OUT.parent / 'runs' / 'pdb'
P53_PROT = '#c9c2b4'
P53_DNA = '#8fa6c4'
CONTACT, STRUCT, Y220, LIGAND = '#c44e52', '#dd8452', '#8172b2', '#4c985f'


def fetch_pdb(code):
    """Baixa o PDB do RCSB uma vez e guarda em runs/pdb (fora do repositorio)."""
    import urllib.request
    PDB_CACHE.mkdir(parents=True, exist_ok=True)
    path = PDB_CACHE / f'{code}.pdb'
    if not path.exists():
        urllib.request.urlretrieve(f'https://files.rcsb.org/download/{code}.pdb', path)
    return path


def _shade(color, back):
    """Clareia o que esta ao fundo (back = 0 frente, 1 fundo), como em make_replica_figures."""
    from matplotlib.colors import to_rgb
    c = np.array(to_rgb(color))
    return tuple(c + (1 - c) * 0.55 * back)


def draw_spheres(ax, P, colors, radii):
    """Projecao ortografica: P em nm com colunas (tela x, tela y, profundidade); fundo desenhado primeiro."""
    from matplotlib.collections import PatchCollection
    d = P[:, 2]
    dn = (d - d.min()) / max(np.ptp(d), 1e-6)
    order = np.argsort(dn)
    patches = [Circle((P[k, 0], P[k, 1]), radii[k]) for k in order]
    faces = [_shade(colors[k], 1 - dn[k]) for k in order]
    ax.add_collection(PatchCollection(patches, facecolors=faces, edgecolors=[_shade(INK, 0.6)] * len(order),
                                      linewidths=0.25, match_original=False))
    ax.set_aspect('equal'); frame_off(ax)
    ax.set_xlim(P[:, 0].min() - 0.4, P[:, 0].max() + 0.4); ax.set_ylim(P[:, 1].min() - 0.4, P[:, 1].max() + 0.4)


def _frame(xyz, axis_x, toward_y):
    """Base ortonormal: x ao longo de axis_x, y o mais proximo de toward_y, z a profundidade."""
    ex = axis_x / np.linalg.norm(axis_x)
    ey = toward_y - ex * (toward_y @ ex); ey /= np.linalg.norm(ey)
    ez = np.cross(ex, ey)
    return np.c_[xyz @ ex, xyz @ ey, xyz @ ez]


def _label(ax, xy, text, xytext, color):
    ax.annotate(text, xy=xy, xytext=xytext, fontsize=8.5, color=INK, ha='center', va='center',
                arrowprops=dict(arrowstyle='-', color=color, lw=1.0, shrinkA=2, shrinkB=2),
                bbox=dict(boxstyle='round,pad=0.2', fc='white', ec=color, lw=0.8))


def fig_p53_structure():
    import mdtraj as md
    t = md.load_pdb(str(fetch_pdb('1TSR')))
    top = t.topology
    prot = top.select('chainid 3 and not element H')            # cadeia B, ligada ao sitio consenso
    dna = top.select('(chainid 0 or chainid 1) and not element H')
    xyz = t.xyz[0]
    hot = {175: STRUCT, 245: STRUCT, 249: STRUCT, 282: STRUCT, 248: CONTACT, 273: CONTACT, 220: Y220}
    names = {175: 'R175', 245: 'G245', 249: 'R249', 282: 'R282', 248: 'R248', 273: 'R273', 220: 'Y220'}
    sel = np.r_[prot, dna]
    cols, rad = [], []
    for i in sel:
        a = top.atom(i)
        rn = a.residue.resSeq
        if i in set(prot) and rn in hot:
            cols.append(hot[rn]); rad.append(0.2)
        else:
            cols.append(P53_PROT if i in set(prot) else P53_DNA); rad.append(0.17)
    # DNA na horizontal, proteina acima
    dc = xyz[dna] - xyz[dna].mean(0)
    axis = np.linalg.svd(dc, full_matrices=False)[2][0]
    up = xyz[prot].mean(0) - xyz[dna].mean(0)
    P = _frame(xyz[sel] - xyz[sel].mean(0), axis, up)
    fig, ax = plt.subplots(figsize=(9.5, 7.2))
    draw_spheres(ax, P, cols, rad)
    sel_index = {a: k for k, a in enumerate(sel)}
    spots = {}
    for rn in hot:
        idx = [sel_index[i] for i in prot if top.atom(i).residue.resSeq == rn]
        spots[rn] = P[idx, :2].mean(0)
    hl = [k for k, i in enumerate(sel) if i in set(prot) and top.atom(i).residue.resSeq in hot]
    draw_spheres_top = [(P[k, 0], P[k, 1], cols[k]) for k in hl]
    for x, y, col in draw_spheres_top:          # destacados por cima, para nao sumirem atras da proteina
        ax.add_patch(Circle((x, y), 0.2, fc=col, ec=_shade(INK, 0.6), lw=0.25, zorder=3))
    xmin, xmax = P[:, 0].min(), P[:, 0].max()
    mid = (xmin + xmax) / 2
    for side, sign in ((xmin - 1.3, -1), (xmax + 1.3, 1)):
        group = sorted([r for r in spots if (spots[r][0] < mid) == (sign < 0)], key=lambda r: spots[r][1])
        ys, last = [], -1e9
        for rn in group:                         # espaco minimo de 0,8 nm entre rotulos
            y = max(spots[rn][1], last + 0.8); ys.append(y); last = y
        for rn, y in zip(group, ys):
            _label(ax, spots[rn], names[rn], (side, y), hot[rn])
    ax.set_xlim(xmin - 2.0, xmax + 2.0)
    from matplotlib.lines import Line2D
    hs = [Line2D([], [], ls='none', marker='o', ms=9, mfc=c, mec=INK, mew=0.5) for c in (CONTACT, STRUCT, Y220, P53_PROT, P53_DNA)]
    ax.legend(hs, ['contato com o DNA', 'estrutural', 'Y220 (cavidade quando vira Cys)', 'domínio de ligação ao DNA', 'DNA'],
              frameon=False, loc='lower center', ncol=3, fontsize=8.5, bbox_to_anchor=(0.5, -0.12))
    ax.set_title('p53 ligada ao DNA (PDB 1TSR, átomos pesados) e os resíduos mais mutados em câncer', fontsize=11)
    save(fig, 'fund_p53_estrutura.png')


def fig_p53_y220c():
    import mdtraj as md
    t = md.load_pdb(str(fetch_pdb('2VUK')))
    top = t.topology
    lig = top.select('resname P83 and chainid 3')
    if not len(lig):
        lig = top.select('resname P83')
    chain = top.atom(lig[0]).residue.chain.index
    prot_chain = [c.index for c in top.chains if c.chain_id == top.atom(lig[0]).residue.chain.chain_id and c.n_residues > 50][0]
    prot = top.select(f'chainid {prot_chain} and protein and not element H')
    xyz = t.xyz[0]
    sel = np.r_[prot, lig]
    cols = [Y220 if (i in set(prot) and top.atom(i).residue.resSeq == 220) else (LIGAND if i in set(lig) else P53_PROT) for i in sel]
    rad = [0.2 if c != P53_PROT else 0.17 for c in cols]
    view = xyz[lig].mean(0) - xyz[prot].mean(0)                  # olhar de fora para dentro da cavidade
    ref = np.array([0, 0, 1.0]) if abs(view[2]) < 0.9 * np.linalg.norm(view) else np.array([1.0, 0, 0])
    ex = np.cross(ref, view); ex /= np.linalg.norm(ex)
    ey = np.cross(view / np.linalg.norm(view), ex)
    c0 = xyz[sel].mean(0)
    P = np.c_[(xyz[sel] - c0) @ ex, (xyz[sel] - c0) @ ey, (xyz[sel] - c0) @ (view / np.linalg.norm(view))]
    fig, (a, b) = plt.subplots(1, 2, figsize=(12, 5.6))
    draw_spheres(a, P, cols, rad)
    a.set_title('Y220C com o ligante PhiKan083 na cavidade\n(PDB 2VUK, vista de fora para dentro)', fontsize=10)
    # recorte: so atomos a menos de 1,2 nm do ligante
    dmin = np.min(np.linalg.norm(xyz[sel][:, None, :] - xyz[lig][None, :, :], axis=2), axis=1)
    near = dmin < 1.2
    draw_spheres(b, P[near], [c for c, k in zip(cols, near) if k], [r for r, k in zip(rad, near) if k])
    b.set_title('Detalhe: átomos a menos de 1,2 nm do ligante', fontsize=10)
    from matplotlib.lines import Line2D
    hs = [Line2D([], [], ls='none', marker='o', ms=9, mfc=c, mec=INK, mew=0.5) for c in (Y220, LIGAND, P53_PROT)]
    fig.legend(hs, ['Cys220 (no selvagem, Tyr)', 'PhiKan083', 'p53'], frameon=False, loc='lower center', ncol=3, fontsize=8.5,
               bbox_to_anchor=(0.5, -0.02))
    save(fig, 'fund_p53_y220c.png')
    _ = chain


# ---------------------------------------------------------------- quatro efeitos de uma mutacao pontual (esquema em esferas)
def _chain(ax, pts, color=P53_PROT, r=0.28, z=2):
    for k in range(len(pts) - 1):
        ax.plot(*np.c_[pts[k], pts[k + 1]], color=GREY, lw=1.2, zorder=z - 1)
    for p in pts:
        ax.add_patch(Circle(p, r, fc=color, ec=_shade(INK, 0.6), lw=0.5, zorder=z))


def _hairpin(x0=0.0, y0=0.0, n=7, gap=1.1, step=0.62):
    up = [(x0, y0 + step * k) for k in range(n)]
    down = [(x0 + gap, y0 + step * (n - 1 - k)) for k in range(n)]
    return up + [(x0 + gap / 2, y0 + step * n)] + down


def fig_mutation_effects():
    fig, axs = plt.subplots(1, 4, figsize=(15.5, 4.6))
    for ax in axs:
        frame_off(ax); ax.set_aspect('equal'); ax.set_xlim(-1.4, 4.2); ax.set_ylim(-1.2, 5.6)
    a, b, c, d = axs
    # estrutura: cadeia lateral grande preenche o bolsao; a menor deixa um buraco
    for x0, big, lab in ((-0.6, True, 'selvagem: Tyr'), (1.9, False, 'mutante: Cys')):
        pts = _hairpin(x0, 0.4, n=6, gap=1.6)
        _chain(a, pts)
        a.add_patch(Circle((x0 + 0.8, 2.0), 0.55 if big else 0.25, fc=Y220, ec=INK, lw=0.5, zorder=3))
        if not big:
            a.add_patch(Circle((x0 + 0.8, 2.6), 0.42, fc='white', ec=Y220, lw=1, ls='--', zorder=3))
        a.text(x0 + 0.8, -0.6, lab, ha='center', fontsize=8.5)
    a.set_title('Estrutura: abre uma cavidade', fontsize=10)
    # enovelamento: dois pocos
    xx = np.linspace(0, 3.6, 300)
    for depth, col, ls, lab in ((3.0, C[0], '-', 'selvagem'), (1.2, C[7], (0, (6, 3)), 'mutante')):
        yy = 3.6 - depth * np.exp(-((xx - 1.0) / 0.45) ** 2) + 0.4 * np.exp(-((xx - 2.6) / 0.6) ** 2)
        b.plot(xx, yy, color=col, lw=2, ls=ls, label=lab)
    b.text(1.0, -0.3, 'enovelada', ha='center', fontsize=8.5); b.text(2.6, 4.45, 'desenovelada', ha='center', fontsize=8.5)
    b.legend(frameon=False, loc='upper left', fontsize=8.5)
    b.set_title('Enovelamento: o vale nativo fica raso', fontsize=10)
    # flexibilidade: quadros sobrepostos de uma alca
    rng = np.random.default_rng(2)
    for x0, amp, lab in ((-0.6, 0.08, 'selvagem'), (1.9, 0.35, 'mutante')):
        base = np.array(_hairpin(x0, 0.4, n=6, gap=1.2))
        for f in range(6):
            w = base + rng.normal(0, amp, base.shape) * (base[:, 1:2] / 4.0)
            c.plot(w[:, 0], w[:, 1], color=C[0], lw=0.9, alpha=0.45)
        _chain(c, base[[0, -1]], r=0.22)
        c.text(x0 + 0.6, -0.6, lab, ha='center', fontsize=8.5)
    c.set_title('Dinâmica: a alça oscila mais', fontsize=10)
    # funcao: arginina toca o DNA
    for x0, touch, lab in ((-0.6, True, 'selvagem: Arg248'), (1.9, False, 'mutante: Gln248')):
        for k in range(5):
            d.add_patch(Circle((x0 + 0.25 * (k % 2), 0.2 + 0.55 * k), 0.24, fc=P53_DNA, ec=_shade(INK, 0.6), lw=0.5))
        pts = [(x0 + 1.6, 4.6), (x0 + 1.6, 3.9), (x0 + 1.5, 3.2)]
        _chain(d, pts)
        tip = (x0 + 0.75, 2.5) if touch else (x0 + 1.2, 2.8)
        d.plot([pts[-1][0], tip[0]], [pts[-1][1], tip[1]], color=GREY, lw=1.2)
        d.add_patch(Circle(tip, 0.28, fc=CONTACT if touch else STRUCT, ec=INK, lw=0.5, zorder=3))
        if touch:
            d.plot([tip[0] - 0.3, x0 + 0.45], [tip[1], 2.4], color=CONTACT, lw=1, ls=':')
        d.text(x0 + 0.9, -0.6, lab, ha='center', fontsize=8.5)
    d.set_title('Função: perde o contato com o DNA', fontsize=10)
    fig.suptitle('Quatro efeitos de uma mutação pontual (esquema)', fontsize=11, y=1.0)
    save(fig, 'fund_mutacao_efeitos.png')


# ---------------------------------------------------------------- 11. termos do campo de forca e diedro
ATOM, OXY, HYD = '#c9c2b4', '#c44e52', '#ffffff'


def _atom(ax, xy, r=0.32, fc=ATOM, text=None, z=4, fs=9, lw=0.8, alpha=1.0):
    ax.add_patch(Circle(xy, r, fc=fc, ec=INK, lw=lw, zorder=z, alpha=alpha))
    if text:
        ax.text(xy[0], xy[1], text, ha='center', va='center', fontsize=fs, zorder=z + 1)


def _stick(ax, p, q, color=GREY, lw=2.2, z=2, ls='-'):
    ax.plot([p[0], q[0]], [p[1], q[1]], color=color, lw=lw, zorder=z, ls=ls, solid_capstyle='round')


def _spring(ax, p, q, n=7, amp=0.16, color=GREY):
    p, q = np.array(p, float), np.array(q, float)
    t = np.linspace(0, 1, 2 * n + 3)
    d = q - p; nrm = np.array([-d[1], d[0]]) / np.linalg.norm(d)
    off = np.zeros_like(t); off[2:-2] = amp * (-1) ** np.arange(len(t) - 4)
    pts = p + np.outer(t, d) + np.outer(off, nrm)
    ax.plot(pts[:, 0], pts[:, 1], color=color, lw=1.6, zorder=2)


def _dim(ax, p, q, text, y_text=None, color=INK):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle='<|-|>', mutation_scale=10, color=color, lw=1.0))
    ax.text((p[0] + q[0]) / 2, y_text if y_text is not None else p[1] - 0.25, text,
            ha='center', va='top', fontsize=9)


def fig_force_field_schemes():
    fig, axs = plt.subplots(1, 5, figsize=(17, 3.6))
    for ax in axs:
        frame_off(ax); ax.set_aspect('equal'); ax.set_xlim(-1.6, 3.6); ax.set_ylim(-1.7, 2.7)
    a, b, c, d, e = axs
    # ligacao: mola entre dois atomos
    _spring(a, (0, 0.5), (2, 0.5))
    _atom(a, (0, 0.5), 0.42, text='C'); _atom(a, (2, 0.5), 0.42, text='C')
    _dim(a, (0, -0.35), (2, -0.35), 'r', y_text=-0.55)
    a.text(1, 1.5, 'esticar ou comprimir\ncusta energia', ha='center', va='bottom', fontsize=8.5)
    a.set_title('Ligação: dois átomos', fontsize=10)
    # angulo: tres atomos, theta no atomo central
    o = np.array([1.0, -0.6]); L = 1.8; h = np.radians(109.5 / 2)
    p1, p3 = o + L * np.array([-np.sin(h), np.cos(h)]), o + L * np.array([np.sin(h), np.cos(h)])
    _stick(b, o, p1); _stick(b, o, p3)
    th = np.linspace(90 - 54.75, 90 + 54.75, 60)
    b.plot(o[0] + 0.85 * np.cos(np.radians(th)), o[1] + 0.85 * np.sin(np.radians(th)), color=C[0], lw=1.4)
    b.text(o[0], o[1] + 1.05, 'θ', ha='center', va='bottom', fontsize=11, color=C[0])
    for pt, k in ((p1, '1'), (o, '2'), (p3, '3')):
        _atom(b, pt, 0.36, text=k)
    b.text(1.0, -1.35, 'θ₀ = 109,5° (carbono sp³)', ha='center', va='center', fontsize=8.5)
    b.set_title('Ângulo: três átomos', fontsize=10)
    # diedro: quatro atomos, rotacao em torno da ligacao 2-3
    p1, p2, p3, p4 = (-0.6, 1.3), (0.0, 0.3), (1.6, 0.3), (2.2, -0.7)
    _stick(c, p1, p2); _stick(c, p2, p3, lw=3.0); _stick(c, p3, p4)
    _stick(c, p3, (2.2, 1.3), ls=(0, (2, 2)), lw=1.4)
    _atom(c, (2.2, 1.3), 0.36, fc='white', alpha=0.6)
    for pt, k in ((p1, '1'), (p2, '2'), (p3, '3'), (p4, '4')):
        _atom(c, pt, 0.36, text=k)
    th = np.radians(np.linspace(110, 410, 120))
    cx, rx, ry = 0.8, 0.28, 0.6
    c.plot(cx + rx * np.cos(th[:-6]), 0.3 + ry * np.sin(th[:-6]), color=C[0], lw=1.4, zorder=5)
    arrow(c, (cx + rx * np.cos(th[-6]), 0.3 + ry * np.sin(th[-6])),
          (cx + rx * np.cos(th[-1]), 0.3 + ry * np.sin(th[-1])), color=C[0], ms=10, z=5)
    c.text(cx, 1.15, 'φ', ha='center', va='bottom', fontsize=11, color=C[0])
    c.text(1.0, -1.35, 'giro em torno da ligação 2-3', ha='center', va='center', fontsize=8.5)
    c.set_title('Diedro: quatro átomos', fontsize=10)
    # Lennard-Jones: tres distancias entre dois oxigenios
    for y, gap, lab, push in ((2.0, 0.42, 'perto demais: repele', 1), (0.65, 0.72, 'r ≈ 0,35 nm: mínimo', 0),
                              (-0.7, 2.0, 'longe: U ≈ 0', 0)):
        x0 = 1.0 - gap / 2
        _atom(d, (x0, y), 0.34, fc=OXY, alpha=0.85); _atom(d, (x0 + gap, y), 0.34, fc=OXY, alpha=0.85)
        if push:
            arrow(d, (x0 - 0.38, y), (x0 - 0.9, y), color=C[0], ms=9)
            arrow(d, (x0 + gap + 0.38, y), (x0 + gap + 0.9, y), color=C[0], ms=9)
        d.text(1.0, y - 0.42, lab, ha='center', va='top', fontsize=8.5)
    d.set_title('Lennard-Jones: O e O da água', fontsize=10)
    # Coulomb: atracao entre cargas opostas, repulsao entre iguais
    for y, qa, qb, fa, attract, lab in ((1.3, 'O⁻', 'H⁺', OXY, True, 'atrai'), (-0.6, 'H⁺', 'H⁺', HYD, False, 'repele')):
        _atom(e, (-0.6, y), 0.42 if fa == OXY else 0.3, fc=fa, text=qa, alpha=0.9)
        _atom(e, (1.6, y), 0.3, fc=HYD, text=qb)
        if attract:
            arrow(e, (-0.1, y), (0.4, y), color=C[0], ms=9); arrow(e, (1.2, y), (0.7, y), color=C[0], ms=9)
        else:
            arrow(e, (-1.0, y), (-1.5, y), color=C[0], ms=9); arrow(e, (2.0, y), (2.5, y), color=C[0], ms=9)
        e.text(0.5, y - 0.6, lab, ha='center', va='top', fontsize=8.5)
    e.set_title('Coulomb: cargas da água', fontsize=10)
    fig.suptitle('O que cada termo do campo de força descreve (esquema)', fontsize=11, y=1.02)
    save(fig, 'fund_termos_esquema.png')


def _newman(ax, phi, title, sub_col=C[0], back_col=C[5]):
    """Projecao de Newman: olhando a ligacao 2-3 de frente; atomo 1 na frente, atomo 4 atras."""
    ax.set_aspect('equal', adjustable='datalim'); frame_off(ax); ax.set_xlim(-1.75, 1.75); ax.set_ylim(-1.75, 1.95)
    R, L = 0.55, 1.25
    shift = 8 if phi % 360 == 0 else 0           # eclipsada: afasta um pouco para os dois aparecerem
    for k in range(3):                           # atras: comecam na borda do disco
        ang = np.radians(90 - phi - shift + 120 * k)
        u = np.array([np.cos(ang), np.sin(ang)])
        _stick(ax, R * u, L * u, lw=2.0, z=1)
        _atom(ax, L * u, 0.24 if k else 0.3, fc=back_col if k == 0 else HYD, text='4' if k == 0 else None, z=1.5)
    ax.add_patch(Circle((0, 0), R, fc='white', ec=INK, lw=1.2, zorder=3))
    ax.text(0, -0.02, '2', ha='center', va='center', fontsize=8, color=GREY, zorder=6)
    for k in range(3):                           # frente: saem do centro
        ang = np.radians(90 + 120 * k)
        u = np.array([np.cos(ang), np.sin(ang)])
        _stick(ax, 0.12 * u, L * u, lw=2.0, z=4)
        _atom(ax, L * u, 0.3 if k == 0 else 0.24, fc=sub_col if k == 0 else HYD, text='1' if k == 0 else None, z=5)
    if phi % 360:
        th = np.radians(np.linspace(90, 90 - phi, 60))
        ax.plot(0.85 * np.cos(th), 0.85 * np.sin(th), color=INK, lw=1.0, zorder=6)
        m = np.radians(90 - phi / 2)
        ax.text(1.05 * np.cos(m), 1.05 * np.sin(m), 'φ', ha='center', va='center', fontsize=10, zorder=6)
    ax.set_title(title, fontsize=10)


def fig_dihedral():
    fig = plt.figure(figsize=(18, 4.4))
    gs = fig.add_gridspec(1, 5, width_ratios=[1.35, 1, 1, 1, 1.35], wspace=0.25)
    # painel 1: quatro atomos em 3D, dois planos e o angulo entre eles
    ax = fig.add_subplot(gs[0]); ax.set_aspect('equal', adjustable='datalim'); frame_off(ax)
    phi = np.radians(120)
    b = 2.2; s = 1.4
    X = np.array([[-0.5, s, 0.0], [0, 0, 0], [b, 0, 0], [b + 0.5, s * np.cos(phi), s * np.sin(phi)]])
    def Rx(a):
        a = np.radians(a); return np.array([[1, 0, 0], [0, np.cos(a), -np.sin(a)], [0, np.sin(a), np.cos(a)]])
    def Ry(a):
        a = np.radians(a); return np.array([[np.cos(a), 0, np.sin(a)], [0, 1, 0], [-np.sin(a), 0, np.cos(a)]])
    M = Rx(25) @ Ry(-40)
    P = lambda v: (M @ np.asarray(v, float))[:2]
    v1, v4 = np.array([0, s, 0]), np.array([0, s * np.cos(phi), s * np.sin(phi)])
    for v, col in ((v1, C[0]), (v4, C[5])):
        quad = [P(X[1] + [-0.4, 0, 0]), P(X[2] + [0.4, 0, 0]), P(X[2] + [0.4, 0, 0] + 1.15 * v), P(X[1] + [-0.4, 0, 0] + 1.15 * v)]
        ax.add_patch(Polygon(quad, closed=True, fc=col, alpha=0.14, ec=col, lw=0.8, zorder=1))
    mid = (X[1] + X[2]) / 2
    arc = np.array([P(mid + 0.8 * (np.cos(t) * v1 / s + np.sin(t) * np.array([0, 0, 1]))) for t in np.linspace(0, phi, 50)])
    ax.plot(arc[:, 0], arc[:, 1], color=INK, lw=1.0, zorder=6)
    for k, v in ((0, v1), (1, v4)):
        _stick(ax, P(mid), P(mid + 0.8 * v / s), color=INK, lw=0.6, ls=(0, (2, 2)), z=5)
    k = np.argmin(arc[:, 0])
    ax.text(arc[k, 0] - 0.12, arc[k, 1], 'φ', ha='right', va='center', fontsize=11, zorder=6)
    pts = [P(x) for x in X]
    _stick(ax, pts[0], pts[1]); _stick(ax, pts[1], pts[2], lw=3.0); _stick(ax, pts[2], pts[3])
    for k, (pt, col) in enumerate(zip(pts, (C[0], ATOM, ATOM, C[5]))):
        _atom(ax, pt, 0.22, fc=col, text=str(k + 1), fs=8.5, z=7)
    allp = np.array([P(X[i] + [d, 0, 0] + 1.15 * v) for v in (v1, v4) for i, d in ((1, -0.4), (2, 0.4))] + pts)
    c1 = P(X[2] + [0.4, 0, 0] + 1.15 * v1); c4 = P(X[2] + [0.4, 0, 0])
    ax.text(c1[0] + 0.1, c1[1], 'plano 1-2-3', color=C[0], fontsize=8.5, ha='left', va='center')
    ax.text(c4[0] + 0.1, c4[1], 'plano 2-3-4', color=C[5], fontsize=8.5, ha='left', va='center')
    x0, x1 = allp[:, 0].min() - 0.3, max(c1[0], c4[0]) + 1.6
    y0, y1 = allp[:, 1].min() - 0.4, allp[:, 1].max() + 0.4
    bb = ax.get_position(); cell = (bb.width * fig.get_figwidth()) / (bb.height * fig.get_figheight())
    if (x1 - x0) / (y1 - y0) < cell:             # alarga o eixo curto para o painel ocupar a celula inteira
        pad = (cell * (y1 - y0) - (x1 - x0)) / 2; x0, x1 = x0 - pad, x1 + pad
    else:
        pad = ((x1 - x0) / cell - (y1 - y0)) / 2; y0, y1 = y0 - pad, y1 + pad
    ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
    ax.set_title('Ângulo entre dois planos', fontsize=10)
    # paineis 2-4: projecoes de Newman
    for k, (ph, t) in enumerate(((180, 'φ = 180°: anti, mínimo'), (60, 'φ = 60°: gauche, mínimo'),
                                  (0, 'φ = 0°: eclipsada, máximo'))):
        _newman(fig.add_subplot(gs[k + 1]), ph, t)
    # painel 5: energia do diedro com os tres estados marcados
    ax = fig.add_subplot(gs[4])
    x = np.linspace(-180, 180, 400)
    U = lambda f: 5 * (1 + np.cos(np.radians(3 * np.asarray(f))))
    ax.plot(x, U(x), color=C[0])
    for f, lab, dy, ha in ((180, 'anti', -0.7, 'right'), (60, 'gauche', -0.7, 'center'), (0, 'eclipsada', 0.6, 'center')):
        ax.plot(f, U(f), 'o', color=C[5], ms=6, zorder=5)
        ax.text(f, U(f) + dy, lab, ha=ha, fontsize=8.5, va='bottom' if dy > 0 else 'top')
    ax.set(xlabel='φ (graus)', ylabel='U (kJ/mol)', xticks=[-180, -60, 60, 180], ylim=(-2.6, 12.8))
    ax.set_title('U(φ) = k(1 + cos 3φ)', fontsize=10)
    fig.suptitle('O diedro: rotação em torno da ligação central (esquema)', fontsize=11, y=1.03)
    save(fig, 'fund_diedro.png')


MEM_REP = Path(__file__).resolve().parents[1] / 'runs' / 'aula_martini' / 'rep1'
DENS_COL = {'NC3': '#4c72b0', 'PO4': '#c44e52', 'GL': '#55a868', 'Caudas': '#8c8c8c', 'W': WATER}


def _density_profile():
    """Perfil real da replica 1 (segunda metade da equilibracao), como no notebook,
    e a profundidade media |z| de cada bead do DSPC para posicionar o esquema."""
    import mdtraj as md
    t = md.load(MEM_REP / 'mem_eq_whole.xtc', top=MEM_REP / 'mem_rot.gro')
    last = t[t.n_frames // 2:]
    tq = last.topology
    tails = tq.select('name =~ "C[1-4][AB]"')
    groups = {'NC3': tq.select('name NC3'), 'PO4': tq.select('name PO4'), 'GL': tq.select('name =~ "GL[12]"'),
              'Caudas': tails, 'W': tq.select('resname W')}
    Lz = last.unitcell_lengths[:, 2]
    z = last.xyz[:, :, 2]
    ang = 2 * np.pi * z[:, tails] / Lz[:, None]
    center = (np.arctan2(np.sin(ang).mean(1), np.cos(ang).mean(1)) / (2 * np.pi)) * Lz
    zc = (z - center[:, None] + Lz[:, None] / 2) % Lz[:, None] - Lz[:, None] / 2
    bins = np.linspace(-Lz.mean() / 2, Lz.mean() / 2, 81)
    mid = 0.5 * (bins[1:] + bins[:-1])
    prof = {k: np.histogram(zc[:, i].ravel(), bins=bins)[0] / last.n_frames for k, i in groups.items()}
    depth = {n: float(np.abs(zc[:, tq.select(f'name {n}')]).mean())
             for n in ('NC3', 'PO4', 'GL1', 'GL2', 'C1A', 'C2A', 'C3A', 'C4A', 'C1B', 'C2B', 'C3B', 'C4B')}
    return mid, prof, depth


def fig_density_scheme():
    """Bicamada de DSPC em beads deitada sobre o eixo do perfil de densidade: cada curva
    do grafico e a contagem, fatia a fatia, de um tipo de bead do esquema de cima."""
    mid, prof, dep = _density_profile()
    lim = mid[-1] + 0.1
    H = 5.8                                            # altura do esquema (nm), mesma escala do eixo x
    left, wid, bot, hb, gap, top = 0.85, 6.6, 0.7, 3.4, 0.75, 0.8
    ht = wid * H / (2 * lim)
    W, Ht = left + wid + 1.3, bot + hb + gap + ht + top
    fig = plt.figure(figsize=(W, Ht))
    b = fig.add_axes([left / W, bot / Ht, wid / W, hb / Ht])
    a = fig.add_axes([left / W, (bot + hb + gap) / Ht, wid / W, ht / Ht], sharex=b)
    rng = np.random.default_rng(3)
    r = 0.17
    rows = np.arange(0.45, H, 0.96)                    # posicao de cada lipidio ao longo do plano
    tail_y = (-0.21, 0.21)                             # cadeias A e B lado a lado
    occupied = []
    for side in (-1, 1):
        for y in rows:
            j = rng.normal(0, 0.07)                    # flutuacao de cada lipidio ao longo da normal
            X = lambda n: side * (dep[n] + j)  # noqa: E731
            pa = [(X('GL1'), y + tail_y[0])] + [(X(f'C{k}A'), y + tail_y[0]) for k in range(1, 5)]
            pb = [(X('GL2'), y + tail_y[1])] + [(X(f'C{k}B'), y + tail_y[1]) for k in range(1, 5)]
            po4 = (X('PO4'), y)
            nc3 = (X('NC3'), y + 0.36)                 # a colina deita ao lado do fosfato
            for p, q in [(nc3, po4), (po4, pa[0]), (po4, pb[0]), (pa[0], pb[0])] + list(zip(pa, pa[1:])) + list(zip(pb, pb[1:])):
                a.plot(*np.c_[p, q], color=GREY, lw=1.0, zorder=1)
            beads = [(nc3, DENS_COL['NC3']), (po4, DENS_COL['PO4']), (pa[0], DENS_COL['GL']), (pb[0], DENS_COL['GL'])]
            beads += [(p, DENS_COL['Caudas']) for p in pa[1:] + pb[1:]]
            for p, col in beads:
                a.add_patch(Circle(p, r, fc=col, ec=INK, lw=0.5, zorder=3))
                occupied.append(p)
    # Agua sorteada pela densidade real de W: some no nucleo, entra pouco entre as cabecas
    w_prof = prof['W'] / prof['W'].max()
    occ = np.array(occupied)
    water = []
    for _ in range(40000):
        p = np.array([rng.uniform(-lim + r, lim - r), rng.uniform(r, H - r)])
        dens = np.interp(p[0], mid, w_prof)
        if dens < 0.02 or rng.random() > dens:
            continue
        pts = np.r_[occ, np.array(water).reshape(-1, 2)]
        if np.min(np.hypot(*(pts - p).T)) < 2 * r + 0.03:
            continue
        water.append(p)
    for p in water:
        a.add_patch(Circle(p, r, fc=WATER, ec='#6f9fb8', lw=0.5, zorder=2))
    a.set_xlim(-lim, lim); a.set_ylim(0, H)
    a.set_yticks([]); a.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    for x, lab in ((-3.7, 'água'), (-2.0, 'cabeças'), (0, 'núcleo de caudas'), (2.0, 'cabeças'), (3.7, 'água')):
        a.text(x, H + 0.12, lab, ha='center', va='bottom', fontsize=9, color=INK)
    a.set_title('Bicamada de DSPC vista de lado: a normal é o eixo horizontal (esquema)', fontsize=11, pad=22)

    for k, y in prof.items():
        b.plot(mid, y, color=DENS_COL[k], lw=1.6, label=k)
    i_l = np.argmax(np.where(mid < 0, prof['PO4'], 0)); i_r = np.argmax(np.where(mid > 0, prof['PO4'], 0))
    for ax in (a, b):
        for x in (mid[i_l], 0, mid[i_r]):
            ax.axvline(x, color=INK, lw=0.7, ls=(0, (3, 3)), zorder=0)
    ytop = 44
    b.add_patch(FancyArrowPatch((mid[i_l], ytop), (mid[i_r], ytop), arrowstyle='<|-|>', mutation_scale=10, color=DENS_COL['PO4'], lw=1.1))
    b.text(0, ytop + 1.2, f'$D_{{PP}}$ = {mid[i_r] - mid[i_l]:.2f} nm'.replace('.', ','), ha='center', va='bottom', fontsize=10, color=INK)
    b.set_ylim(-1, 50); b.set_xlim(-lim, lim)
    b.set_xlabel('Distância ao centro da bicamada (nm)'); b.set_ylabel('beads por fatia')
    b.set_title('Perfil de densidade ao longo da normal (réplica 1)', fontsize=11)
    fig.legend(*b.get_legend_handles_labels(), loc='center left', bbox_to_anchor=((left + wid + 0.15) / W, (bot + hb / 2) / Ht), frameon=False, fontsize=9)
    save(fig, 'martini_perfil_esquema.png')


def main():
    OUT.mkdir(exist_ok=True)
    fig_photo_film(); fig_sickle(); fig_misfolding(); fig_milestones()
    fig_ubiquitin_tag(); fig_oil_surfactant(); fig_replica_divergence()
    fig_p53_structure(); fig_p53_y220c(); fig_mutation_effects()
    fig_force_field_schemes(); fig_dihedral(); fig_density_scheme()


if __name__ == '__main__':
    main()
