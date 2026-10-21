"""Gera as ilustracoes do capitulo Fundamentos (imagens/fund_*.png).

Estilo de figura de artigo, o mesmo de viz/replica_plots.py. Esquemas
marcados como ilustrativos usam funcoes analiticas ou modelos de brinquedo;
as coordenadas moleculares sao reais (1UBQ).

Uso: python tools/make_figures.py --aa ubq_aa.gro --ua ubq_ua.gro --cg ubq_cg.pdb
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Rectangle

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'viz'))
from pure import palette, mpl as pure_mpl  # noqa: E402
from replica_plots import STYLE, INK, BAND  # noqa: E402
from replica_views import SEM as SEM_  # noqa: E402
sys.path.insert(0, str(ROOT / 'tools'))

C = palette.CATEGORICAL
OUT = ROOT / 'imagens'
DPI = 200
pure_mpl.use('light')          # fontes e colormaps pure
plt.rcParams.update(STYLE)
plt.rcParams.update({'font.size': 10})


def save(fig, name):
    fig.savefig(OUT / name, dpi=DPI, bbox_inches='tight')
    plt.close(fig)
    print('  ', name)


def frame_off(ax):
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)


# ---------------------------------------------------------------- 1. ciclo da MD
def _arc(ax, t0, t1, R, color, head=True):
    """Arco do anel entre os angulos t0 e t1 (graus), com ponta de seta no fim."""
    th = np.radians(np.linspace(t0, t1, 80))
    x, y = R * np.cos(th), R * np.sin(th)
    ax.plot(x[:-2], y[:-2], color=color, lw=2.2, solid_capstyle='round', zorder=2)
    if head:
        ax.add_patch(FancyArrowPatch((x[-6], y[-6]), (x[-1], y[-1]), arrowstyle='-|>',
                                     mutation_scale=18, color=color, lw=0, zorder=3))


def fig_md_loop():
    fig, ax = plt.subplots(figsize=(11, 6.8))
    frame_off(ax); ax.set_aspect('equal')
    ax.set_xlim(-7.4, 7.4); ax.set_ylim(-3.6, 3.5)
    R = 2.55
    ring = plt.Circle((0, 0), R, fill=False, ec=BAND, lw=7, alpha=0.35, zorder=1)
    ax.add_patch(ring)
    # passo: (angulo, numero, nome, equacao, cor)
    steps = [(90, '1', 'Posições', '$\\mathbf{r}(t)$', C[0]),
             (0, '2', 'Forças', '$\\mathbf{F}_i = -\\nabla_i\\, U(\\mathbf{r})$', C[1]),
             (-90, '3', 'Velocidades', '$\\mathbf{v}(t+\\frac{\\Delta t}{2}) = \\mathbf{v}(t-\\frac{\\Delta t}{2}) + \\frac{\\mathbf{F}}{m}\\Delta t$', C[2]),
             (180, '4', 'Novas posições', '$\\mathbf{r}(t+\\Delta t) = \\mathbf{r}(t) + \\mathbf{v}(t+\\frac{\\Delta t}{2})\\,\\Delta t$', C[3])]
    # folga angular de cada lado do no (caixas largas no topo e na base)
    gap = {90: 22, 0: 20, -90: 40, 180: 20}
    for i, (ang, *_rest) in enumerate(steps):
        nxt = steps[(i + 1) % 4][0]
        end = nxt if nxt < ang else nxt - 360
        _arc(ax, ang - gap[ang], end + gap[nxt], R, INK)
    fig.canvas.draw()
    inv = ax.transData.inverted()
    for ang, num, name, eq, col in steps:
        x, y = R * np.cos(np.radians(ang)), R * np.sin(np.radians(ang))
        t1 = ax.text(x, y + 0.24, name, ha='center', va='center', fontsize=11.5, color=INK,
                     weight='bold', zorder=5)
        t2 = ax.text(x, y - 0.3, eq, ha='center', va='center', fontsize=11, color=INK, zorder=5)
        r = fig.canvas.get_renderer()
        ext = [inv.transform(tx.get_window_extent(r)) for tx in (t1, t2)]
        if ang in (0, 180):        # cartoes laterais por fora do anel, borda interna sobre ele
            w = max(e[1, 0] - e[0, 0] for e in ext) / 2 + 0.28
            dx = (R - 0.3 + w - x) if ang == 0 else (-R + 0.3 - w - x)
            for tx in (t1, t2):
                tx.set_x(tx.get_position()[0] + dx)
            ext = [e + [dx, 0] for e in ext]
        x0 = min(e[0, 0] for e in ext) - 0.28; x1 = max(e[1, 0] for e in ext) + 0.28
        y0 = min(e[0, 1] for e in ext) - 0.2; y1 = max(e[1, 1] for e in ext) + 0.2
        ax.add_patch(FancyBboxPatch((x0, y0), x1 - x0, y1 - y0, boxstyle='round,pad=0,rounding_size=0.22',
                                    fc='white', ec=col, lw=1.6, zorder=4))
        ax.add_patch(plt.Circle((x0, y1), 0.2, fc=col, ec='white', lw=1.5, zorder=6))
        ax.text(x0, y1, num, ha='center', va='center', fontsize=9, color='white', weight='bold', zorder=7)
    # centro: uma particula avancando um passo, com a forca que a desvia
    s = np.linspace(-1.1, 1.1, 60)
    ax.plot(s, -0.25 * s ** 2 + 0.25, color=BAND, lw=1.2, ls=(0, (3, 3)), zorder=2)
    for xx, al in ((-0.9, 0.3), (0, 1.0), (0.9, 0.3)):
        ax.add_patch(plt.Circle((xx, -0.25 * xx ** 2 + 0.25), 0.17, fc=C[0], ec='white', lw=1, alpha=al, zorder=4))
    ax.add_patch(FancyArrowPatch((0, 0.25), (0.75, 0.25), arrowstyle='-|>', mutation_scale=12,
                                 color=C[3], lw=1.6, shrinkA=6, zorder=5))
    ax.add_patch(FancyArrowPatch((0, 0.25), (0, -0.45), arrowstyle='-|>', mutation_scale=12,
                                 color=C[1], lw=1.6, shrinkA=6, zorder=5))
    ax.text(0.42, 0.45, '$\\mathbf{v}$', color=C[3], fontsize=11, ha='center')
    ax.text(0.17, -0.5, '$\\mathbf{F}$', color=C[1], fontsize=11, ha='left')
    ax.text(-0.9, 0.62, '$t-\\Delta t$', color=INK, fontsize=8.5, ha='center', alpha=0.7)
    ax.text(0, 0.62, '$t$', color=INK, fontsize=8.5, ha='center')
    ax.text(0.9, 0.62, '$t+\\Delta t$', color=INK, fontsize=8.5, ha='center', alpha=0.7)
    ax.text(0, -1.05, 'repete a cada passo $\\Delta t$', ha='center', va='center', fontsize=10, color=INK)
    ax.text(0, -1.42, 'AA 2 fs  ·  Martini 20 fs', ha='center', va='center', fontsize=9, color=INK, alpha=0.7)
    ax.set_title('Ciclo de integração da dinâmica molecular (leap-frog)')
    save(fig, 'fund_ciclo_md.png')


# ---------------------------------------------------------------- 2. termos do campo de forca
def fig_force_field():
    fig, axes = plt.subplots(1, 5, figsize=(17, 3.4))
    r = np.linspace(0.08, 0.2, 200); b0, kb = 0.14, 3e5
    axes[0].plot(r, 0.5 * kb * (r - b0) ** 2, color=C[0])
    axes[0].set(title='Ligação', xlabel='r (nm)', ylabel='U (kJ/mol)', ylim=(0, 400))
    th = np.linspace(80, 150, 200); t0, kt = 109.5, 0.1
    axes[1].plot(th, 0.5 * kt * (th - t0) ** 2, color=C[0])
    axes[1].set(title='Ângulo', xlabel='θ (graus)')
    ph = np.linspace(-180, 180, 400)
    axes[2].plot(ph, 5 * (1 + np.cos(np.radians(3 * ph))), color=C[0])
    axes[2].set(title='Diedro', xlabel='φ (graus)', xticks=[-180, -60, 60, 180])
    rr = np.linspace(0.28, 1.2, 400); sg, e = 0.315, 0.636
    axes[3].plot(rr, 4 * e * ((sg / rr) ** 12 - (sg / rr) ** 6), color=C[0])
    axes[3].axhline(0, color=INK, lw=0.6)
    axes[3].set(title='Lennard-Jones (O-O, TIP3P)', xlabel='r (nm)', ylim=(-1.2, 2.5))
    rq = np.linspace(0.1, 1.2, 400)
    axes[4].plot(rq, 138.935 * (-0.834 * 0.417) / rq, color=C[0], label='O⁻ e H⁺ (atrai)')
    axes[4].plot(rq, 138.935 * (0.417 * 0.417) / rq, color=C[0], ls=(0, (6, 3)), label='H⁺ e H⁺ (repele)')
    axes[4].axhline(0, color=INK, lw=0.6)
    axes[4].set(title='Coulomb (cargas TIP3P, vácuo)', xlabel='r (nm)', ylim=(-500, 250))
    axes[4].legend(loc='lower right', fontsize=8)
    fig.suptitle('Termos de um campo de força clássico', y=1.05, fontsize=11)
    fig.tight_layout()
    save(fig, 'fund_campo_de_forca.png')


# ---------------------------------------------------------------- 3. passo de tempo
def fig_timestep():
    period = 11.0  # fs, estiramento C-H (~3000 cm-1)
    t = np.linspace(0, 60, 2000)
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.4), sharey=True)
    for ax, dt in zip(axes, (1, 2, 10)):
        ts = np.arange(0, 60 + 1e-9, dt)
        l1, = ax.plot(t, np.cos(2 * np.pi * t / period), color=INK, lw=0.8)
        l2, = ax.plot(ts, np.cos(2 * np.pi * ts / period), 'o-', color=C[0], ms=4, lw=1.2)
        ax.set(title=f'Δt = {dt} fs = {dt}×10⁻⁶ ns', xlabel='Tempo (fs; 1 fs = 10⁻⁶ ns)', ylim=(-1.15, 1.15))
    axes[0].set_ylabel('Deslocamento (rel.)')
    fig.legend([l1, l2], ['vibração C-H real (período ~11 fs)', 'o que o integrador enxerga'],
               loc='upper center', ncol=2, bbox_to_anchor=(0.5, 1.02), fontsize=9)
    fig.suptitle('O passo de tempo precisa ser muito menor que a vibração mais rápida', y=1.1, fontsize=11)
    fig.tight_layout()
    save(fig, 'fund_passo_de_tempo.png')


# ---------------------------------------------------------------- 4. PBC e cutoff
def fig_pbc():
    rng = np.random.default_rng(3)
    pts = rng.uniform(0.05, 0.95, (22, 2))
    fig, ax = plt.subplots(figsize=(6.4, 6.9))
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            central = i == 0 and j == 0
            ax.add_patch(Rectangle((i, j), 1, 1, fill=False, ec=INK, lw=1.6 if central else 0.6,
                                   ls='-' if central else (0, (4, 3))))
            ax.scatter(pts[:, 0] + i, pts[:, 1] + j, s=28, color=C[0], alpha=1 if central else 0.3, lw=0)
    p = pts[0]
    ax.add_patch(Circle(p, 0.38, fill=False, ec=C[1], lw=1.6))
    ax.scatter(*p, s=60, color=C[1], zorder=3)
    ax.annotate('esfera de corte $r_c$: a partícula interage\ncom vizinhos da caixa ao lado (imagens)',
                xy=(p[0] - 0.36, p[1] - 0.12), xytext=(-0.95, -0.85),
                arrowprops=dict(arrowstyle='-', color=C[1]), color=C[1], fontsize=9,
                bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='none'))
    ax.text(0.5, 1.03, 'caixa central', ha='center', va='bottom', color=INK)
    ax.set_xlim(-1, 2); ax.set_ylim(-1, 2.1); ax.set_aspect('equal'); frame_off(ax)
    ax.set_title('Condições periódicas de contorno')
    save(fig, 'fund_pbc_cutoff.png')


# ---------------------------------------------------------------- 5. mesma molecula em QM, AA, UA e CG
def read_res(path, resid, fmt):
    rows = []
    for l in open(path):
        if fmt == 'gro':
            if len(l) < 44 or not l[:5].strip().isdigit():
                continue
            if int(l[:5]) != resid:
                continue
            name = l[10:15].strip(); xyz = [float(l[20:28]), float(l[28:36]), float(l[36:44])]
            rows.append((name, np.array(xyz) * 10))
        else:
            if not l.startswith('ATOM') or int(l[22:26]) != resid:
                continue
            name = l[12:16].strip(); xyz = [float(l[30:38]), float(l[38:46]), float(l[46:54])]
            rows.append((name, np.array(xyz)))
    return rows


def project(rows, ref):
    X = np.array([x for _, x in rows]) - ref.mean(0)
    u, s, vt = np.linalg.svd(ref - ref.mean(0))
    return X @ vt[:2].T


def bonds(P, names, cut):
    out = []
    for i in range(len(P)):
        for j in range(i + 1, len(P)):
            if np.linalg.norm(P[i] - P[j]) < cut:
                out.append((i, j))
    return out


def fig_resolutions(aa, ua, cg, resid=45):
    A, U, G = read_res(aa, resid, 'gro'), read_res(ua, resid, 'gro'), read_res(cg, resid, 'pdb')
    ref = np.array([x for n, x in A if not n.startswith('H')])
    PA, PU, PG = project(A, ref), project(U, ref), project(G, ref)
    is_h = lambda n: n.startswith('H')
    fig, axes = plt.subplots(1, 4, figsize=(15, 4.4))
    # QM: densidade eletronica esquematica (soma de gaussianas)
    gx, gy = np.meshgrid(np.linspace(-6, 6, 300), np.linspace(-6, 6, 300))
    dens = np.zeros_like(gx)
    for (n, _), p in zip(A, PA):
        w, s = (0.4, 0.45) if is_h(n) else (1.0, 0.6)
        dens += w * np.exp(-((gx - p[0]) ** 2 + (gy - p[1]) ** 2) / (2 * s ** 2))
    axes[0].contourf(gx, gy, dens, levels=12, cmap='Blues')
    axes[0].scatter(PA[:, 0], PA[:, 1], s=6, color=INK)
    titles = [f'QM: núcleos + elétrons', f'All-atom: {len(A)} átomos', f'United-atom: {len(U)} partículas',
              f'Martini 3: {len(G)} beads']
    # AA e UA: atomos e ligacoes
    for ax, (rows, P) in zip(axes[1:3], ((A, PA), (U, PU))):
        nm = [n for n, _ in rows]
        cut = [1.25 if (is_h(nm[i]) or is_h(nm[j])) else 1.75 for i in range(len(P)) for j in range(len(P))]
        for i, j in bonds(P, nm, 1.75):
            lim = 1.25 if (is_h(nm[i]) or is_h(nm[j])) else 1.75
            if np.linalg.norm(np.array([x for _, x in rows][i]) - np.array([x for _, x in rows][j])) < lim:
                ax.plot(*P[[i, j]].T, color=INK, lw=1.2, zorder=1)
        for (n, _), p in zip(rows, P):
            ax.scatter(*p, s=40 if is_h(n) else 110, color=C[0] if is_h(n) else INK, zorder=2,
                       edgecolor='white', lw=0.8)
    # CG: beads com raio aproximado (R ~2,35 A, S ~2,05 A, T ~1,7 A)
    for (n, _), p in zip(G, PG):
        r = 2.35 if n == 'BB' else 2.05
        axes[3].add_patch(Circle(p, r, color=C[1] if n == 'BB' else C[2], alpha=0.75, lw=0))
        axes[3].text(*p, n, ha='center', va='center', fontsize=8, color='white')
    axes[3].scatter(PA[[i for i, (n, _) in enumerate(A) if not is_h(n)], 0],
                    PA[[i for i, (n, _) in enumerate(A) if not is_h(n)], 1], s=8, color=INK, zorder=3)
    for ax, t in zip(axes, titles):
        ax.set_xlim(-6, 6); ax.set_ylim(-6, 6); ax.set_aspect('equal'); frame_off(ax)
        ax.set_title(t)
    axes[1].scatter([], [], s=40, color=C[0], label='hidrogênio'); axes[1].scatter([], [], s=110, color=INK, label='átomo pesado')
    axes[1].legend(loc='lower left', fontsize=8)
    fig.suptitle('Fenilalanina 45 da ubiquitina em quatro resoluções (coordenadas reais; densidade QM esquemática)',
                 y=0.93, fontsize=11)
    save(fig, 'fund_resolucoes_molecula.png')


# ---------------------------------------------------------------- 6. escalas de tamanho e tempo
def fig_scales():
    # ordens de grandeza tipicas: (nome, tamanho nm min-max, tempo ns min-max)
    regions = [('QM', (0.1, 3), (1e-7, 1e-2)), ('All-atom', (1, 30), (1e-6, 1e4)),
               ('United-atom', (1, 50), (1e-6, 1e4)), ('Coarse-grained', (3, 500), (1e-4, 1e6))]
    fig, axes = plt.subplots(1, 2, figsize=(13, 3.6), sharey=True)
    for k, (name, size, time) in enumerate(regions):
        axes[0].barh(k, size[1] - size[0], left=size[0], color=C[k], height=0.55)
        axes[1].barh(k, time[1] - time[0], left=time[0], color=C[k], height=0.55)
    axes[0].set(xscale='log', xlim=(0.05, 2000), xlabel='Tamanho do sistema (nm)', title='Tamanho acessível',
                yticks=range(len(regions)), yticklabels=[r[0] for r in regions])
    axes[0].invert_yaxis()
    events = [('vibração', 1e-5), ('cadeia lateral', 1e-1), ('alças', 1e1), ('bicamada', 1e2),
              ('enovelamento', 1e4), ('fusão', 1e6)]
    for lab, tt in events:
        axes[1].axvline(tt, color=INK, lw=0.6, ls=(0, (3, 3)))
        axes[1].text(tt, -0.75, lab, rotation=45, ha='left', va='bottom', fontsize=8, color=INK)
    axes[1].set(xscale='log', xlim=(1e-7, 1e8), xlabel='Tempo acessível (ns)')
    axes[1].set_title('Tempo acessível e processos', pad=70)
    fig.suptitle('Alcance típico de cada resolução (ordens de grandeza)', y=1.12, fontsize=11)
    fig.tight_layout()
    save(fig, 'fund_escalas.png')


# ---------------------------------------------------------------- 7. mapeamento e LJ
def fig_mapping_lj():
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.4))
    ax = axes[0]; rng = np.random.default_rng(7)
    centers = np.array([[0, 0], [0.28, 0.05], [0.1, 0.27], [-0.18, 0.2]])
    for c0 in centers:
        ang = rng.uniform(0, 2 * np.pi)
        h1 = c0 + 0.096 * np.array([np.cos(ang), np.sin(ang)])
        h2 = c0 + 0.096 * np.array([np.cos(ang + 1.82), np.sin(ang + 1.82)])
        for h in (h1, h2):
            ax.plot(*np.array([c0, h]).T, color=INK, lw=1)
            ax.add_patch(Circle(h, 0.03, color=C[0]))
        ax.add_patch(Circle(c0, 0.055, color=INK))
    com = centers.mean(0)
    ax.add_patch(Circle(com, 0.235, fc=C[6], alpha=0.22, ec=C[6], lw=1.5))
    ax.text(com[0], com[1] - 0.29, '1 bead W (σ = 0,47 nm)', ha='center', color=C[6])
    ax.set_xlim(-0.5, 0.6); ax.set_ylim(-0.4, 0.5); ax.set_aspect('equal'); frame_off(ax)
    ax.set_title('Mapeamento 4:1: quatro H₂O viram um bead W')
    ax = axes[1]; r = np.linspace(0.25, 1.2, 500)
    for (lab, s, e, col, ls) in (('O-O, TIP3P (σ 0,315 nm; ε 0,64 kJ/mol)', 0.315, 0.636, C[0], '-'),
                                 ('W-W, Martini 3 (σ 0,47 nm; ε 4,65 kJ/mol)', 0.47, 4.65, C[6], (0, (6, 3)))):
        ax.plot(r, 4 * e * ((s / r) ** 12 - (s / r) ** 6), color=col, ls=ls, label=lab)
    ax.axhline(0, color=INK, lw=0.6)
    ax.set(ylim=(-5.5, 5), xlabel='Distância (nm)', ylabel='U (kJ/mol)', title='Lennard-Jones: átomo × bead')
    ax.legend(loc='lower right', fontsize=8)
    save(fig, 'fund_mapeamento_lj.png')


# ---------------------------------------------------------------- 8. paisagem 2D
def landscape(x, y, rough):
    base = 3 * (x ** 2 - 1) ** 2 + 2.5 * y ** 2
    return base + rough * (0.8 * np.sin(9 * x) * np.cos(8 * y) + 0.5 * np.sin(17 * x + 1) * np.sin(15 * y))


def fig_landscape_2d():
    g = np.linspace(-1.7, 1.7, 300)
    X, Y = np.meshgrid(g, g)
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.3))
    for ax, (name, rough) in zip(axes, (('All-atom', 1.0), ('United-atom', 0.55), ('Coarse-grained', 0.0))):
        Z = landscape(X, Y, rough) * (0.6 if rough == 0 else 1)
        cs = ax.contourf(X, Y, Z - Z.min(), levels=18, cmap='pure_r')
        ax.contour(X, Y, Z - Z.min(), levels=18, colors='white', linewidths=0.3)
        ax.set(title=name, xlabel='Coordenada 1', ylabel='Coordenada 2' if ax is axes[0] else None)
        ax.set_aspect('equal')
    fig.suptitle('Paisagem de energia livre (esquema ilustrativo): o CG integra a rugosidade e mantém os vales',
                 y=1.03, fontsize=11)
    save(fig, 'fund_paisagem_2d.png')


# ---------------------------------------------------------------- 9. relogios diferentes (modelo de brinquedo)
def fig_clocks():
    """Dinamica de Langevin superamortecida em 1D: mesma barreira, rugosidade diferente."""
    rng = np.random.default_rng(11)
    kT, dt, n = 1.0, 2e-4, 400000

    def U(x, rough):
        return 3 * (x ** 2 - 1) ** 2 * (0.6 if rough == 0 else 1) + rough * (0.9 * np.sin(23 * x) + 0.6 * np.sin(41 * x + 1.3))

    def dU(x, rough, h=1e-5):
        return (U(x + h, rough) - U(x - h, rough)) / (2 * h)

    def run(rough):
        x = np.empty(n); x[0] = -1.0
        noise = rng.standard_normal(n) * np.sqrt(2 * kT * dt)
        for i in range(1, n):
            x[i] = x[i - 1] - dU(x[i - 1], rough) * dt + noise[i]
        return x

    def crossings(x):
        side = np.sign(x); side[side == 0] = 1
        # conta transicoes entre os vales com histerese em +-0.6
        state, c = -1, 0
        for v in x:
            if state < 0 and v > 0.6: state, c = 1, c + 1
            elif state > 0 and v < -0.6: state, c = -1, c + 1
        return c

    xs = np.linspace(-1.7, 1.7, 600)
    traj = {'All-atom (rugosa)': (1.0, C[0]), 'Coarse-grained (lisa)': (0.0, C[1])}
    fig = plt.figure(figsize=(13, 6.2))
    gs = fig.add_gridspec(2, 2, width_ratios=[1, 2.4], hspace=0.45, wspace=0.18)
    for k, (name, (rough, col)) in enumerate(traj.items()):
        axl = fig.add_subplot(gs[k, 0]); axr = fig.add_subplot(gs[k, 1])
        u = U(xs, rough)
        axl.plot(xs, u - u.min(), color=col)
        axl.set(title=name, xlabel='x', ylabel='U / kT', ylim=(0, 8))
        x = run(rough)
        steps = np.arange(n)[::40]
        axr.plot(steps, x[::40], color=col, lw=0.5)
        axr.axhline(0, color=INK, lw=0.5, ls=(0, (4, 3)))
        axr.set(ylabel='x', xlabel='Passos de integração' if k == 1 else None, ylim=(-1.8, 1.8),
                title=f'{crossings(x)} travessias da barreira em {n:,} passos'.replace(',', ' '))
    fig.suptitle('Mesmo número de passos, relógios diferentes (modelo de brinquedo, dinâmica de Langevin)',
                 y=1.0, fontsize=11)
    save(fig, 'fund_relogios.png')


# ---------------------------------------------------------------- 10. mapeamento de varios residuos
def fig_residue_mapping(aa, cg, picks=(('GLY', 10), ('ALA', 28), ('LEU', 8), ('LYS', 6), ('ARG', 42), ('PHE', 45), ('TYR', 59))):
    fig, axes = plt.subplots(1, len(picks), figsize=(2.6 * len(picks), 3.6))
    for ax, (name, resid) in zip(axes, picks):
        A = [(n, x) for n, x in read_res(aa, resid, 'gro')]
        G = read_res(cg, resid, 'pdb')
        heavy = np.array([x for n, x in A if not n.startswith('H')])
        PA, PG = project(A, heavy), project(G, heavy)
        for (n, _), p in zip(G, PG):
            ax.add_patch(Circle(p, 2.35 if n == 'BB' else 2.0, color=C[1] if n == 'BB' else C[2], alpha=0.35, lw=0))
        for (n, _), p in zip(A, PA):
            ax.scatter(*p, s=12 if n.startswith('H') else 34, color=C[0] if n.startswith('H') else INK, zorder=3, lw=0)
        ax.set_xlim(-6.5, 6.5); ax.set_ylim(-5, 5); ax.set_aspect('equal'); frame_off(ax)
        nh = sum(1 for n, _ in A if not n.startswith('H'))
        ax.set_title(f'{name}{resid}\n{nh} pesados, {len(A)} átomos\n{len(G)} bead' + ('s' if len(G) > 1 else ''), fontsize=9)
    fig.suptitle('Mapeamento Martini 3 de resíduos da ubiquitina: amarelo BB, magenta cadeia lateral', y=0.98, fontsize=11)
    save(fig, 'fund_mapeamento_residuos.png')


# ---------------------------------------------------------------- 11. minimizacao x dinamica
def fig_em_vs_md():
    rng = np.random.default_rng(5)
    U = lambda x: 3 * (x ** 2 - 1) ** 2 + 0.9 * np.sin(23 * x) + 0.6 * np.sin(41 * x + 1.3)
    dU = lambda x, h=1e-5: (U(x + h) - U(x - h)) / (2 * h)
    xs = np.linspace(-1.7, 1.7, 800)
    x0 = -0.35
    em = [x0]
    for _ in range(400):
        em.append(em[-1] - 0.002 * dU(em[-1]))
    md = [x0]
    for _ in range(60000):
        md.append(md[-1] - dU(md[-1]) * 2e-4 + rng.standard_normal() * np.sqrt(2 * 2e-4))
    md = np.array(md)
    fig, axes = plt.subplots(1, 2, figsize=(12, 3.8))
    for ax in axes:
        ax.plot(xs, U(xs) - U(xs).min(), color=INK, lw=1.2)
        ax.set(xlabel='Coordenada', ylabel='U / kT', ylim=(-0.5, 9))
    axes[0].plot(em, U(np.array(em)) - U(xs).min(), 'o', color=C[0], ms=3)
    axes[0].annotate('para no mínimo\nlocal mais próximo', xy=(em[-1], U(em[-1]) - U(xs).min()), xytext=(0.4, 6),
                     arrowprops=dict(arrowstyle='->', color=C[0]), color=C[0])
    axes[0].set_title('Minimização (steepest descent): só desce')
    h, edges = np.histogram(md, bins=80, range=(-1.7, 1.7))
    axes[1].plot(md[::30], U(md[::30]) - U(xs).min(), '.', color=C[1], ms=2, alpha=0.4)
    axes[1].set_title('Dinâmica a temperatura finita: visita vales e cruza barreiras')
    fig.suptitle('Mesma paisagem, dois algoritmos (esquema ilustrativo)', y=1.04, fontsize=11)
    save(fig, 'fund_minimizacao_vs_md.png')


# ---------------------------------------------------------------- 12. custo teorico
def fig_cost(V=295.5, n_w=2084, n_cg=2304, n_aa=1231, n_ua=762, n_ions=54):
    rows = [('All-atom', n_aa + 12 * n_w + n_ions, 1.2, 2), ('United-atom', n_ua + 12 * n_w + n_ions, 1.4, 2),
            ('Martini 3', n_cg, 1.1, 20)]
    names = [r[0] for r in rows]
    N = np.array([r[1] for r in rows]); rc = np.array([r[2] for r in rows]); dt = np.array([r[3] for r in rows])
    nviz = N / V * 4 / 3 * np.pi * rc ** 3
    steps = 1e6 / dt
    cost = N * nviz * steps
    fig, axes = plt.subplots(1, 4, figsize=(15, 3.4), sharey=True)
    for ax, vals, title in zip(axes, (N, nviz, steps, cost / cost[-1]),
                               ('Partículas', 'Vizinhos no cutoff', 'Passos por ns', 'Custo por ns ÷ Martini')):
        ax.barh(names, vals, color=C[0], height=0.55)
        ax.set_xscale('log'); ax.set_title(title)
        for i, v in enumerate(vals):
            ax.text(v * 1.15, i, f'{v:,.0f}'.replace(',', ' '), va='center', fontsize=8, color=INK)
        ax.set_xlim(vals.min() / 3, vals.max() * 8)
    axes[0].invert_yaxis()
    fig.suptitle('Custo teórico da ubiquitina solvatada: os três fatores se multiplicam', y=1.05, fontsize=11)
    fig.tight_layout()
    save(fig, 'fund_custo_teorico.png')


# ---------------------------------------------------------------- 13. automontagem em quadros
def fig_self_assembly(xtc, gro, rep_dir=None, times=(0, 2, 5, 8, 15)):
    """Quadros reais da automontagem: normal final na vertical, bicamada centrada em cada quadro."""
    import mdtraj as md
    from make_replica_figures import draw, ORDER
    from replica_views import wrap_molecules, center_along, itp_bonds, system_bonds
    t = md.load(xtc, top=gro)
    top = t.topology
    rep = Path(rep_dir) if rep_dir else Path(gro).parent
    bonds = system_bonds(top, {'DSPC': itp_bonds(rep / 'martini_v3.0.0_phospholipids_PC_v2.itp', 'DSPC')})
    lip = top.select('resname DSPC'); wat = top.select('resname W')
    tails = top.select('name =~ "C[1-4][AB]"')
    names = [a.name for a in top.atoms]
    o = np.load(rep / 'mem_order.npz')
    axis = int(np.argmax(np.abs(o['normal'])))
    perm = list(ORDER[axis])                          # normal final -> z
    idx = {int(a): k for k, a in enumerate(lip)}
    lb = [(idx[i], idx[j]) for i, j in bonds if i in idx and j in idx]
    lipnames = [names[i] for i in lip]
    L = max(t.unitcell_lengths[0])
    fig, axes = plt.subplots(1, len(times), figsize=(3.2 * len(times), 4.0))
    for ax, tn in zip(axes, times):
        i = int(np.argmin(np.abs(t.time / 1000 - tn)))
        box = t.unitcell_lengths[i][perm]
        xyz = t.xyz[i][:, perm]
        xyz = center_along(xyz, box, tails, 2)
        xyz = wrap_molecules(xyz, box, top)
        draw(ax, xyz[lip], lipnames, lb, 0, 2, 1, box_hv=(box[0], box[2]), water=xyz[wat],
             extent=((-0.5, L + 0.5), (-0.5, L + 0.5)))
        ax.set_title(f'{t.time[i] / 1000:.0f} ns')
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], ls='none', marker='o', ms=9, mfc=SEM_[k], mec=INK, mew=0.4)
               for k in ('head_cat', 'po4', 'glyc', 'tail', 'water')]
    fig.legend(handles, ['Colina NC3', 'Fosfato PO4', 'Glicerol', 'Caudas', 'Água W'], loc='upper center',
               ncol=5, frameon=False, bbox_to_anchor=(0.5, 0.06))
    fig.suptitle('Automontagem de 128 DSPC, réplica 1 (simulação real, vista lateral, normal final na vertical)',
                 y=0.99, fontsize=11)
    fig.tight_layout(rect=(0, 0.07, 1, 0.95))
    fig.savefig(OUT / 'fund_automontagem_quadros.png', dpi=DPI, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print('   fund_automontagem_quadros.png')

# ---------------------------------------------------------------- 14. funil de energia
def fig_funnel():
    """Funil de enovelamento em 4 resolucoes (esquema ilustrativo, superficies analiticas)."""
    r = np.linspace(0, 1, 120); th = np.linspace(0, 2 * np.pi, 160)
    R, T = np.meshgrid(r, th)
    X, Y = R * np.cos(T), R * np.sin(T)
    base = 6 * R ** 1.4                                   # funil: energia cai rumo ao nativo (centro)
    rough = lambda a, k: a * (np.sin(k * R * 6) * np.cos(5 * T) + 0.6 * np.sin(k * 11 * R + 3 * T))
    panels = [('QM', base + rough(1.0, 1.0), 'mesma superfície, mais exata;\nacessível só na região em destaque'),
              ('All-atom', base + rough(1.0, 1.0), 'rugosa: muitos mínimos locais'),
              ('United-atom', base + rough(0.55, 0.9), 'rugosidade menor'),
              ('Martini 3 + rede elástica', 40 * R ** 2, 'poço único em torno da estrutura\ninicial: sem enovelamento')]
    fig = plt.figure(figsize=(16, 4.6))
    for k, (name, Z, note) in enumerate(panels):
        ax = fig.add_subplot(1, 4, k + 1, projection='3d')
        Zp = np.where(Z <= 8, Z, np.nan)
        fc = plt.get_cmap('pure_r')(np.clip(np.nan_to_num(Zp, nan=8) / 8, 0, 1))
        if name == 'QM':
            fc[R < 0.3] = matplotlib.colors.to_rgba(C[1])   # regiao acessivel ao QM
        ax.plot_surface(X, Y, Zp, facecolors=fc, linewidth=0, antialiased=True, rstride=2, cstride=2, shade=True)
        ax.set_zlim(0, 8); ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
        ax.set_zlabel('Energia livre', labelpad=-8, fontsize=9)
        ax.view_init(elev=24, azim=-60)
        ax.set_title(name, pad=0)
        ax.text2D(0.5, -0.04, note, transform=ax.transAxes, ha='center', va='top', fontsize=8.5, color=INK)
    fig.suptitle('Funil de energia do enovelamento em quatro resoluções (esquema ilustrativo). '
                 'Largura: número de conformações; profundidade: estabilidade', y=0.92, fontsize=11)
    save(fig, 'fund_funil.png')


# ---------------------------------------------------------------- 15. membranas de celulas reais
LIPID = {  # nome: (cor, carga da cabeca, rotulo)
    'PC':   ('#4c72b0', 0, 'PC (fosfatidilcolina)'),
    'SM':   ('#64b5cd', 0, 'SM (esfingomielina)'),
    'PE':   ('#8172b2', 0, 'PE (fosfatidiletanolamina)'),
    'PS':   ('#c44e52', -1, 'PS (fosfatidilserina)'),
    'PG':   ('#dd8452', -1, 'PG (fosfatidilglicerol)'),
    'CL':   ('#937860', -2, 'CL (cardiolipina)'),
    'CHOL': ('#ccb974', 0, 'Colesterol'),
}


def _leaflet(ax, comp, y, up, n=24, seed=0):
    """Desenha uma monocamada: n lipidios sorteados pela composicao (fracoes)."""
    rng = np.random.default_rng(seed)
    names = list(comp)
    p = np.array([comp[k] for k in names], float); p /= p.sum()
    counts = np.round(p * n).astype(int)
    seq = [k for k, c in zip(names, counts) for _ in range(c)]
    rng.shuffle(seq)
    for i, k in enumerate(seq[:n]):
        x = i + 0.5
        color, q, _ = LIPID[k]
        d = -1 if up else 1
        if k == 'CHOL':
            ax.add_patch(Rectangle((x - 0.18, y + d * 0.3 - (0.9 if d < 0 else 0)), 0.36, 0.9, color=color, lw=0))
            ax.add_patch(Circle((x, y), 0.16, color=color, ec=INK, lw=0.4))
            continue
        for dx in (-0.12, 0.12):
            ax.plot([x + dx, x + dx], [y + d * 0.35, y + d * 1.55], color='#9a9a9a', lw=1.4, solid_capstyle='round')
        ax.add_patch(Circle((x, y), 0.36, color=color, ec=INK, lw=0.5))
        if q:
            ax.text(x, y, '−' * abs(q), ha='center', va='center', color='white', fontsize=8, fontweight='bold')


def fig_cell_membranes():
    panels = [
        ('Célula humana normal', 'superfície externa neutra',
         {'PC': 45, 'SM': 25, 'CHOL': 30}, {'PE': 40, 'PS': 15, 'PC': 15, 'CHOL': 30}),
        ('Célula de câncer', 'PS exposta: superfície externa aniônica',
         {'PC': 40, 'SM': 20, 'PS': 10, 'CHOL': 30}, {'PE': 40, 'PS': 10, 'PC': 20, 'CHOL': 30}),
        ('Bactéria (membrana interna de E. coli)', 'PG e CL: superfície aniônica, sem colesterol',
         {'PE': 75, 'PG': 20, 'CL': 5}, {'PE': 75, 'PG': 20, 'CL': 5}),
    ]
    fig, axes = plt.subplots(3, 1, figsize=(12, 9.6))
    for k, (ax, (title, sub, outer, inner)) in enumerate(zip(axes, panels)):
        ax.axhspan(2.2, 3.6, color='#a6cfe3', alpha=0.35, lw=0)
        ax.axhspan(-3.6, -2.2, color='#a6cfe3', alpha=0.35, lw=0)
        ax.text(24.4, 2.9, 'periplasma' if k == 2 else 'meio externo', ha='left', va='center', fontsize=8, color=INK)
        ax.text(24.4, -2.9, 'citoplasma', ha='left', va='center', fontsize=8, color=INK)
        _leaflet(ax, outer, 1.9, up=True, seed=k)
        _leaflet(ax, inner, -1.9, up=False, seed=k + 10)
        ax.set_xlim(0, 27); ax.set_ylim(-3.6, 3.6); ax.set_aspect('equal'); frame_off(ax)
        ax.set_title(f'{title}: {sub}', fontsize=10.5)
    used = ['PC', 'SM', 'PE', 'PS', 'PG', 'CL', 'CHOL']
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], ls='none', marker='o', ms=10, mfc=LIPID[n][0], mec=INK, mew=0.5) for n in used]
    fig.legend(handles, [LIPID[n][2] + (' −' * abs(LIPID[n][1]) if LIPID[n][1] else '') for n in used],
               loc='lower center', ncol=4, frameon=False, fontsize=8.5, bbox_to_anchor=(0.5, -0.01))
    fig.suptitle('Composição das membranas por monocamada (esquema; proporções aproximadas)', y=0.995, fontsize=11)
    fig.tight_layout(rect=(0, 0.06, 1, 0.98))
    fig.savefig(OUT / 'fund_membranas_celulas.png', dpi=DPI, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print('   fund_membranas_celulas.png')


def main():
    ap = argparse.ArgumentParser(description='Gera imagens/fund_*.png')
    ap.add_argument('--aa', required=True); ap.add_argument('--ua', required=True); ap.add_argument('--cg', required=True)
    ap.add_argument('--mem-xtc'); ap.add_argument('--mem-gro')
    a = ap.parse_args()
    OUT.mkdir(exist_ok=True)
    fig_md_loop(); fig_force_field(); fig_timestep(); fig_pbc()
    fig_resolutions(a.aa, a.ua, a.cg); fig_scales(); fig_mapping_lj(); fig_landscape_2d(); fig_clocks()
    fig_residue_mapping(a.aa, a.cg); fig_em_vs_md(); fig_cost(); fig_funnel(); fig_cell_membranes()
    if a.mem_xtc:
        fig_self_assembly(a.mem_xtc, a.mem_gro)


if __name__ == '__main__':
    main()
