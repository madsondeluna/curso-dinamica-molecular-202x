"""Figuras estaticas das replicas para o README (imagens/martini_vistas_*.png).

Projecao ortografica em matplotlib: mesma escala (nm) em todos os paineis,
caixa periodica desenhada, esferas com o raio do bead ordenadas por
profundidade e bastoes nas ligacoes CG. Cores do vmd_view.sh do AMP-evo-MD.

Uso: python tools/make_replica_figures.py runs/aula_martini
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np
from matplotlib.collections import LineCollection, PatchCollection
from matplotlib.colors import to_rgb
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Rectangle

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'viz'))
from replica_plots import STYLE, INK  # noqa: E402
from replica_views import SEM, itp_bonds, system_bonds, wrap_molecules, center_along  # noqa: E402

OUT = ROOT / 'imagens'
plt.rcParams.update(STYLE)

# raio de desenho (nm) e cor de cada bead
BEAD = {'NC3': (0.17, SEM['head_cat']), 'PO4': (0.16, SEM['po4']),
        'GL1': (0.14, SEM['glyc']), 'GL2': (0.14, SEM['glyc'])}
TAIL = (0.09, SEM['tail'])
ORDER = {0: (1, 2, 0), 1: (2, 0, 1), 2: (0, 1, 2)}   # permutacao ciclica que leva o eixo da normal a z


def shade(color, depth01):
    """Clareia o que esta ao fundo (depth01 = 0 frente, 1 fundo)."""
    c = np.array(to_rgb(color))
    return tuple(c + (1 - c) * 0.55 * depth01)


def draw(ax, xyz, names, bonds, h, v, d, box_hv=None, water=None, extent=None):
    """Desenha uma projecao: eixos de tela h (horizontal) e v (vertical), profundidade d."""
    if water is not None and len(water):
        ax.scatter(water[:, h], water[:, v], s=3, color=SEM['water'], alpha=0.35, lw=0, zorder=0)
    depth = xyz[:, d]
    dn = (depth - depth.min()) / max(np.ptp(depth), 1e-6)    # 0 fundo .. 1 frente
    if bonds:
        segs = [[(xyz[i, h], xyz[i, v]), (xyz[j, h], xyz[j, v])] for i, j in bonds]
        bd = np.array([(dn[i] + dn[j]) / 2 for i, j in bonds])
        cols = [shade(SEM['tail'], 1 - x) for x in bd]
        ax.add_collection(LineCollection(segs, colors=cols, linewidths=1.3, zorder=1, capstyle='round'))
    order = np.argsort(dn)                                    # fundo primeiro
    patches, faces = [], []
    for k in order:
        r, c = BEAD.get(names[k], TAIL)
        patches.append(Circle((xyz[k, h], xyz[k, v]), r))
        faces.append(shade(c, 1 - dn[k]))
    ax.add_collection(PatchCollection(patches, facecolors=faces, edgecolors=[shade(INK, 0.6)] * len(patches),
                                      linewidths=0.25, zorder=2, match_original=False))
    if box_hv is not None:
        ax.add_patch(Rectangle((0, 0), box_hv[0], box_hv[1], fill=False, ec=INK, lw=0.9, zorder=3))
    if extent is not None:
        ax.set_xlim(extent[0]); ax.set_ylim(extent[1])
    ax.set_aspect('equal'); ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)


def permute(xyz, axis):
    return xyz[:, list(ORDER[axis])]


def membrane_figure(base):
    reps = [base / f'rep{n}' for n in (1, 2, 3)]
    panels = []
    for rd in reps:
        top_t = md.load(str(rd / 'mem_em.gro')).topology
        lip_bonds = itp_bonds(rd / 'martini_v3.0.0_phospholipids_PC_v2.itp', 'DSPC')
        bonds = system_bonds(top_t, {'DSPC': lip_bonds})
        lip = top_t.select('resname DSPC'); wat = top_t.select('resname W')
        tails = top_t.select('name =~ "C[1-4][AB]"')
        names = [a.name for a in top_t.atoms]
        traj = md.load(str(rd / 'mem_whole.xtc'), top=str(rd / 'mem_em.gro'))
        o = np.load(rd / 'mem_order.npz')
        formed = (rd / 'mem_eq_whole.xtc').exists()
        axis = int(np.argmax(np.abs(o['normal']))) if formed else 2
        # inicio: primeiro quadro, mesma rotacao do fim
        x0 = wrap_molecules(traj.xyz[0], traj.unitcell_lengths[0], top_t)
        b0 = traj.unitcell_lengths[0]
        x0, b0 = permute(x0, axis), b0[list(ORDER[axis])]
        if formed:
            te = md.load(str(rd / 'mem_eq_whole.xtc'), top=str(rd / 'mem_rot.gro'))
            x1, b1 = te.xyz[-1], te.unitcell_lengths[-1]
            x1 = center_along(x1, b1, tails, 2)
        else:
            x1, b1 = traj.xyz[-1], traj.unitcell_lengths[-1]
            x1, b1 = permute(x1, axis), b1[list(ORDER[axis])]
        x1 = wrap_molecules(x1, b1, top_t)
        panels.append((x0, b0, x1, b1, names, bonds, lip, wat, formed))

    big = max(max(p[1].max(), p[3].max()) for p in panels)
    ext = ((-0.6, big + 0.6), (-0.6, big + 0.6))
    fig, axes = plt.subplots(3, 3, figsize=(12, 12.6))
    titles = ['Início: lateral (x-z)', 'Fim: lateral (x-z)', 'Fim: topo (x-y)']
    for r, (x0, b0, x1, b1, names, bonds, lip, wat, formed) in enumerate(panels):
        lipnames = [names[i] for i in lip]
        idx = {int(o): k for k, o in enumerate(lip)}
        lb = [(idx[i], idx[j]) for i, j in bonds if i in idx and j in idx]
        for c, (xyz, box, h, v, d) in enumerate(((x0, b0, 0, 2, 1), (x1, b1, 0, 2, 1), (x1, b1, 0, 1, 2))):
            draw(axes[r, c], xyz[lip], lipnames, lb, h, v, d, box_hv=(box[h], box[v]),
                 water=xyz[wat] if c < 2 else None, extent=ext)
        state = 'bicamada' if formed else 'sem bicamada (fase contínua)'
        axes[r, 0].text(-0.08, 0.5, f'Réplica {r + 1}\n{state}', transform=axes[r, 0].transAxes,
                        rotation=90, ha='center', va='center', fontsize=10, color=INK)
    for c, t in enumerate(titles):
        axes[0, c].set_title(t)
    handles = [Line2D([], [], ls='none', marker='o', ms=10, mfc=SEM[k], mec=INK, mew=0.4) for k in ('head_cat', 'po4', 'glyc', 'tail', 'water')]
    fig.legend(handles, ['Colina NC3', 'Fosfato PO4', 'Glicerol GL1, GL2', 'Caudas', 'Água W'],
               loc='upper center', ncol=5, frameon=False, bbox_to_anchor=(0.5, 0.025))
    fig.suptitle('Automontagem de 128 DSPC em três réplicas (Martini 3, projeção ortográfica, caixa periódica em preto)',
                 y=0.995, fontsize=11)
    fig.tight_layout(rect=(0.02, 0.03, 1, 0.98))
    fig.savefig(OUT / 'martini_vistas_membrana.png', dpi=200, facecolor='white', bbox_inches='tight', pad_inches=0.2)
    plt.close(fig)
    print('  martini_vistas_membrana.png')


def protein_figure(base):
    reps = [base / f'rep{n}' for n in (1, 2, 3)]
    ref = md.load(str(reps[0] / 'ubq_cg.pdb'))
    names = [a.name for a in ref.topology.atoms]
    bonds = itp_bonds(reps[0] / 'molecule_0.itp')
    # orientacao comum: eixos principais da estrutura inicial
    x0 = ref.xyz[0] - ref.xyz[0].mean(0)
    _, _, vt = np.linalg.svd(x0, full_matrices=False)
    if np.linalg.det(vt) < 0:
        vt[2] *= -1
    rot = lambda x: (x - x.mean(0)) @ vt.T
    bb = ref.topology.select('name BB')
    rows = []
    for rd in reps:
        t = md.load(str(rd / 'md_prot_fit.xtc'), top=str(rd / 'ubq_cg.pdb'))
        t.superpose(ref, atom_indices=bb)            # mesma orientacao para inicio, fim e replicas
        rows.append((rot(t.xyz[0]), rot(t.xyz[-1])))
    lim = max(np.abs(np.concatenate([np.concatenate(r) for r in rows])).max(), 1) + 0.4
    ext = ((-lim, lim), (-lim, lim))
    pc = {n: (0.23, SEM['protein']) for n in ['BB']}
    fig, axes = plt.subplots(3, 3, figsize=(11, 11.2))
    for r, (a, b) in enumerate(rows):
        for c, (xyz, h, v, d) in enumerate(((a, 0, 1, 2), (b, 0, 1, 2), (b, 2, 1, 0))):
            ax = axes[r, c]
            draw_protein(ax, xyz, names, bonds, h, v, d, ext)
        axes[r, 0].text(-0.08, 0.5, f'Réplica {r + 1}', transform=axes[r, 0].transAxes, rotation=90,
                        ha='center', va='center', fontsize=10, color=INK)
    for c, t in enumerate(['Início (0 ns)', 'Fim (15 ns): frente', 'Fim (15 ns): lateral, 90°']):
        axes[0, c].set_title(t)
    handles = [Line2D([], [], ls='none', marker='o', ms=10, mfc=SEM[k], mec=INK, mew=0.4) for k in ('protein', 'protein_sc')]
    fig.legend(handles, ['Backbone (BB)', 'Cadeia lateral (SC)'], loc='upper center', ncol=2, frameon=False,
               bbox_to_anchor=(0.5, 0.025))
    fig.suptitle('Ubiquitina em três réplicas (Martini 3, projeção ortográfica nos eixos principais da estrutura inicial)',
                 y=0.995, fontsize=11)
    fig.tight_layout(rect=(0.02, 0.03, 1, 0.98))
    fig.savefig(OUT / 'martini_vistas_proteina.png', dpi=200, facecolor='white', bbox_inches='tight', pad_inches=0.2)
    plt.close(fig)
    print('  martini_vistas_proteina.png')


def draw_protein(ax, xyz, names, bonds, h, v, d, ext):
    depth = xyz[:, d]
    dn = (depth - depth.min()) / max(np.ptp(depth), 1e-6)
    bbset = {k for k, n in enumerate(names) if n == 'BB'}
    for sel, lw, col in (([b for b in bonds if not (b[0] in bbset and b[1] in bbset)], 1.4, SEM['protein_sc']),
                         ([b for b in bonds if b[0] in bbset and b[1] in bbset], 3.2, SEM['protein'])):
        segs = [[(xyz[i, h], xyz[i, v]), (xyz[j, h], xyz[j, v])] for i, j in sel]
        cols = [shade(col, 1 - (dn[i] + dn[j]) / 2) for i, j in sel]
        ax.add_collection(LineCollection(segs, colors=cols, linewidths=lw, zorder=1, capstyle='round'))
    patches, faces = [], []
    for k in np.argsort(dn):
        bb = k in bbset
        patches.append(Circle((xyz[k, h], xyz[k, v]), 0.1 if bb else 0.075))
        faces.append(shade(SEM['protein'] if bb else SEM['protein_sc'], 1 - dn[k]))
    ax.add_collection(PatchCollection(patches, facecolors=faces, edgecolors=[shade(INK, 0.6)] * len(patches),
                                      linewidths=0.3, zorder=2, match_original=False))
    ax.set_xlim(ext[0]); ax.set_ylim(ext[1]); ax.set_aspect('equal')
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)


if __name__ == '__main__':
    base = Path(sys.argv[1] if len(sys.argv) > 1 else ROOT / 'runs' / 'aula_martini')
    membrane_figure(base)
    protein_figure(base)
