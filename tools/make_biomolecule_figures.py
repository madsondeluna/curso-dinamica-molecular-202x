"""Gera as figuras da introducao sobre biomoleculas (imagens/bio_*.png).

Coordenadas reais: estruturas do PDB e coordenadas ideais do Chemical Component
Dictionary (CCD), baixadas para uma pasta temporaria. Os esquemas de ABO, LDL e
escudo de glicanos sao ilustrativos e dizem isso no titulo.

Uso: uv run --no-project --with matplotlib --with numpy python tools/make_biomolecule_figures.py [--cache DIR]
"""
from __future__ import annotations

import argparse
import math
import shlex
import sys
import tempfile
import urllib.request
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Polygon, Rectangle, RegularPolygon

sys.path.insert(0, str(Path(__file__).resolve().parent))
from make_figures import C, INK, BAND, ROOT, LIPID, save, frame_off, _lip, _chol  # noqa: E402

GREY = '#9a9a9a'
# cores CPK (Corey e Pauling, 1953; Koltun, 1965), um pouco atenuadas
CPK = {'C': '#8c8c8c', 'H': '#ffffff', 'N': '#3b5bdb', 'O': '#d9433b', 'P': '#ef8a17', 'S': '#e7c32a'}
STICK_H = '#cfcfcf'
# raios de van der Waals de Bondi (1964), em angstrom
VDW = {'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52, 'P': 1.80, 'S': 1.80}
SS_COL = {'H': C[3], 'E': C[0], '-': GREY}

PDB_IDS = ['1UBQ', '1BNA', '4HHB', '1ANF', '1CAG', '1N4E', '3GM7', '1EHZ']
CCD_IDS = ['PHE', 'GLC', 'PC1', 'CLR', 'DA', 'DT', 'U', 'ALA']
CACHE = None


# ---------------------------------------------------------------- leitura
def fetch(name, url):
    path = CACHE / name
    if not path.exists():
        urllib.request.urlretrieve(url, path)
    return path


def read_pdb(code, chains=None, hetero=False, model=1):
    """Atomos pesados (altloc vazio ou A) do primeiro modelo, sem agua."""
    path = fetch(f'{code}.pdb', f'https://files.rcsb.org/download/{code}.pdb')
    atoms, cur = [], 1
    for line in path.read_text().splitlines():
        if line.startswith('MODEL'):
            cur = int(line.split()[1])
        if cur != model or not line.startswith(('ATOM', 'HETATM')):
            continue
        if line.startswith('HETATM') and not hetero:
            continue
        resn = line[17:20].strip()
        if resn == 'HOH' or line[16] not in ' A':
            continue
        if chains and line[21] not in chains:
            continue
        elem = line[76:78].strip() or line[12:16].strip()[0]
        atoms.append(dict(name=line[12:16].strip(), resn=resn, chain=line[21], resi=int(line[22:26]), elem=elem.capitalize(),
                          xyz=np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])))
    return atoms


def read_course_pdb(path):
    """PDB coarse-grained do proprio curso (results/viz), primeiro modelo."""
    atoms = []
    for line in Path(path).read_text().splitlines():
        if line.startswith('ENDMDL'):
            break
        if line.startswith(('ATOM', 'HETATM')):
            atoms.append(dict(name=line[12:16].strip(), resn=line[17:21].strip(), chain=line[21], resi=int(line[22:26]),
                              xyz=np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])))
    return atoms


def read_ccd(code):
    """Atomos (nome, elemento, coordenada ideal) e ligacoes (i, j, ordem) do CCD."""
    path = fetch(f'{code}.cif', f'https://files.rcsb.org/ligands/view/{code}.cif')
    lines = path.read_text().splitlines()

    def loop(prefix):
        keys, rows, i = [], [], 0
        while i < len(lines) and not lines[i].startswith(prefix):
            i += 1
        while i < len(lines) and lines[i].startswith(prefix):
            keys.append(lines[i].split('.')[1].strip()); i += 1
        while i < len(lines) and not lines[i].startswith('#'):
            rows.append(dict(zip(keys, shlex.split(lines[i])))); i += 1
        return rows

    atoms = [dict(name=r['atom_id'], elem=r['type_symbol'].capitalize(),
                  xyz=np.array([float(r[f'pdbx_model_Cartn_{k}_ideal']) for k in 'xyz'])) for r in loop('_chem_comp_atom.')]
    index = {a['name']: k for k, a in enumerate(atoms)}
    order = {'SING': 1, 'DOUB': 2, 'TRIP': 3}
    bonds = [(index[r['atom_id_1']], index[r['atom_id_2']], order.get(r['value_order'], 1)) for r in loop('_chem_comp_bond.')]
    return atoms, bonds


def dist_bonds(atoms, cut=1.95):
    X = np.array([a['xyz'] for a in atoms])
    out = []
    for i in range(len(X)):
        d = np.linalg.norm(X[i + 1:] - X[i], axis=1)
        for j in np.nonzero(d < cut)[0]:
            out.append((i, i + 1 + j, 1))
    return out


def sec_struct(code):
    """Helices e fitas dos registros HELIX e SHEET: {(cadeia, resid): 'H' ou 'E'}."""
    ss = {}
    for line in fetch(f'{code}.pdb', '').read_text().splitlines():
        if line.startswith('HELIX'):
            for r in range(int(line[21:25]), int(line[33:37]) + 1):
                ss[(line[19], r)] = 'H'
        elif line.startswith('SHEET'):
            for r in range(int(line[22:26]), int(line[33:37]) + 1):
                ss[(line[21], r)] = 'E'
    return ss


# ---------------------------------------------------------------- geometria e desenho
def frame(X, rot=0.0, flip=False, axes=None):
    """Eixos principais de X: (centro, base 3x3). x = maior variancia, z = profundidade."""
    c = X.mean(axis=0)
    if axes is None:
        _, _, vt = np.linalg.svd(X - c)
        axes = vt
    a = math.radians(rot)
    R = np.array([[math.cos(a), -math.sin(a), 0], [math.sin(a), math.cos(a), 0], [0, 0, 1]])
    B = R @ axes
    if flip:
        B = B * np.array([[1], [-1], [-1]])
    return c, B


def best_view(atoms, scale=1.0, n=1500, seed=0, normal=None, max_tilt=40):
    """Rotacao 3D em que as esferas (raios do ballstick) menos se sobrepoem na projecao."""
    X = np.array([a['xyz'] for a in atoms])
    c = X.mean(0)
    r = np.array([(0.22 if a['elem'] == 'H' else 0.36) * scale for a in atoms])
    iu = np.triu_indices(len(X), 1)
    rsum = (r[:, None] + r[None, :])[iu]
    rng = np.random.default_rng(seed)
    best, bestB = -1e9, None
    for _ in range(n):
        q = rng.normal(size=4); q /= np.linalg.norm(q)
        w, x, y, z = q
        B = np.array([[1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
                      [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
                      [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)]])
        if normal is not None and abs(B[2] @ normal) < math.cos(math.radians(max_tilt)):
            continue
        P = (X - c) @ B.T
        d = np.linalg.norm(P[:, None, :2] - P[None, :, :2], axis=-1)[iu] - rsum
        s = np.sort(d)[:4].sum()          # as quatro piores sobreposicoes
        if s > best:
            best, bestB = s, B
    # maior extensao na horizontal
    P = (X - c) @ bestB.T
    _, _, vt = np.linalg.svd(P[:, :2] - P[:, :2].mean(0))
    a = math.atan2(vt[0, 1], vt[0, 0])
    Rz = np.array([[math.cos(-a), -math.sin(-a), 0], [math.sin(-a), math.cos(-a), 0], [0, 0, 1]])
    return c, Rz @ bestB


def frame_y(X, yvec):
    """Frame com y ao longo de `yvec` e x na maior variancia restante."""
    y = yvec / np.linalg.norm(yvec)
    c = X.mean(0)
    Q = (X - c) - np.outer((X - c) @ y, y)
    _, _, vt = np.linalg.svd(Q)
    x = vt[0] - (vt[0] @ y) * y; x /= np.linalg.norm(x)
    return c, np.array([x, y, np.cross(x, y)])


def upright(atoms, X, top, bottom):
    """Frame com o atomo `top` em cima e a media dos atomos `bottom` embaixo."""
    idx = {a['name']: k for k, a in enumerate(atoms)}
    y = X[idx[top]] - X[[idx[b] for b in bottom]].mean(0)
    y /= np.linalg.norm(y)
    c = X.mean(0)
    _, _, vt = np.linalg.svd((X - c) - np.outer((X - c) @ y, y))
    x = vt[0] - (vt[0] @ y) * y; x /= np.linalg.norm(x)
    return c, np.array([x, y, np.cross(x, y)])


def _rotate(V, axis, ang):
    axis = axis / np.linalg.norm(axis)
    return V * math.cos(ang) + np.cross(axis, V) * math.sin(ang) + np.outer(V @ axis, axis) * (1 - math.cos(ang))


def straight_tails(atoms):
    """DSPC do CCD com as duas caudas giradas como corpos rigidos, em torno do carbono da carbonila,
    ate ficarem paralelas ao eixo cabeca-glicerol. Comprimentos de ligacao preservados."""
    X = np.array([a['xyz'] for a in atoms])
    idx = {a['name']: k for k, a in enumerate(atoms)}
    axis = X[idx['C2']] - X[idx['P']]
    axis /= np.linalg.norm(axis)
    for s in '23':
        carbons = [n for n in idx if n.startswith('C' + s) and len(n) == 3 and n != f'C{s}1']
        tail = set(idx[n] for n in carbons)
        tail |= {idx[n] for n in idx if n.startswith('H' + s) and len(n) == 4 and n[:3] != f'H{s}1'}
        piv = X[idx[f'C{s}1']]
        d = X[idx[f'C{s}I']] - piv
        d /= np.linalg.norm(d)
        rot_axis = np.cross(d, axis)
        ang = math.asin(min(1.0, np.linalg.norm(rot_axis)))
        if d @ axis < 0:
            ang = math.pi - ang
        tl = sorted(tail)
        X[tl] = _rotate(X[tl] - piv, rot_axis, ang) + piv
    return X


def project(X, fr):
    c, B = fr
    return (X - c) @ B.T


def kabsch(P, Q):
    """Rotacao e translacao que levam P sobre Q."""
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    return lambda X: (X - pc) @ R.T + qc


def draw_molecule(ax, atoms, bonds, P, style='ballstick', scale=1.0, hide_h=False, z0=0):
    """Desenha com ordenacao por profundidade. style: sticks, ballstick, spacefill, lines."""
    keep = [k for k, a in enumerate(atoms) if not (hide_h and a['elem'] == 'H')]
    keepset = set(keep)
    z = P[:, 2]
    zr = {k: r for r, k in enumerate(sorted(keep, key=lambda k: z[k]))}
    n = len(keep) + 1
    if style == 'spacefill':
        for k in sorted(keep, key=lambda k: z[k]):
            e = atoms[k]['elem']
            ax.add_patch(Circle(P[k, :2], VDW.get(e, 1.7) * scale, fc=CPK.get(e, GREY), ec='#4a4a4a', lw=0.5,
                                zorder=z0 + 2 + zr[k] / n))
        return
    lw = {'sticks': 4.2, 'ballstick': 1.8, 'lines': 1.1}[style] * scale
    for i, j, order in bonds:
        if i not in keepset or j not in keepset:
            continue
        m = (P[i, :2] + P[j, :2]) / 2
        zz = z0 + 2 + min(zr[i], zr[j]) / n - 0.001
        for k, end in ((i, m), (j, m)):
            e = atoms[k]['elem']
            col = (STICK_H if e == 'H' else CPK.get(e, GREY)) if style != 'ballstick' else '#5a5a5a'
            if style == 'ballstick' and order == 2:
                d = P[j, :2] - P[i, :2]
                off = np.array([-d[1], d[0]]) / (np.linalg.norm(d) + 1e-9) * 0.13 * scale
                for s in (1, -1):
                    ax.plot([P[k, 0] + s * off[0], end[0] + s * off[0]], [P[k, 1] + s * off[1], end[1] + s * off[1]],
                            color=col, lw=lw * 0.8, solid_capstyle='round', zorder=zz)
            else:
                ax.plot([P[k, 0], end[0]], [P[k, 1], end[1]], color=col, lw=lw, solid_capstyle='round', zorder=zz)
    if style == 'ballstick':
        for k in keep:
            e = atoms[k]['elem']
            r = (0.22 if e == 'H' else 0.36) * scale
            ax.add_patch(Circle(P[k, :2], r, fc=CPK.get(e, GREY), ec='#4a4a4a', lw=0.6, zorder=z0 + 2 + zr[k] / n))


def draw_skeletal(ax, atoms, bonds, P, scale=1.0):
    """Formula de esqueleto: carbonos sao vertices, heteroatomos rotulados, H so em heteroatomo."""
    heavy = [k for k, a in enumerate(atoms) if a['elem'] != 'H']
    hset = set(heavy)
    nH = {k: 0 for k in heavy}
    for i, j, _ in bonds:
        if i in hset and j not in hset:
            nH[i] += 1
        elif j in hset and i not in hset:
            nH[j] += 1
    hetero = {k for k in heavy if atoms[k]['elem'] != 'C'}
    if len(heavy) == len(atoms):    # estrutura sem H (cristal): valencia completa com H implicito
        deg = {k: 0 for k in heavy}
        for i, j, _ in bonds:
            deg[i] += 1; deg[j] += 1
        val = {'O': 2, 'N': 3, 'S': 2}
        nH = {k: max(0, val.get(atoms[k]['elem'], 0) - deg[k]) if k in hetero else 0 for k in heavy}
    for i, j, order in bonds:
        if i not in hset or j not in hset:
            continue
        a, b = P[i, :2].copy(), P[j, :2].copy()
        d = b - a
        u = d / np.linalg.norm(d)
        if i in hetero:
            a = a + u * 0.38 * scale
        if j in hetero:
            b = b - u * 0.38 * scale
        ax.plot(*np.c_[a, b], color=INK, lw=1.5, solid_capstyle='round', zorder=2)
        if order == 2:
            off = np.array([-u[1], u[0]]) * 0.17 * scale
            ax.plot(*np.c_[a + off + u * 0.12, b + off - u * 0.12], color=INK, lw=1.5, zorder=2)
    for k in hetero:
        e = atoms[k]['elem']
        lab = e + ('H' if nH[k] == 1 else 'H' + '₀₁₂₃'[nH[k]] if nH[k] > 1 else '')
        ax.text(*P[k, :2], lab, ha='center', va='center', fontsize=9 * scale, color=CPK[e] if e != 'C' else INK,
                fontweight='bold', zorder=3)


def fit(ax, P, pad=1.2):
    lo, hi = P[:, :2].min(0) - pad, P[:, :2].max(0) + pad
    c, half = (lo + hi) / 2, (hi - lo).max() / 2
    ax.set_xlim(c[0] - half, c[0] + half); ax.set_ylim(c[1] - half, c[1] + half)
    ax.set_aspect('equal'); frame_off(ax)


def fit_rect(ax, P, pad=1.2):
    lo, hi = P[:, :2].min(0) - pad, P[:, :2].max(0) + pad
    ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1]); ax.set_aspect('equal'); frame_off(ax)


def tag(ax, xy, text, xytext=None, color=INK, fs=8.5, box=True):
    """Registra uma etiqueta; a posicao final sai de place_tags, depois do enquadramento."""
    pref = None if xytext is None else np.subtract(xytext, xy)
    ax.__dict__.setdefault('_tags', []).append((np.asarray(xy, float), text, pref, color, fs, box))


def _obstacles(ax):
    """Circulos (centro e raio em pixels) de tudo o que esta desenhado no eixo."""
    T = ax.transData
    obs = []
    for p in ax.patches:
        if isinstance(p, Circle):
            c = T.transform(p.center)
            r = abs(T.transform((p.center[0] + p.radius, p.center[1]))[0] - c[0])
            obs.append((c, r))
        else:
            ext = p.get_window_extent()
            if ext.width < 400 and ext.height < 400:
                for xx in np.linspace(ext.x0, ext.x1, 6):
                    for yy in np.linspace(ext.y0, ext.y1, 6):
                        obs.append((np.array([xx, yy]), 3.0))
    for ln in ax.lines:
        x, y = ln.get_data()
        if len(x) < 2:
            continue
        pts = T.transform(np.c_[x, y])
        r = max(ln.get_linewidth(), 1.0) * ax.figure.dpi / 72 / 2
        for a, b in zip(pts[:-1], pts[1:]):
            n = max(2, int(np.linalg.norm(b - a) / 3))
            for s in np.linspace(0, 1, n):
                obs.append((a + s * (b - a), r))
    for col in ax.collections:
        for c in col.get_offsets():
            obs.append((T.transform(c), 3.5))
    return obs


def _seg_clear(a, b, obs, skip, near=0.5):
    """Segmento livre de obstaculos, exceto os que estao a menos de `near` px da ancora (atomo e vizinhos)."""
    d = b - a
    L2 = d @ d + 1e-9
    for c, r in obs:
        if np.linalg.norm(c - skip) < near:
            continue
        s = np.clip((c - a) @ d / L2, 0, 1)
        if np.linalg.norm(a + s * d - c) < r + 1.5:
            return False
    return True


def place_tags(ax):
    tags = ax.__dict__.get('_tags', [])
    if not tags:
        return
    fig = ax.figure
    rend = fig.canvas.get_renderer()
    obs = _obstacles(ax)
    T, Ti = ax.transData, ax.transData.inverted()
    axbox = ax.get_window_extent()
    placed = []
    for xy, text, pref, color, fs, box in tags:
        anchor = T.transform(xy)
        probe = ax.text(0, 0, text, fontsize=fs, bbox=dict(boxstyle='round,pad=0.25'))
        ext = probe.get_window_extent(rend)
        probe.remove()
        w, h = ext.width + 8, ext.height + 8
        # atomo ancora: o obstaculo cujo centro coincide com o ponto
        own = min(obs, key=lambda o: np.linalg.norm(o[0] - anchor))
        own_r = own[1] if np.linalg.norm(own[0] - anchor) < 1 else 0
        start_ang = math.atan2(pref[1], pref[0]) if pref is not None else math.pi / 2
        best = None
        for dist in range(int(own_r + h / 2 + 10), 420, 8):
            for k in range(24):
                ang = start_ang + (k + 1) // 2 * (-1) ** k * math.pi / 12
                cx = anchor[0] + math.cos(ang) * (dist + w / 2 * abs(math.cos(ang)))
                cy = anchor[1] + math.sin(ang) * (dist + h / 2 * abs(math.sin(ang)))
                x0, x1, y0, y1 = cx - w / 2, cx + w / 2, cy - h / 2, cy + h / 2
                if x0 < axbox.x0 - 30 or x1 > axbox.x1 + 30 or y0 < axbox.y0 - 10 or y1 > axbox.y1 + 10:
                    continue
                if any(x0 - r < c[0] < x1 + r and y0 - r < c[1] < y1 + r for c, r in obs):
                    continue
                if any(not (x1 < q[0] or x0 > q[1] or y1 < q[2] or y0 > q[3]) for q in placed):
                    continue
                # ponto da caixa mais perto da ancora: fim da linha guia
                end = np.array([np.clip(anchor[0], x0, x1), np.clip(anchor[1], y0, y1)])
                u = (end - anchor) / (np.linalg.norm(end - anchor) + 1e-9)
                if not _seg_clear(anchor + u * (own_r + 1), end, obs, anchor, near=max(own_r * 2.2, 6)):
                    continue
                if any(_seg_cross_box(anchor, end, q) for q in placed):
                    continue
                best = (cx, cy, x0, x1, y0, y1, end, u)
                break
            if best:
                break
        if best is None:
            raise RuntimeError(f'sem lugar livre para a etiqueta {text!r}')
        cx, cy, x0, x1, y0, y1, end, u = best
        placed.append((x0, x1, y0, y1))
        ax.text(*Ti.transform((cx, cy)), text, ha='center', va='center', fontsize=fs, color=color, zorder=60,
                fontweight='normal' if box else 'bold',
                bbox=dict(boxstyle='round,pad=0.25', fc='white', ec=color if box else 'none', lw=0.8))
        ax.plot(*np.c_[Ti.transform(anchor + u * (own_r + 1)), Ti.transform(end)], color=color if box else GREY,
                lw=0.9 if box else 0.7, zorder=59)
        obs.extend((p, 2.0) for p in np.linspace(anchor, end, 12))
    ax._tags = []


def _seg_cross_box(a, b, q):
    x0, x1, y0, y1 = q
    for s in np.linspace(0, 1, 25):
        p = a + s * (b - a)
        if x0 <= p[0] <= x1 and y0 <= p[1] <= y1:
            return True
    return False


def done(fig, name):
    fig.canvas.draw()
    for ax in fig.axes:
        place_tags(ax)
    save(fig, name)


def element_legend(ax, elems, loc='lower center', ncol=None, anchor=None):
    hs = [plt.Line2D([], [], marker='o', ls='', ms=8, mfc=CPK[e], mec='#4a4a4a', label=NAMES[e]) for e in elems]
    kw = dict(handles=hs, loc=loc, ncol=ncol or len(elems), frameon=False, fontsize=8.5, handletextpad=0.3, columnspacing=1.0)
    if anchor:
        kw['bbox_to_anchor'] = anchor
    ax.legend(**kw)


NAMES = {'H': 'hidrogênio', 'C': 'carbono', 'N': 'nitrogênio', 'O': 'oxigênio', 'P': 'fósforo', 'S': 'enxofre'}


def trace(ax, pts, cols, lw=3.0, z=2):
    for k in range(len(pts) - 1):
        ax.plot(*pts[k:k + 2, :2].T, color=cols[k], lw=lw, solid_capstyle='round', zorder=z + pts[k, 2] * 1e-3)


# ---------------------------------------------------------------- 1. os atomos
def fig_atoms():
    prot = read_pdb('1UBQ')
    dna = read_pdb('1BNA')
    cnt = lambda at, e: sum(a['elem'] == e for a in at)
    info = [('H', 1, ['em todas as biomoléculas']),
            ('C', 4, ['esqueleto de todas', 'as biomoléculas']),
            ('N', 3, ['grupo amino,', 'bases nitrogenadas']),
            ('O', 2, ['água, açúcares,', 'carbonilas, fosfatos']),
            ('P', 4, ['fosfatos do DNA, RNA,', 'fosfolipídios e ATP']),
            ('S', 2, ['cisteína e metionina,', 'pontes dissulfeto'])]
    W, GAP, S = 4.0, 0.45, 0.78          # largura do cartao, espaco entre cartoes, escala das esferas
    fig, ax = plt.subplots(figsize=(15.5, 5.6))
    frame_off(ax); ax.set_aspect('equal')
    ax.set_xlim(-W / 2 - 0.1, 5 * (W + GAP) + W / 2 + 0.1); ax.set_ylim(-6.9, 2.15)
    for k, (e, val, where) in enumerate(info):
        x = k * (W + GAP)
        ax.add_patch(FancyBboxPatch((x - W / 2, -6.75), W, 8.75, boxstyle='round,pad=0,rounding_size=0.25',
                                    fc='white', ec=BAND, lw=1.0, zorder=0))
        ax.add_patch(Circle((x, 0.5), VDW[e] * S, fc=CPK[e], ec='#4a4a4a', lw=0.9, zorder=2))
        ax.text(x, 0.5, e, ha='center', va='center', fontsize=17, fontweight='bold',
                color=INK if e in 'HS' else 'white', zorder=3)
        ax.text(x, -1.55, NAMES[e], ha='center', va='center', fontsize=11, fontweight='bold')
        rows = [('vizinhos', f'{val}'), ('raio de vdW', f'{VDW[e]:.2f} Å'.replace('.', ','))]
        if e == 'H':
            rows.append(('ubiquitina', '629 de 1 231'))
            rows.append(('nos cristais', 'não aparece'))
        else:
            rows.append(('1UBQ', str(cnt(prot, e))))
            rows.append(('1BNA', str(cnt(dna, e))))
        for r, (lab, v) in enumerate(rows):
            y = -2.35 - 0.52 * r
            ax.text(x - W / 2 + 0.3, y, lab, ha='left', va='center', fontsize=9, color='#5a5a5a')
            ax.text(x + W / 2 - 0.3, y, v, ha='right', va='center', fontsize=9, fontweight='bold')
        ax.plot([x - W / 2 + 0.3, x + W / 2 - 0.3], [-4.55, -4.55], color=BAND, lw=0.8)
        ax.text(x, -5.15, '\n'.join(where), ha='center', va='top', fontsize=9, linespacing=1.35)
    fig.suptitle('Os seis elementos das biomoléculas: cor CPK, quantos vizinhos costumam ter, raio de van der Waals '
                 '(Bondi, 1964) e quantos há nos cristais 1UBQ e 1BNA', fontsize=11, y=0.96)
    done(fig, 'bio_atomos.png')


# ---------------------------------------------------------------- 2. as quatro classes
def fig_classes():
    fig, axes = plt.subplots(2, 4, figsize=(15, 7.6), gridspec_kw=dict(height_ratios=[1, 1.15]))
    mono = [('ALA', 'Aminoácido: alanina', 0), ('DA', 'Nucleotídeo: dAMP', 0), ('PC1', 'Fosfolipídio: DSPC\n(caudas alinhadas, sem H)', 0),
            ('GLC', 'Monossacarídeo: glicose', 0)]
    for ax, (code, title, rot) in zip(axes[0], mono):
        atoms, bonds = read_ccd(code)
        X = np.array([a['xyz'] for a in atoms]) if code != 'PC1' else straight_tails(atoms)
        P = project(X, best_view(atoms)) if code != 'PC1' else project(X, upright(atoms, X, 'P', ('C2I', 'C3I')))
        sc = 0.75 if code == 'PC1' else 1.0
        draw_molecule(ax, atoms, bonds, P, 'ballstick', scale=1.2 if code == 'PC1' else sc, hide_h=code == 'PC1')
        fit(ax, P, pad=0.9)
        ax.set_title(title, fontsize=10.5)
    # polimeros e agregados
    ub = [a for a in read_pdb('1UBQ') if a['name'] == 'CA']
    ss = sec_struct('1UBQ')
    X = np.array([a['xyz'] for a in ub]); P = project(X, frame(X))
    trace(axes[1, 0], P, [SS_COL[ss.get(('A', a['resi']), '-')] for a in ub], lw=3.2)
    fit(axes[1, 0], P, pad=2); axes[1, 0].set_title('Proteína: ubiquitina\n76 resíduos (1UBQ)', fontsize=10.5)
    hs = [plt.Line2D([], [], color=SS_COL[k], lw=3, label=l) for k, l in (('H', 'hélice α'), ('E', 'fita β'), ('-', 'alça'))]
    axes[1, 0].legend(handles=hs, loc='lower center', bbox_to_anchor=(0.5, -0.12), ncol=3, frameon=False, fontsize=8.5)

    dna = read_pdb('1BNA')
    X = np.array([a['xyz'] for a in dna]); P = project(X, frame(X, rot=90))
    draw_molecule(axes[1, 1], dna, dist_bonds(dna), P, 'lines', scale=1.0)
    fit(axes[1, 1], P, pad=1.5); axes[1, 1].set_title('Ácido nucleico: DNA\n12 pares de bases (1BNA)', fontsize=10.5)

    mem = read_course_pdb(ROOT / 'results/viz/rep1_mem_final.pdb')
    X = np.array([a['xyz'] for a in mem])
    head = {'NC3': C[0], 'PO4': C[5], 'GL1': C[4], 'GL2': C[4]}
    # a normal da bicamada e o eixo y; fatia fina em z, vista de lado (x contra y)
    sel = X[:, 2] > np.percentile(X[:, 2], 75)
    for r in {a['resi'] for a in mem}:
        idx = [k for k, a in enumerate(mem) if a['resi'] == r]
        if np.ptp(X[idx, 0]) > 25 or np.ptp(X[idx, 1]) > 35:
            sel[idx] = False
    heads = np.array([n['name'] in head for n in mem])
    lo, hi = X[sel & heads, 0].min(), X[sel & heads, 0].max()
    sel &= (X[:, 0] >= lo - 1) & (X[:, 0] <= hi + 1)
    order = np.argsort(X[:, 2])
    for k in order:
        if not sel[k]:
            continue
        axes[1, 2].add_patch(Circle((X[k, 0], X[k, 1]), 2.3, fc=head.get(mem[k]['name'], '#d9c98f'), ec='#6a6a6a',
                                    lw=0.3, zorder=2 + X[k, 2] * 1e-3))
    fit(axes[1, 2], X[sel][:, :2], pad=4)
    axes[1, 2].set_title('Agregado: bicamada de DSPC\n(Parte B desta aula, Martini)', fontsize=10.5)
    hs = [plt.Line2D([], [], marker='o', ls='', ms=8, mfc=c, mec='#6a6a6a', label=l)
          for c, l in ((C[0], 'colina'), (C[5], 'fosfato'), (C[4], 'glicerol'), ('#d9c98f', 'caudas'))]
    axes[1, 2].legend(handles=hs, loc='lower center', bbox_to_anchor=(0.5, -0.12), ncol=4, frameon=False, fontsize=8.5,
                      handletextpad=0.2, columnspacing=0.8)

    mal = read_pdb('1ANF', chains='B', hetero=True)
    mal = [a for a in mal if a['resn'] == 'GLC']
    X = np.array([a['xyz'] for a in mal]); P = project(X, frame(X))
    draw_molecule(axes[1, 3], mal, dist_bonds(mal, 1.7), P, 'ballstick', scale=1.0)
    o4 = next(k for k, a in enumerate(mal) if a['name'] == 'O4' and a['resi'] == 1)
    tag(axes[1, 3], P[o4, :2], 'ligação glicosídica α(1-4)', P[o4, :2] + np.array([0, -3.6]), color=C[3])
    fit(axes[1, 3], P, pad=2.2); axes[1, 3].set_title('Carboidrato: maltose\nduas glicoses (1ANF)', fontsize=10.5)
    element_legend(axes[0, 1], ['H', 'C', 'N', 'O', 'P'], anchor=(1.04, -0.1))
    fig.suptitle('Quatro classes de biomoléculas: a unidade (coordenadas ideais do CCD) e a molécula que ela forma '
                 '(estruturas do PDB e simulação desta aula)', fontsize=11, y=0.99)
    fig.subplots_adjust(hspace=0.32, wspace=0.08)
    done(fig, 'bio_classes.png')


# ---------------------------------------------------------------- 3. niveis de estrutura da proteina
def fig_protein_levels():
    ub = read_pdb('1UBQ')
    ss = sec_struct('1UBQ')
    fig, axes = plt.subplots(1, 4, figsize=(16, 5.0), gridspec_kw=dict(width_ratios=[1.15, 1, 1, 1.1]))
    # primaria: Met1-Gln2-Ile3, ligacoes peptidicas
    tri = [a for a in ub if a['resi'] <= 3]
    X = np.array([a['xyz'] for a in tri]); P = project(X, frame(X, rot=0))
    b = dist_bonds(tri)
    draw_molecule(axes[0], tri, b, P, 'ballstick', scale=1.0)
    dpep = []
    for r in (1, 2):
        i = next(k for k, a in enumerate(tri) if a['resi'] == r and a['name'] == 'C')
        j = next(k for k, a in enumerate(tri) if a['resi'] == r + 1 and a['name'] == 'N')
        axes[0].plot(*P[[i, j], :2].T, color=C[3], lw=11, alpha=0.55, zorder=1.5, solid_capstyle='round')
        dpep.append(np.linalg.norm(X[i] - X[j]))
    ctr = P[:, :2].mean(0)
    for k, a in enumerate(tri):
        if a['name'] == 'CA':
            tag(axes[0], P[k, :2], f'{a["resn"].title()}{a["resi"]}', P[k, :2] + (P[k, :2] - ctr), box=False, fs=9)
    fit(axes[0], P, pad=2.4)
    seq = 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'
    axes[0].set_title('Primária: a sequência dos 76 resíduos', fontsize=9.5, pad=34)
    axes[0].text(0.5, 1.02, seq[:38] + '\n' + seq[38:], transform=axes[0].transAxes, ha='center', va='bottom',
                 fontsize=8, family='monospace')
    pep = plt.Line2D([], [], color=C[3], lw=8, alpha=0.55,
                     label=f'ligação peptídica C-N: {dpep[0]:.2f} e {dpep[1]:.2f} Å'.replace('.', ','))
    hs = [plt.Line2D([], [], marker='o', ls='', ms=8, mfc=CPK[e], mec='#4a4a4a', label=NAMES[e]) for e in 'CNOS']
    lg = axes[0].legend(handles=hs, loc='upper center', bbox_to_anchor=(0.5, 0.02), ncol=4, frameon=False,
                        fontsize=8.5, handletextpad=0.3, columnspacing=0.9)
    axes[0].add_artist(lg)
    axes[0].legend(handles=[pep], loc='upper center', bbox_to_anchor=(0.5, -0.06), frameon=False, fontsize=8.5,
                   handlelength=2.2)
    # secundaria: helice 23-34 com ligacoes de hidrogenio O(i)-N(i+4)
    hel = [a for a in ub if 23 <= a['resi'] <= 34 and a['name'] in ('N', 'CA', 'C', 'O')]
    X = np.array([a['xyz'] for a in hel]); P = project(X, frame(X))
    draw_molecule(axes[1], hel, dist_bonds(hel), P, 'ballstick', scale=1.0)
    nhb = []
    for i, a in enumerate(hel):
        if a['name'] != 'O':
            continue
        for j, c in enumerate(hel):
            if c['name'] == 'N' and c['resi'] == a['resi'] + 4:
                d = np.linalg.norm(X[i] - X[j])
                if d < 3.5:
                    nhb.append(d)
                    axes[1].plot(*P[[i, j], :2].T, color=C[0], lw=1.2, ls=(0, (2, 2)), zorder=1)
    fit(axes[1], P, pad=1.6)
    axes[1].set_title(f'Secundária: hélice α 23-34\n{len(nhb)} ligações de hidrogênio O(i)···N(i+4), '
                      f'{min(nhb):.1f} a {max(nhb):.1f} Å'.replace('.', ','), fontsize=9.5)
    axes[1].plot([], [], color=C[0], lw=1.2, ls=(0, (2, 2)), label='ligação de hidrogênio')
    axes[1].legend(loc='lower center', bbox_to_anchor=(0.5, -0.1), frameon=False, fontsize=8.5)
    # terciaria
    ca = [a for a in ub if a['name'] == 'CA']
    X = np.array([a['xyz'] for a in ca]); P = project(X, frame(X))
    trace(axes[2], P, [SS_COL[ss.get(('A', a['resi']), '-')] for a in ca], lw=3.4)
    fit(axes[2], P, pad=2)
    axes[2].set_title('Terciária: a cadeia inteira dobrada\nubiquitina (1UBQ), traço do Cα', fontsize=9.5)
    hs = [plt.Line2D([], [], color=SS_COL[k], lw=3, label=l) for k, l in (('H', 'hélice α'), ('E', 'fita β'), ('-', 'alça'))]
    axes[2].legend(handles=hs, loc='lower center', bbox_to_anchor=(0.5, -0.1), ncol=3, frameon=False, fontsize=8.5)
    # quaternaria
    hb = read_pdb('4HHB', hetero=True)
    ca = [a for a in hb if a['name'] == 'CA' and a['resn'] != 'HEM']
    hem = [a for a in hb if a['resn'] == 'HEM']
    X = np.array([a['xyz'] for a in ca]); fr = frame(X); P = project(X, fr)
    col = {'A': C[0], 'B': C[3], 'C': C[2], 'D': C[5]}
    for ch in 'ABCD':
        idx = [k for k, a in enumerate(ca) if a['chain'] == ch]
        trace(axes[3], P[idx], [col[ch]] * len(idx), lw=1.8)
    H = project(np.array([a['xyz'] for a in hem]), fr)
    axes[3].scatter(H[:, 0], H[:, 1], s=6, color=INK, zorder=50)
    fit(axes[3], P, pad=2)
    axes[3].set_title('Quaternária: várias cadeias\nhemoglobina (4HHB), α₂β₂ com quatro hemes', fontsize=9.5)
    hs = [plt.Line2D([], [], color=col[k], lw=3, label=l) for k, l in (('A', 'α1'), ('B', 'β1'), ('C', 'α2'), ('D', 'β2'))]
    hs.append(plt.Line2D([], [], marker='o', ls='', ms=4, color=INK, label='heme'))
    axes[3].legend(handles=hs, loc='lower center', bbox_to_anchor=(0.5, -0.1), ncol=5, frameon=False, fontsize=8.5,
                   columnspacing=0.8)
    fig.suptitle('Os quatro níveis de estrutura de uma proteína (coordenadas cristalográficas, átomos pesados)', fontsize=11, y=1.02)
    done(fig, 'bio_proteina_niveis.png')


# ---------------------------------------------------------------- 4. colageno: Gly trocada por Ala
def fig_collagen():
    at = read_pdb('1CAG')
    X = np.array([a['xyz'] for a in at])
    ca = np.array([a['xyz'] for a in at if a['name'] == 'CA'])
    fr = frame(ca, rot=0)
    P = project(X, fr)
    fig, axes = plt.subplots(2, 1, figsize=(14, 7.2), gridspec_kw=dict(height_ratios=[1, 1.2]))
    ax = axes[0]
    col = {'A': C[0], 'B': C[2], 'C': C[5]}
    for ch in 'ABC':
        idx = [k for k, a in enumerate(at) if a['chain'] == ch and a['name'] == 'CA']
        trace(ax, P[idx], [col[ch]] * len(idx), lw=3)
    for k, a in enumerate(at):
        if a['name'] == 'CA' and a['resn'] in ('GLY', 'ALA'):
            ax.add_patch(Circle(P[k, :2], 0.75, fc='white' if a['resn'] == 'GLY' else C[3], ec=INK, lw=0.8, zorder=10))
    fit_rect(ax, P, pad=2.5)
    ax.set_title('Tripla hélice de um peptídeo de colágeno (1CAG): três cadeias (Pro-Hyp-Gly)₁₀, com Ala no lugar de uma Gly',
                 fontsize=10.5)
    hs = [plt.Line2D([], [], color=col[k], lw=3, label=f'cadeia {k}') for k in 'ABC']
    hs += [plt.Line2D([], [], marker='o', ls='', ms=8, mfc='white', mec=INK, label='Gly'),
           plt.Line2D([], [], marker='o', ls='', ms=8, mfc=C[3], mec=INK, label='Ala (substituição)')]
    ax.legend(handles=hs, loc='lower center', bbox_to_anchor=(0.5, -0.2), ncol=5, frameon=False, fontsize=9)
    # distancias N(Gly/Ala) ate o O mais proximo de outra cadeia
    ax = axes[1]
    rows = []
    for a in at:
        if a['name'] == 'N' and a['resn'] in ('GLY', 'ALA'):
            d = min(np.linalg.norm(a['xyz'] - b['xyz']) for b in at if b['name'] == 'O' and b['chain'] != a['chain'])
            rows.append((a['chain'], a['resi'], a['resn'], d))
    for ch_i, ch in enumerate('ABC'):
        rr = [r for r in rows if r[0] == ch]
        pos = [(r[1] - 1) % 30 + 1 for r in rr]
        ax.plot(pos, [r[3] for r in rr], color=col[ch], lw=1.2, marker='o', ms=0, label=f'cadeia {ch}')
        for p, r in zip(pos, rr):
            ax.scatter(p, r[3], s=46, color='white' if r[2] == 'GLY' else C[3], edgecolor=col[ch], lw=1.4, zorder=5)
    ax.axhspan(2.6, 3.2, color=BAND, alpha=0.5, lw=0, zorder=0, label='faixa de ligação de hidrogênio (2,6 a 3,2 Å)')
    ax.set_xlabel('posição na cadeia (1 a 30)'); ax.set_ylabel('N da Gly ou Ala até o O\nmais próximo de outra cadeia (Å)')
    ax.set_xlim(0, 31); ax.set_ylim(2.2, 5.0)
    ax.legend(loc='upper right', frameon=False, fontsize=9)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(True)
    fig.subplots_adjust(hspace=0.45)
    done(fig, 'bio_colageno.png')
    return rows


# ---------------------------------------------------------------- 5. DNA: pares de bases, fosfodiester, dimero de timina
BASE = {'N1', 'C2', 'N3', 'C4', 'C5', 'C6', 'N7', 'C8', 'N9', 'N6', 'O6', 'N2', 'O2', 'O4', 'N4', 'C7', "C1'"}


def _pair(ax, dna, ra, rb, hb, title):
    sel = [a for a in dna if ((a['chain'], a['resi']) in (ra, rb)) and a['name'] in BASE]
    X = np.array([a['xyz'] for a in sel]); P = project(X, frame(X))
    draw_molecule(ax, sel, dist_bonds(sel, 1.6), P, 'ballstick', scale=1.0)
    out = []
    for na, nb in hb:
        i = next(k for k, a in enumerate(sel) if (a['chain'], a['resi']) == ra and a['name'] == na)
        j = next(k for k, a in enumerate(sel) if (a['chain'], a['resi']) == rb and a['name'] == nb)
        d = np.linalg.norm(X[i] - X[j]); out.append(d)
        ax.plot(*P[[i, j], :2].T, color=C[0], lw=1.4, ls=(0, (2, 2)), zorder=1)
        m = P[[i, j], :2].mean(0)
        ax.text(m[0], m[1] + 0.32, f'{d:.2f}'.replace('.', ','), ha='center', va='bottom', fontsize=8, color=C[0],
                zorder=40, bbox=dict(boxstyle='round,pad=0.12', fc='white', ec='none'))
    ctr = P[:, :2].mean(0)
    for r in (ra, rb):
        k = next(k for k, a in enumerate(sel) if (a['chain'], a['resi']) == r and a['name'] == "C1'")
        resn = next(a['resn'] for a in sel if (a['chain'], a['resi']) == r)
        tag(ax, P[k, :2], f'{resn[1]}{r[1]}', P[k, :2] + (P[k, :2] - ctr), box=False, fs=9)
    fit(ax, P, pad=2.4)
    ax.set_title(title, fontsize=10)
    return out


def fig_dna():
    dna = read_pdb('1BNA')
    fig, axes = plt.subplots(1, 4, figsize=(16, 4.8))
    at = _pair(axes[0], dna, ('A', 5), ('B', 20), [('N6', 'O4'), ('N1', 'N3')], 'Par A·T: 2 ligações de hidrogênio')
    gc = _pair(axes[1], dna, ('A', 4), ('B', 21), [('O6', 'N4'), ('N1', 'N3'), ('N2', 'O2')], 'Par G·C: 3 ligações de hidrogênio')
    # fosfodiester entre dois nucleotideos da mesma fita
    nt = [a for a in dna if a['chain'] == 'A' and a['resi'] in (5, 6)]
    X = np.array([a['xyz'] for a in nt]); P = project(X, frame(X, rot=90))
    draw_molecule(axes[2], nt, dist_bonds(nt, 1.75), P, 'ballstick')
    p = next(k for k, a in enumerate(nt) if a['name'] == 'P' and a['resi'] == 6)
    tag(axes[2], P[p, :2], 'ligação fosfodiéster\n(açúcar-fosfato-açúcar)', P[p, :2] + np.array([3.2, 2.6]), color=C[5])
    fit(axes[2], P, pad=2.6)
    axes[2].set_title('Esqueleto: A5 e T6 ligados pelo fosfato', fontsize=10)
    element_legend(axes[2], ['C', 'N', 'O', 'P'], anchor=(0.5, -0.1))
    # dimero de timina (1N4E) contra duas timinas empilhadas normais (1BNA)
    ax = axes[3]
    dim = [a for a in read_pdb('1N4E', chains='B') if a['resi'] in (15, 16) and a['name'] in BASE]
    X = np.array([a['xyz'] for a in dim])
    c15 = X[[k for k, a in enumerate(dim) if a['resi'] == 15]].mean(0)
    c16 = X[[k for k, a in enumerate(dim) if a['resi'] == 16]].mean(0)
    P = project(X, frame_y(X, c16 - c15))
    b = dist_bonds(dim, 1.6)
    draw_molecule(ax, dim, b, P, 'ballstick')
    dd = []
    for nm in ('C5', 'C6'):
        i = next(k for k, a in enumerate(dim) if a['resi'] == 15 and a['name'] == nm)
        j = next(k for k, a in enumerate(dim) if a['resi'] == 16 and a['name'] == nm)
        dd.append(np.linalg.norm(X[i] - X[j]))
        ax.plot(*P[[i, j], :2].T, color=C[3], lw=5, alpha=0.6, zorder=1.5, solid_capstyle='round')
    m = P[[next(k for k, a in enumerate(dim) if a['resi'] == 15 and a['name'] == 'C5')], :2][0]
    tag(ax, m, f'anel de ciclobutano\nC5-C5 {dd[0]:.2f} Å, C6-C6 {dd[1]:.2f} Å'.replace('.', ','), m + np.array([4, 0]),
        color=C[3])
    for r in (15, 16):
        k = next(k for k, a in enumerate(dim) if a['resi'] == r and a['name'] == "C1'")
        tag(ax, P[k, :2], f'T{r}', P[k, :2] + np.array([-3, 0]), box=False, fs=9)
    fit(ax, P, pad=2.6)
    nrm = [a for a in dna if a['chain'] == 'A' and a['resi'] in (7, 8)]
    c5 = [a['xyz'] for a in nrm if a['name'] == 'C5']
    ax.set_title(f'Dímero de timina (1N4E): ligações novas\nentre T15 e T16; no DNA normal (1BNA, T7-T8),\n'
                 f'C5···C5 = {np.linalg.norm(c5[0] - c5[1]):.1f} Å'.replace('.', ','), fontsize=9.5)
    axes[0].plot([], [], color=C[0], lw=1.4, ls=(0, (2, 2)), label='ligação de hidrogênio (Å)')
    axes[0].legend(loc='lower center', bbox_to_anchor=(1.1, -0.12), frameon=False, fontsize=8.5)
    fig.suptitle('DNA: pareamento de bases e esqueleto no dodecâmero 1BNA, e a lesão causada pela luz UV no 1N4E '
                 '(distâncias medidas nos cristais)', fontsize=11, y=1.04)
    done(fig, 'bio_dna.png')
    return at, gc, dd


# ---------------------------------------------------------------- 6. RNA
def fig_rna():
    fig, axes = plt.subplots(1, 4, figsize=(16, 4.8), gridspec_kw=dict(width_ratios=[1, 1, 1.2, 1.2]))
    Ps = []
    for ax, code, title, marks in ((axes[0], 'DT', 'DNA: timidina (dTMP)', [("C2'", 'sem OH em 2\''), ('C7', 'metila')]),
                                   (axes[1], 'U', 'RNA: uridina (UMP)', [("O2'", 'OH em 2\''), ('C5', 'sem metila')])):
        atoms, bonds = read_ccd(code)
        X = np.array([a['xyz'] for a in atoms])
        P = project(X, best_view(atoms))
        draw_molecule(ax, atoms, bonds, P, 'ballstick')
        for (nm, lab), off in zip(marks, (np.array([2.6, -2.4]), np.array([-2.8, 2.2]))):
            k = next(k for k, a in enumerate(atoms) if a['name'] == nm)
            tag(ax, P[k, :2], lab, P[k, :2] + off, color=C[3])
        Ps.append((ax, P))
        ax.set_title(title, fontsize=10)
    half = max(np.ptp(P[:, :2], axis=0).max() for _, P in Ps) / 2 + 1.8
    for ax, P in Ps:
        c = (P[:, :2].min(0) + P[:, :2].max(0)) / 2
        ax.set_xlim(c[0] - half, c[0] + half); ax.set_ylim(c[1] - half, c[1] + half)
        ax.set_aspect('equal'); frame_off(ax)
    element_legend(axes[0], ['H', 'C', 'N', 'O', 'P'], anchor=(1.05, -0.12))
    # tRNA-Phe: fosfatos ligados em sequencia
    t = [a for a in read_pdb('1EHZ', hetero=True) if a['name'] == 'P']
    X = np.array([a['xyz'] for a in t]); P = project(X, frame(X, rot=0))
    cols = plt.get_cmap('viridis')(np.linspace(0, 1, len(t)))
    trace(axes[2], P, cols, lw=2.6)
    fit(axes[2], P, pad=3)
    axes[2].set_title('Saudável: tRNA-Phe de levedura (1EHZ)\numa fita de 76 nucleotídeos dobrada em L', fontsize=10)
    axes[2].text(0.5, -0.06, "traço dos fosfatos, cor de 5' (roxo) a 3' (amarelo)", transform=axes[2].transAxes,
                 ha='center', va='top', fontsize=8.5)
    # (CUG)6: pares U·U destacados
    cug = read_pdb('3GM7', chains='GH')
    X = np.array([a['xyz'] for a in cug]); P = project(X, frame(X, rot=90))
    draw_molecule(axes[3], cug, dist_bonds(cug, 1.75), P, 'lines')
    for k, a in enumerate(cug):
        if a['resn'] == 'U' and a['name'] == 'N3':
            axes[3].add_patch(Circle(P[k, :2], 1.1, fc=C[3], ec=INK, lw=0.6, alpha=0.85, zorder=40))
    fit(axes[3], P, pad=2)
    axes[3].set_title('Patológico: repetições (CUG)₆ (3GM7)\npareadas consigo mesmas, com U·U', fontsize=10)
    axes[3].scatter([], [], s=60, color=C[3], edgecolor=INK, label='uracila do par U·U')
    axes[3].legend(loc='lower center', bbox_to_anchor=(0.5, -0.12), frameon=False, fontsize=8.5)
    fig.suptitle('RNA: a ribose tem OH em 2\' e a uracila não tem a metila da timina; a fita única dobra sobre si mesma',
                 fontsize=11, y=1.04)
    done(fig, 'bio_rna.png')


# ---------------------------------------------------------------- 7. lipidios
def _rbc(ax, x, y, spur, r=0.95):
    """Hemacia vista de cima: disco com centro mais claro, ou com projecoes espinhosas (spur cell)."""
    col = '#c44e52'
    if spur:
        th = np.linspace(0, 2 * np.pi, 400)
        rr = r * (0.82 + 0.22 * np.clip(np.cos(9 * th + 0.4 * np.sin(3 * th)), 0, None) ** 6)
        ax.fill(x + rr * np.cos(th), y + rr * np.sin(th), color=col, alpha=0.85, ec='#8f2f33', lw=1, zorder=3)
    else:
        ax.add_patch(Circle((x, y), r, fc=col, alpha=0.85, ec='#8f2f33', lw=1, zorder=3))
        ax.add_patch(Circle((x, y), r * 0.45, fc='#e8a3a5', ec='none', zorder=4))


def _bilayer(ax, x0, n, pattern):
    """Bicamada de PC e colesterol; pattern: sequencia de 'L' (fosfolipidio) e 'C' (colesterol)."""
    for i in range(n):
        x = x0 + 0.5 + i * 0.62
        for up in (True, False):
            kind = pattern[(i + (0 if up else 1)) % len(pattern)]
            if kind == 'C':
                _chol(ax, x, 1.55 if up else -1.55, up, tlen=0.85)
            else:
                _lip(ax, x, 1.55 if up else -1.55, up, LIPID['PC'][0], tail='mono', tlen=1.1, r=0.26)


def fig_lipids():
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.4), gridspec_kw=dict(width_ratios=[1.2, 0.8, 1.6]))
    atoms, bonds = read_ccd('PC1')
    X = straight_tails(atoms); P = project(X, upright(atoms, X, 'P', ('C2I', 'C3I')))
    draw_molecule(axes[0], atoms, bonds, P, 'ballstick', scale=0.85)
    cx = P[:, 0].mean()
    groups = [('colina (+)', 'N', 1.5), ('fosfato (−)', 'P', 1.5), ('glicerol', 'C2', 0.0),
              ('cauda sn-1, 18 C', 'C3I', -1.0), ('cauda sn-2, 18 C', 'C2I', -1.0)]
    for lab, name, dy in groups:
        k = next(k for k, a in enumerate(atoms) if a['name'] == name)
        side = 1 if P[k, 0] >= cx else -1
        x = P[:, 0].max() + 4.5 if side > 0 else P[:, 0].min() - 4.5
        tag(axes[0], P[k, :2], lab, (x, P[k, 1] + dy), color=INK)
    fit(axes[0], P, pad=7.5)
    axes[0].set_title('DSPC: cabeça polar e duas caudas saturadas\n(caudas alinhadas por rotação rígida)', fontsize=10.5)
    element_legend(axes[0], ['H', 'C', 'N', 'O', 'P'], anchor=(0.5, -0.08))
    atoms, bonds = read_ccd('CLR')
    X = np.array([a['xyz'] for a in atoms]); P = project(X, upright(atoms, X, 'O1', ('C26', 'C27')))
    draw_molecule(axes[1], atoms, bonds, P, 'ballstick', scale=0.9)
    k = next(k for k, a in enumerate(atoms) if a['name'] == 'O1')
    tag(axes[1], P[k, :2], 'OH: a única parte polar', P[k, :2] + np.array([3.0, 1.2]), color=C[3])
    fit(axes[1], P, pad=1.5)
    axes[1].set_title('Colesterol: quatro anéis rígidos\nque ajustam a fluidez da membrana', fontsize=10)
    ax = axes[2]
    frame_off(ax); ax.set_aspect('equal')
    n = 9
    for x0, spur, pat, title, notes in (
            (0.0, False, 'LLC', 'Saudável: hemácia normal', ['colesterol e fosfolipídios', 'em proporção normal']),
            (7.6, True, 'LCC', 'Patológico: spur cell (doença hepática)', ['colesterol em excesso: membrana', 'rígida, célula com espículas,',
                                                                'destruída no baço'])):
        _bilayer(ax, x0, n, pat)
        w = n * 0.62 + 0.4
        for yb in (1.95, -1.95 - 0.55):
            ax.add_patch(Rectangle((x0, yb), w, 0.55, fc='#a6cfe3', alpha=0.45, lw=0, zorder=0))
        xc = x0 + (n * 0.62 + 0.4) / 2
        _rbc(ax, xc, 4.05, spur)
        ax.text(xc, 5.35, title, ha='center', va='bottom', fontsize=9.5, fontweight='bold')
        for k, s in enumerate(notes):
            ax.text(xc, -2.95 - 0.5 * k, s, ha='center', va='center', fontsize=8.5)
    ax.set_xlim(-0.3, 13.9); ax.set_ylim(-4.6, 6.2)
    hs = [plt.Line2D([], [], marker='o', ls='', ms=9, mfc=LIPID['PC'][0], mec=INK, label='fosfatidilcolina'),
          plt.Line2D([], [], marker='s', ls='', ms=9, mfc=LIPID['CHOL'][0], mec='none', label='colesterol')]
    ax.legend(handles=hs, loc='upper center', bbox_to_anchor=(0.5, -0.02), ncol=2, frameon=False, fontsize=8.5)
    ax.set_title('Colesterol na membrana da hemácia\n(Cooper, 1969; esquema)', fontsize=10)
    fig.suptitle('Lipídios: moléculas anfipáticas (DSPC e colesterol, coordenadas ideais do CCD)', fontsize=11, y=1.03)
    done(fig, 'bio_lipidios.png')


# ---------------------------------------------------------------- 8. carboidratos
SNFG = {'Glc': ('circle', '#0072bc'), 'Gal': ('circle', '#ffd400'), 'GalNAc': ('square', '#ffd400'),
        'GlcNAc': ('square', '#0072bc'), 'Fuc': ('triangle', '#ed1c24'), 'Man': ('circle', '#00a651')}


def snfg(ax, x, y, kind, s=0.42, z=5):
    shape, col = SNFG[kind]
    if shape == 'circle':
        ax.add_patch(Circle((x, y), s, fc=col, ec=INK, lw=0.8, zorder=z))
    elif shape == 'square':
        ax.add_patch(Rectangle((x - s, y - s), 2 * s, 2 * s, fc=col, ec=INK, lw=0.8, zorder=z))
    else:
        ax.add_patch(RegularPolygon((x, y), 3, radius=s * 1.25, fc=col, ec=INK, lw=0.8, zorder=z))


def fig_carbs():
    fig, axes = plt.subplots(1, 4, figsize=(16, 4.9), gridspec_kw=dict(width_ratios=[0.9, 1, 1.1, 1.2]))
    atoms, bonds = read_ccd('GLC')
    X = np.array([a['xyz'] for a in atoms])
    ring = [k for k, a in enumerate(atoms) if a['name'] in ('C1', 'C2', 'C3', 'C4', 'C5', 'O5')]
    nrm = frame(X[ring])[1][2]
    P = project(X, best_view(atoms, n=6000, normal=nrm, max_tilt=45))      # anel quase de frente
    draw_molecule(axes[0], atoms, bonds, P, 'ballstick')
    ctr = P[:, :2].mean(0)
    for k, a in enumerate(atoms):
        if a['elem'] == 'C':
            tag(axes[0], P[k, :2], a['name'], P[k, :2] + (P[k, :2] - ctr), box=False, fs=8.5)
    fit(axes[0], P, pad=2.6)
    axes[0].set_title('α-D-glicose: anel de 5 C e 1 O', fontsize=10)
    element_legend(axes[0], ['H', 'C', 'O'], anchor=(0.5, -0.1))
    # skeletal da maltose
    mal = [a for a in read_pdb('1ANF', chains='B', hetero=True) if a['resn'] == 'GLC']
    X = np.array([a['xyz'] for a in mal]); P = project(X, frame(X))
    draw_skeletal(axes[1], mal, dist_bonds(mal, 1.7), P)
    o4 = next(k for k, a in enumerate(mal) if a['name'] == 'O4' and a['resi'] == 1)
    tag(axes[1], P[o4, :2], 'α(1-4)', P[o4, :2] + np.array([0, -3.0]), color=C[3])
    fit(axes[1], P, pad=1.8)
    axes[1].set_title('Maltose (1ANF), fórmula de esqueleto:\nduas glicoses e uma ligação glicosídica', fontsize=10)
    # ABO em SNFG
    ax = axes[2]
    frame_off(ax); ax.set_aspect('equal'); ax.set_xlim(-0.8, 9.4); ax.set_ylim(-3.7, 3.6)
    for k, (name, extra) in enumerate((('O (antígeno H)', None), ('A', 'GalNAc'), ('B', 'Gal'))):
        x = 1 + 3.3 * k
        ax.plot([x, x], [-1.6, 1.0], color=INK, lw=1.0, zorder=1)
        snfg(ax, x, -1.6, 'GlcNAc'); snfg(ax, x, -0.3, 'Gal')
        ax.plot([x, x - 0.9], [-0.3, 0.6], color=INK, lw=1.0, zorder=1); snfg(ax, x - 0.9, 0.6, 'Fuc')
        if extra:
            snfg(ax, x, 1.0, extra)
        ax.text(x, -2.35, name, ha='center', va='top', fontsize=9.5, fontweight='bold')
        ax.plot([x - 1.0, x + 1.0], [-2.1, -2.1], color='#c44e52', lw=3)
    ax.text(4.3, -3.15, 'presos a lipídios e proteínas da membrana da hemácia (linha vermelha)', ha='center',
            va='top', fontsize=8, color='#5a5a5a')
    hs = [plt.Line2D([], [], marker={'circle': 'o', 'square': 's', 'triangle': '^'}[SNFG[k][0]], ls='', ms=9,
                     mfc=SNFG[k][1], mec=INK, label=k) for k in ('GlcNAc', 'Gal', 'Fuc', 'GalNAc')]
    ax.legend(handles=hs, loc='upper center', bbox_to_anchor=(0.5, 1.04), ncol=4, frameon=False, fontsize=8.5,
              handletextpad=0.2, columnspacing=0.8)
    ax.set_title('Saudável: grupos sanguíneos ABO\n(símbolos SNFG, esquema)', fontsize=10, pad=16)
    # escudo de glicanos
    ax = axes[3]
    frame_off(ax); ax.set_aspect('equal'); ax.set_xlim(-4.6, 5.9); ax.set_ylim(-3.6, 4.4)
    ax.add_patch(Rectangle((-4.4, -3.4), 8.8, 0.6, fc='#e6d9a8', ec='#b9a86a', lw=0.8))
    ax.add_patch(Polygon([(-1.2, -2.8), (1.2, -2.8), (2.2, 2.6), (0, 3.6), (-2.2, 2.6)], closed=True,
                         fc=C[0], alpha=0.25, ec=C[0], lw=1.4))
    rng = np.random.default_rng(3)
    for x0, y0, sgn in ((-1.4, -1.8, -1), (1.4, -1.8, 1), (-1.9, 0.2, -1), (1.9, 0.2, 1), (-2.0, 2.0, -1),
                        (2.0, 2.0, 1), (0.0, 3.6, 0), (-0.5, 0.8, -1), (0.6, -0.4, 1)):
        dx = 0.55 * sgn if sgn else 0
        dy = 0.0 if sgn else 0.55
        snfg(ax, x0 + dx, y0 + dy, 'GlcNAc', s=0.2)
        snfg(ax, x0 + 2 * dx, y0 + 2 * dy, 'GlcNAc', s=0.2)
        snfg(ax, x0 + 3 * dx, y0 + 3 * dy + (0.3 if sgn else 0), 'Man', s=0.2)
        snfg(ax, x0 + 3 * dx, y0 + 3 * dy - (0.3 if sgn else 0) + (0 if sgn else 0.5), 'Man', s=0.2)
        ax.plot([x0, x0 + 3 * dx], [y0, y0 + 3 * dy], color=INK, lw=0.7, zorder=1)
    # anticorpo em Y, com os bracos voltados para a proteina
    ab = C[3]
    ax.plot([5.5, 4.75], [1.2, 1.2], color=ab, lw=3, solid_capstyle='round')
    ax.plot([4.75, 4.15], [1.2, 1.75], color=ab, lw=3, solid_capstyle='round')
    ax.plot([4.75, 4.15], [1.2, 0.65], color=ab, lw=3, solid_capstyle='round')
    ax.text(4.85, 0.2, 'anticorpo\nnão alcança\na proteína', ha='center', va='top', fontsize=8.5, color=ab)
    ax.text(-4.4, -3.55, 'membrana do vírus', fontsize=8, color='#5a5a5a', va='top')
    hs = [plt.Line2D([], [], marker=m, ls='', ms=8, mfc=SNFG[k][1], mec=INK, label=k) for k, m in (('GlcNAc', 's'), ('Man', 'o'))]
    ax.legend(handles=hs, loc='upper left', bbox_to_anchor=(-0.02, 1.04), ncol=1, frameon=False, fontsize=8.5)
    ax.set_title('Patológico: glicanos cobrem a espícula\ndo SARS-CoV-2 (esquema)', fontsize=10, pad=16)
    fig.suptitle('Carboidratos: anéis de açúcar ligados por ligações glicosídicas, sozinhos ou presos a proteínas e lipídios',
                 fontsize=11, y=1.05)
    done(fig, 'bio_carboidratos.png')


# ---------------------------------------------------------------- 9. representacoes de uma molecula pequena
def fig_repr_small():
    atoms, bonds = read_ccd('PHE')
    X = np.array([a['xyz'] for a in atoms]); P = project(X, best_view(atoms))
    fig, axes = plt.subplots(1, 4, figsize=(16, 4.4))
    draw_skeletal(axes[0], atoms, bonds, P)
    draw_molecule(axes[1], atoms, bonds, P, 'sticks')
    draw_molecule(axes[2], atoms, bonds, P, 'ballstick')
    draw_molecule(axes[3], atoms, bonds, P, 'spacefill')
    titles = ['Fórmula de esqueleto: carbonos nos\nvértices, H só em heteroátomos',
              'Bastões: cada ligação é um cilindro\nna cor dos dois átomos',
              'Esferas e bastões: átomos menores\nque o tamanho real, ligações visíveis',
              'Preenchimento espacial: esferas com\no raio de van der Waals (Bondi, 1964)']
    for ax, t in zip(axes, titles):
        fit(ax, P, pad=2.0); ax.set_title(t, fontsize=10)
    element_legend(axes[1], ['H', 'C', 'N', 'O'], anchor=(1.1, -0.08))
    fig.suptitle('A mesma fenilalanina (coordenadas ideais do CCD, com hidrogênios) em quatro representações',
                 fontsize=11, y=1.04)
    done(fig, 'bio_representacoes_atomo.png')


# ---------------------------------------------------------------- 11. escalas de tempo dos movimentos
def fig_timescales():
    rows = [('vibração de ligação X-H', 1.1e-14, None, C[7]),
            ('rotação de cadeias laterais', 1e-12, 1e-9, C[0]),
            ('movimento de alças', 1e-9, 1e-6, C[0]),
            ('movimento de domínios e enovelamento', 1e-6, 1e0, C[0]),
            ('bicamada de DSPC se forma (Parte B)', 4e-9, 6.3e-9, C[2]),
            ('lipídio troca de monocamada (flip-flop)', 6.5 * 3600, None, C[5]),
            ('simulações desta aula: 20 fs a 15 ns', 20e-15, 15e-9, C[1])]
    fig, ax = plt.subplots(figsize=(13, 5.0))
    n = len(rows)
    for k, (lab, t0, t1, col) in enumerate(rows):
        y = n - 1 - k
        if t1 is None:
            ax.scatter(t0, y, s=120, marker='D', color=col, edgecolor=INK, lw=0.6, zorder=3)
        else:
            ax.plot([t0, t1], [y, y], color=col, lw=11, solid_capstyle='butt', zorder=3)
    ax.axhline(0.5, color=INK, lw=0.6)
    ax.set_yticks(range(n)); ax.set_yticklabels([r[0] for r in rows][::-1], fontsize=9.5)
    ax.tick_params(axis='y', length=0)
    ticks = [1e-15, 1e-12, 1e-9, 1e-6, 1e-3, 1, 3600]
    for x in ticks:
        ax.axvline(x, color=BAND, lw=0.8, zorder=0)
    ax.set_xscale('log'); ax.set_xlim(5e-16, 1e5); ax.set_ylim(-0.7, n - 0.4)
    ax.set_xticks(ticks); ax.set_xticklabels(['1 fs', '1 ps', '1 ns', '1 µs', '1 ms', '1 s', '1 h'])
    ax.minorticks_off()
    ax.tick_params(axis='x', direction='out', top=False)
    for sp in ('top', 'right', 'left'):
        ax.spines[sp].set_visible(False)
    ax.set_xlabel('tempo (escala logarítmica)')
    ax.text(11e-15 * 4, n - 1, '~11 fs', va='center', fontsize=8.5)
    ax.text(6.5 * 3600 * 4, n - 6, 'meia-vida de 6,5 h', va='center', fontsize=8.5)
    ax.set_title('Em que escala de tempo cada biomolécula se move (proteínas: Henzler-Wildman e Kern, 2007; '
                 'flip-flop: Kornberg e McConnell, 1971)', fontsize=11, pad=12)
    done(fig, 'bio_escalas_tempo.png')


def main():
    global CACHE
    ap = argparse.ArgumentParser(description='Gera imagens/bio_*.png')
    ap.add_argument('--cache', help='pasta para os arquivos baixados do PDB e do CCD (padrao: temporaria)')
    a = ap.parse_args()
    CACHE = Path(a.cache) if a.cache else Path(tempfile.mkdtemp())
    CACHE.mkdir(parents=True, exist_ok=True)
    for code in PDB_IDS:
        fetch(f'{code}.pdb', f'https://files.rcsb.org/download/{code}.pdb')
    for code in CCD_IDS:
        fetch(f'{code}.cif', f'https://files.rcsb.org/ligands/view/{code}.cif')
    fig_atoms(); fig_classes(); fig_protein_levels(); fig_collagen(); fig_dna(); fig_rna()
    fig_lipids(); fig_carbs(); fig_repr_small(); fig_timescales()


if __name__ == '__main__':
    main()
