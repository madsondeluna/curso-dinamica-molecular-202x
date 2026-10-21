"""Visualizacao das replicas Martini no estilo de tools/vmd_view.sh do AMP-evo-MD.

Cores com significado quimico (paleta seaborn "deep"), ligacoes CG desenhadas a
partir da topologia (sem a rede elastica), agua como camada translucida.

Exporta arquivos PDB com CONECT para o py3Dmol e monta tres visualizadores
lado a lado, sincronizados.
"""
from __future__ import annotations

import os
import re
from pathlib import Path

import numpy as np

# Cores do vmd_view.sh (seaborn "deep"): a mesma convencao no notebook inteiro
SEM = {
    'tail': '#b8b8b8', 'head_cat': '#4c72b0', 'po4': '#c44f52', 'glyc': '#54a869',
    'protein': '#de8552', 'protein_sc': '#eab08a', 'na': '#63b5cc', 'cl': '#d98cc2',
    'water': '#a6cfe3',
}


# ---------------------------------------------------------------- topologia
def itp_bonds(itp, molname=None):
    """Pares de ligacao (base 0) de um moleculetype: [bonds] e [constraints], sem rede elastica."""
    lines = open(itp).read().splitlines()
    pairs, inside, section, skip = set(), molname is None, None, False
    for i, raw in enumerate(lines):
        l = raw.strip()
        if l.startswith('[') and 'moleculetype' in l:
            # o nome vem na primeira linha nao comentada depois do cabecalho
            nxt = next((x.split()[0] for x in lines[i + 1:] if x.strip() and not x.strip().startswith(';')), '')
            inside = molname is None or nxt == molname
            section = None
            continue
        if not inside:
            continue
        if l.startswith('['):
            section = l.strip('[] ').lower(); skip = False; continue
        if l.startswith(';'):
            low = l.lower()
            if section in ('bonds', 'constraints'):
                skip = ('elastic' in low) or ('rubber' in low)
            continue
        if section in ('bonds', 'constraints') and not skip and l and not l.startswith('#'):
            parts = l.split()
            if len(parts) >= 2 and parts[0].isdigit() and parts[1].isdigit() and 'RUBBER' not in l:
                a, b = int(parts[0]) - 1, int(parts[1]) - 1
                pairs.add((min(a, b), max(a, b)))
    return sorted(pairs)


def system_bonds(topology, per_residue):
    """Expande ligacoes por tipo de residuo/molecula para o sistema inteiro.

    per_residue: {resname: [(i, j), ...]} com indices locais da molecula;
    'PROTEIN' usa indices globais da cadeia (comecando no primeiro atomo da proteina).
    """
    bonds = []
    prot = [a.index for a in topology.atoms if a.residue.is_protein]
    if 'PROTEIN' in per_residue and prot:
        off = prot[0]
        bonds += [(off + i, off + j) for i, j in per_residue['PROTEIN']]
    for res in topology.residues:
        if res.name in per_residue:
            idx = [a.index for a in res.atoms]
            bonds += [(idx[i], idx[j]) for i, j in per_residue[res.name] if j < len(idx)]
    return bonds


# ---------------------------------------------------------------- geometria
def wrap_molecules(xyz, box, topology):
    """Coloca o centro de cada molecula dentro da caixa, sem quebrar a molecula."""
    out = xyz.copy()
    for res in topology.residues:
        idx = [a.index for a in res.atoms]
        com = xyz[idx].mean(0)
        out[idx] -= np.floor(com / box) * box
    return out


def center_along(xyz, box, sel, axis):
    """Desloca o sistema para o centro (media circular de `sel`) ficar no meio da caixa em `axis`."""
    L = box[axis]
    ang = 2 * np.pi * xyz[sel, axis] / L
    c = np.arctan2(np.sin(ang).mean(), np.cos(ang).mean()) / (2 * np.pi) * L
    out = xyz.copy()
    out[:, axis] += L / 2 - c
    return out


def side_view(xyz, normal_axis):
    """Rotacao propria que leva a normal da bicamada para o eixo y (vertical na tela)."""
    if normal_axis == 1:
        return xyz
    if normal_axis == 2:   # (x, y, z) -> (x, z, -y)
        return np.stack([xyz[:, 0], xyz[:, 2], -xyz[:, 1]], axis=1)
    return np.stack([-xyz[:, 1], xyz[:, 0], xyz[:, 2]], axis=1)   # x -> y


# ---------------------------------------------------------------- PDB com CONECT
def write_pdb(path, topology, frames_nm, bonds, keep=None):
    """PDB multi-modelo (coordenadas em nm -> A) com CONECT das ligacoes CG."""
    keep = np.arange(topology.n_atoms) if keep is None else np.asarray(keep)
    newid = {int(o): k + 1 for k, o in enumerate(keep)}
    atoms = [topology.atom(int(i)) for i in keep]
    nb = {}
    for i, j in bonds:
        if i in newid and j in newid:
            nb.setdefault(newid[i], []).append(newid[j])
    conect = []
    for i in sorted(nb):
        partners = nb[i]
        for k in range(0, len(partners), 4):
            conect.append('CONECT' + f'{i:5d}' + ''.join(f'{j:5d}' for j in partners[k:k + 4]))
    out = []
    for m, xyz in enumerate(frames_nm, 1):
        out.append(f'MODEL     {m:4d}')
        for k, (a, p) in enumerate(zip(atoms, xyz[keep] * 10), 1):
            # HETATM: o 3Dmol so aplica CONECT a HETATM; para ATOM ele recalcula ligacoes por distancia
            out.append(f'HETATM{k % 100000:5d} {a.name[:4]:<4s} {a.residue.name[:4]:<4s}A{a.residue.resSeq % 10000:4d}    '
                       f'{p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}  1.00  0.00           C')
        out += conect               # CONECT dentro de cada modelo: o 3Dmol associa ao modelo corrente
        out.append('ENDMDL')
    out.append('END')
    Path(path).write_text('\n'.join(out) + '\n')


# ---------------------------------------------------------------- exportacao por replica
def export_replica(rep_dir, n_frames_prot=31, n_frames_mem=16):
    """Gera viz_prot.pdb, viz_mem_final.pdb e viz_mem_anim.pdb dentro de rep_dir."""
    import mdtraj as md
    d = Path(rep_dir)
    prot_bonds = itp_bonds(d / 'molecule_0.itp')
    lip_bonds = itp_bonds(d / 'martini_v3.0.0_phospholipids_PC_v2.itp', 'DSPC')
    out = {}
    # proteina: trajetoria ja ajustada (so a proteina)
    if (d / 'md_prot_fit.xtc').exists() and not (d / 'viz_prot.pdb').exists():
        t = md.load(str(d / 'md_prot_fit.xtc'), top=str(d / 'ubq_cg.pdb'))
        idx = np.linspace(0, t.n_frames - 1, n_frames_prot).astype(int)
        xyz = t.xyz[idx] - t.xyz[idx].mean(1, keepdims=True)
        write_pdb(d / 'viz_prot.pdb', t.topology, xyz, [(i, j) for i, j in prot_bonds])
    out['prot'] = d / 'viz_prot.pdb'
    # membrana: automontagem (so lipidios) e estado final (com agua)
    if (d / 'mem_whole.xtc').exists():
        t = md.load(str(d / 'mem_whole.xtc'), top=str(d / 'mem_em.gro'))
        top = t.topology
        bonds = system_bonds(top, {'DSPC': lip_bonds})
        lip = top.select('resname DSPC')
        if not (d / 'viz_mem_anim.pdb').exists():
            idx = np.linspace(0, t.n_frames - 1, n_frames_mem).astype(int)
            frames = [wrap_molecules(t.xyz[i], t.unitcell_lengths[i], top) for i in idx]
            write_pdb(d / 'viz_mem_anim.pdb', top, frames, bonds, keep=lip)
        if not (d / 'viz_mem_final.pdb').exists():
            src = d / 'mem_eq_whole.xtc'
            if src.exists():   # bicamada formada: estado equilibrado, normal em z
                te = md.load(str(src), top=str(d / 'mem_rot.gro'))
                xyz, box, axis = te.xyz[-1], te.unitcell_lengths[-1], 2
            else:              # sem bicamada: ultimo quadro da automontagem
                xyz, box, axis = t.xyz[-1], t.unitcell_lengths[-1], None
            tails = top.select('name =~ "C[1-4][AB]"')
            if axis is not None:
                xyz = center_along(xyz, box, tails, axis)
            xyz = wrap_molecules(xyz, box, top)
            if axis is not None:
                xyz = side_view(xyz, axis)
            water = top.select('resname W')
            write_pdb(d / 'viz_mem_final.pdb', top, [xyz], bonds, keep=top.select('not resname W'))
            write_pdb(d / 'viz_mem_water.pdb', top, [xyz], [], keep=water)
    out['mem_anim'] = d / 'viz_mem_anim.pdb'
    out['mem_final'] = d / 'viz_mem_final.pdb'
    out['mem_water'] = d / 'viz_mem_water.pdb'
    return out


# ---------------------------------------------------------------- estilos py3Dmol
_TAILS = [f'C{i}{s}' for i in range(1, 5) for s in 'AB']


def style_protein(view, viewer):
    """Bolas e varetas: bastao nas ligacoes CG, esfera menor que a metade da ligacao."""
    view.setStyle({}, {'stick': {'color': SEM['protein'], 'radius': 0.5},
                       'sphere': {'color': SEM['protein_sc'], 'radius': 1.1}}, viewer=viewer)
    view.setStyle({'atom': 'BB'}, {'stick': {'color': SEM['protein'], 'radius': 0.5},
                                   'sphere': {'color': SEM['protein'], 'radius': 1.4}}, viewer=viewer)


_LIPID = ['NC3', 'PO4', 'GL1', 'GL2'] + _TAILS


def style_membrane(view, viewer, water=True):
    """Selecao por nome de bead: o 3Dmol le so 3 caracteres do nome do residuo (DSPC vira DSP)."""
    view.setStyle({}, {}, viewer=viewer)
    view.setStyle({'atom': _TAILS}, {'stick': {'color': SEM['tail'], 'radius': 0.5},
                                     'sphere': {'color': SEM['tail'], 'radius': 0.9}}, viewer=viewer)
    for names, color, r in ((['NC3'], SEM['head_cat'], 1.7), (['PO4'], SEM['po4'], 1.6),
                            (['GL1', 'GL2'], SEM['glyc'], 1.4)):
        view.setStyle({'atom': names}, {'stick': {'color': SEM['tail'], 'radius': 0.5},
                                        'sphere': {'color': color, 'radius': r}}, viewer=viewer)
    view.setStyle({'atom': 'NA'}, {'sphere': {'color': SEM['na'], 'radius': 1.6}}, viewer=viewer)
    view.setStyle({'atom': 'CL'}, {'sphere': {'color': SEM['cl'], 'radius': 1.6}}, viewer=viewer)


def triple_view(paths, kind, labels, height=420, animate=False, water=True):
    """Tres visualizadores lado a lado, camera compartilhada; animacoes em sincronia."""
    import py3Dmol
    view = py3Dmol.view(width='100%', height=height, viewergrid=(1, len(paths)), linked=True)
    for col, path in enumerate(paths):
        v = (0, col)
        if path is None or not Path(path).exists():
            continue
        text = Path(path).read_text()
        if animate:
            view.addModelsAsFrames(text, 'pdb', viewer=v)
        else:
            view.addModel(text, 'pdb', viewer=v)
        (style_protein if kind == 'protein' else lambda vw, viewer: style_membrane(vw, viewer, water))(view, v)
        wfile = Path(str(path).replace('mem_final', 'mem_water'))   # viz_mem_final.pdb ou repN_mem_final.pdb
        if kind == 'membrane' and water and not animate and wfile.exists():
            # agua em modelo proprio: no 3Dmol a opacidade vale por modelo
            view.addModel(wfile.read_text(), 'pdb', viewer=v)
            view.setStyle({'atom': 'W'}, {'sphere': {'color': SEM['water'], 'radius': 2.0, 'opacity': 0.35}}, viewer=v)
        view.setBackgroundColor('white', 0, viewer=v)
    view.zoomTo()
    if animate:
        view.animate({'loop': 'forward', 'reps': 0, 'interval': 180})
    head = ''.join(f'<li class="compare-ff">{l}</li>' for l in labels)
    return f'<ul class="compare-labels">{head}</ul>{view._make_html()}'
