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

# Barra sob cada visualizador: escolher o quadro e baixar imagem (PNG) ou estrutura (PDB) do quadro.
# Funciona com o py3Dmol, que declara o visualizador como variavel global viewer_<id> no HTML.
_CONTROLS = """
<div id="ctl_UID" style="font:13px sans-serif;color:#111827;padding:6px 8px;display:flex;flex-wrap:wrap;gap:8px;align-items:center">
  <span id="ctl_UID_anim" style="display:none;align-items:center;gap:8px">
    <button type="button" id="ctl_UID_play">Pausar</button>
    <input type="range" id="ctl_UID_frame" min="0" max="0" value="0" style="width:220px">
    <span id="ctl_UID_label"></span>
  </span>
  <button type="button" id="ctl_UID_png">Baixar imagem (PNG)</button>
  <button type="button" id="ctl_UID_pdb">Baixar estrutura do quadro (PDB)</button>
</div>
<script>
(function () {
  var name = 'NAME';
  function v() { return window['viewer_UID']; }
  function el(s) { return document.getElementById('ctl_UID_' + s); }
  function save(href, file) {
    var a = document.createElement('a'); a.href = href; a.download = file;
    document.body.appendChild(a); a.click(); a.remove();
  }
  var playing = true;
  function label() {
    var n = v().getNumFrames(), f = v().getFrame();
    el('label').textContent = 'quadro ' + (f + 1) + ' de ' + n;
    el('frame').value = f;
  }
  var wait = setInterval(function () {
    if (!v()) return;
    clearInterval(wait);
    var n = v().getNumFrames();
    if (n > 1) {
      el('anim').style.display = 'inline-flex';
      el('frame').max = n - 1;
      setInterval(function () { if (playing) label(); }, 250);
      label();
    }
  }, 200);
  el('play').onclick = function () {
    if (playing) { v().pauseAnimate(); this.textContent = 'Continuar'; }
    else { v().resumeAnimate(); this.textContent = 'Pausar'; }
    playing = !playing;
  };
  el('frame').oninput = function () {
    if (playing) { v().pauseAnimate(); el('play').textContent = 'Continuar'; playing = false; }
    var f = parseInt(this.value);
    Promise.resolve(v().setFrame(f)).then(function () { v().render(); label(); });
  };
  function tag() { return v().getNumFrames() > 1 ? '_quadro' + (v().getFrame() + 1) : ''; }
  el('png').onclick = function () { v().render(); save(v().pngURI(), name + tag() + '.png'); };
  el('pdb').onclick = function () {
    var atoms = v().selectedAtoms({}), out = [], s = 0;
    function pad(x, n, left) { x = String(x); while (x.length < n) x = left ? x + ' ' : ' ' + x; return x.slice(0, n); }
    atoms.forEach(function (a) {
      s += 1;
      var nm = a.atom || a.elem || 'X';
      out.push('HETATM' + pad(s % 100000, 5) + ' ' + pad(nm.length < 4 ? ' ' + nm : nm, 4, true) + ' ' +
               pad(a.resn || nm, 3) + ' ' + (a.chain || 'A').slice(0, 1) + pad((a.resi || 1) % 10000, 4) + '    ' +
               pad(a.x.toFixed(3), 8) + pad(a.y.toFixed(3), 8) + pad(a.z.toFixed(3), 8) + '  1.00  0.00');
    });
    out.push('END');
    var blob = new Blob([out.join('\\n') + '\\n'], {type: 'chemical/x-pdb'});
    save(URL.createObjectURL(blob), name + tag() + '.pdb');
  };
})();
</script>
"""


def _controls(inner, name='estrutura'):
    """Barra de quadro e downloads para o visualizador contido em `inner` (vazia se nao houver um)."""
    import re as _re
    m = _re.search(r'viewer_(\w+)\s*=\s*\$3Dmol\.createViewer', inner)
    if not m:
        return ''
    return _CONTROLS.replace('UID', m.group(1)).replace('NAME', name)


def white_box(inner, name='estrutura'):
    """Envolve o HTML do visualizador num bloco branco (o Colab em tema escuro pinta o fundo de preto)
    e acrescenta a barra de quadro e downloads."""
    return ('<div style="background:#ffffff;color:#111827;color-scheme:light;'
            f'display:inline-block;max-width:100%;border-radius:4px">{inner}{_controls(inner, name)}</div>')


def show_white(view, name='estrutura'):
    """Mostra o visualizador py3Dmol sempre sobre fundo branco opaco, com a barra de quadro e downloads."""
    from IPython.display import HTML, display
    view.setBackgroundColor('#ffffff', 1)
    display(HTML(white_box(view._make_html(), name)))

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


def center_protein(xyz, box, prot):
    """Proteina inteira no centro da caixa; agua e ions reembrulhados dentro dela."""
    out = xyz.copy()
    for k in range(1, len(prot)):            # desembrulha bead a bead ao longo da cadeia
        d = out[prot[k]] - out[prot[k - 1]]
        out[prot[k]] = out[prot[k - 1]] + d - box * np.round(d / box)
    out += box / 2 - out[prot].mean(axis=0)
    rest = np.setdiff1d(np.arange(len(out)), prot)
    out[rest] -= box * np.floor(out[rest] / box)
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


# ---------------------------------------------------------------- filme da automontagem
def _write_movie(d, t, top, bonds, lip, n_frames, name):
    """Filme com a normal final da bicamada na vertical e cada quadro centrado nas caudas."""
    tails_all = top.select('name =~ "C[1-4][AB]"')
    order = d / 'mem_order.npz'
    normal_axis = int(np.argmax(np.abs(np.load(order)['normal']))) if order.exists() else 2
    perm = {0: [1, 2, 0], 1: [2, 0, 1], 2: [0, 1, 2]}[normal_axis]
    idx = np.linspace(0, t.n_frames - 1, n_frames).astype(int)
    frames = []
    for i in idx:
        box = t.unitcell_lengths[i][perm]
        xyz = center_along(t.xyz[i][:, perm], box, tails_all, 2)
        frames.append(side_view(wrap_molecules(xyz, box, top), 2))
    write_pdb(d / f'{name}.pdb', top, frames, bonds, keep=lip)
    write_xyz(d / f'{name}_water.xyz', frames, top.select('resname W'))


def export_protein_movie(rep_dir, n_frames=31, name='viz_prot_movie'):
    """Filme completo da proteina com o sistema inteiro (agua em arquivo separado)."""
    import mdtraj as md
    d = Path(rep_dir)
    if not (d / f'{name}.pdb').exists():
        t = md.load(str(d / 'md_prot_pbc.xtc'), top=str(d / 'ubq_ions.gro'))
        top = t.topology
        idx = np.linspace(0, t.n_frames - 1, n_frames).astype(int)
        prot = top.select('not resname W and not name NA CL')
        frames = [center_protein(t.xyz[i], t.unitcell_lengths[i], prot) for i in idx]
        write_pdb(d / f'{name}.pdb', top, frames, itp_bonds(d / 'molecule_0.itp'), keep=top.select('not resname W'))
        write_xyz(d / f'{name}_water.xyz', frames, top.select('resname W'))
    return d / f'{name}.pdb', d / f'{name}_water.xyz'


def export_movie(rep_dir, n_frames=31, name='viz_mem_movie'):
    """Filme completo da automontagem (lipidios e agua em arquivos separados)."""
    import mdtraj as md
    d = Path(rep_dir)
    if not (d / f'{name}.pdb').exists():
        t = md.load(str(d / 'mem_whole.xtc'), top=str(d / 'mem_em.gro'))
        top = t.topology
        bonds = system_bonds(top, {'DSPC': itp_bonds(d / 'martini_v3.0.0_phospholipids_PC_v2.itp', 'DSPC')})
        _write_movie(d, t, top, bonds, top.select('resname DSPC'), n_frames, name)
    return d / f'{name}.pdb', d / f'{name}_water.xyz'


def write_xyz(path, frames_nm, keep, label='W'):
    """XYZ multi-quadro (A) para a agua: sem ligacoes, ~3 vezes menor que o PDB."""
    out = []
    for xyz in frames_nm:
        p = xyz[keep] * 10
        out.append(str(len(p)))
        out.append('agua W')
        out += [f'{label} {a:.1f} {b:.1f} {c:.1f}' for a, b, c in p]
    Path(path).write_text('\n'.join(out) + '\n')


# ---------------------------------------------------------------- exportacao por replica
def export_replica(rep_dir, n_frames_prot=16, n_frames_mem=16):
    """Gera viz_prot.pdb, viz_mem_final.pdb e viz_mem_anim.pdb dentro de rep_dir."""
    import mdtraj as md
    d = Path(rep_dir)
    prot_bonds = itp_bonds(d / 'molecule_0.itp')
    lip_bonds = itp_bonds(d / 'martini_v3.0.0_phospholipids_PC_v2.itp', 'DSPC')
    out = {}
    # proteina: sistema inteiro centrado na proteina (agua em arquivo separado, translucida)
    if (d / 'md_prot_pbc.xtc').exists() and not (d / 'viz_prot.pdb').exists():
        t = md.load(str(d / 'md_prot_pbc.xtc'), top=str(d / 'ubq_ions.gro'))
        top = t.topology
        idx = np.linspace(0, t.n_frames - 1, n_frames_prot).astype(int)
        prot = top.select('not resname W and not name NA CL')
        frames = [center_protein(t.xyz[i], t.unitcell_lengths[i], prot) for i in idx]
        write_pdb(d / 'viz_prot.pdb', top, frames, list(prot_bonds), keep=top.select('not resname W'))
        write_xyz(d / 'viz_prot_water.xyz', frames, top.select('resname W'))
    out['prot'] = d / 'viz_prot.pdb'
    out['prot_water'] = d / 'viz_prot_water.xyz'
    # membrana: automontagem (so lipidios) e estado final (com agua)
    if (d / 'mem_whole.xtc').exists():
        t = md.load(str(d / 'mem_whole.xtc'), top=str(d / 'mem_em.gro'))
        top = t.topology
        bonds = system_bonds(top, {'DSPC': lip_bonds})
        lip = top.select('resname DSPC')
        if not (d / 'viz_mem_anim.pdb').exists():
            _write_movie(d, t, top, bonds, lip, n_frames_mem, 'viz_mem_anim')
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
            write_xyz(d / 'viz_mem_water.xyz', [xyz], water)
    out['mem_anim'] = d / 'viz_mem_anim.pdb'
    out['mem_anim_water'] = d / 'viz_mem_anim_water.xyz'
    out['mem_final'] = d / 'viz_mem_final.pdb'
    out['mem_water'] = d / 'viz_mem_water.xyz'
    return out


# ---------------------------------------------------------------- estilos py3Dmol
_TAILS = [f'C{i}{s}' for i in range(1, 5) for s in 'AB']


def style_protein(view, viewer=None):
    """Bolas e varetas: bastao nas ligacoes CG, esfera menor que a metade da ligacao."""
    kw = {'viewer': viewer} if viewer is not None else {}
    view.setStyle({}, {'stick': {'color': SEM['protein'], 'radius': 0.5},
                       'sphere': {'color': SEM['protein_sc'], 'radius': 1.1}}, **kw)
    view.setStyle({'atom': 'BB'}, {'stick': {'color': SEM['protein'], 'radius': 0.5},
                                   'sphere': {'color': SEM['protein'], 'radius': 1.4}}, **kw)
    view.setStyle({'atom': 'NA'}, {'sphere': {'color': SEM['na'], 'radius': 1.6}}, **kw)
    view.setStyle({'atom': 'CL'}, {'sphere': {'color': SEM['cl'], 'radius': 1.6}}, **kw)




def style_membrane(view, viewer=None, water=True):
    """Selecao por nome de bead: o 3Dmol le so 3 caracteres do nome do residuo (DSPC vira DSP)."""
    kw = {'viewer': viewer} if viewer is not None else {}
    view.setStyle({}, {}, **kw)
    view.setStyle({'atom': _TAILS}, {'stick': {'color': SEM['tail'], 'radius': 0.5},
                                     'sphere': {'color': SEM['tail'], 'radius': 0.9}}, **kw)
    for names, color, r in ((['NC3'], SEM['head_cat'], 1.7), (['PO4'], SEM['po4'], 1.6),
                            (['GL1', 'GL2'], SEM['glyc'], 1.4)):
        view.setStyle({'atom': names}, {'stick': {'color': SEM['tail'], 'radius': 0.5},
                                        'sphere': {'color': color, 'radius': r}}, **kw)
    view.setStyle({'atom': 'NA'}, {'sphere': {'color': SEM['na'], 'radius': 1.6}}, **kw)
    view.setStyle({'atom': 'CL'}, {'sphere': {'color': SEM['cl'], 'radius': 1.6}}, **kw)


def single_view(path, kind, animate=False, water_path=None, height=420, width=820):
    """Um visualizador por estrutura (fundo branco); agua translucida em modelo separado."""
    import py3Dmol
    if path is None or not Path(path).exists():
        print('arquivo de visualizacao ausente:', path)
        return
    view = py3Dmol.view(width=width, height=height)
    text = Path(path).read_text()
    if animate:
        view.addModelsAsFrames(text, 'pdb')
    else:
        view.addModel(text, 'pdb')
    if kind == 'protein':
        style_protein(view, None)
    else:
        style_membrane(view, None, water=False)
    if water_path and Path(water_path).exists():
        wtext = Path(water_path).read_text()
        wfmt = 'xyz' if str(water_path).endswith('.xyz') else 'pdb'
        if animate:
            view.addModelsAsFrames(wtext, wfmt)         # agua animada junto, em modelo proprio
        else:
            view.addModel(wtext, wfmt)
        # agua como lamina continua (efeito QuickSurf): esferas grandes, sobrepostas e translucidas
        for sel in ({'atom': 'W'}, {'elem': 'W'}):
            view.setStyle(sel, {'sphere': {'color': '#5b9bd5', 'radius': 3.0, 'opacity': 0.3}})
    if animate:
        view.animate({'loop': 'forward', 'reps': 0, 'interval': 180})
    view.zoomTo()
    show_white(view, Path(path).stem)

