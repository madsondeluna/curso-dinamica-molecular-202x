"""Gera MD_Martini.ipynb a partir das celulas definidas abaixo."""
import json, sys

cells = []

def md(s):
    cells.append({"cell_type": "markdown", "metadata": {}, "source": s.strip("\n")})

def code(s):
    if not s.strip():
        return
    cells.append({"cell_type": "code", "metadata": {}, "execution_count": None,
                  "outputs": [], "source": s.strip("\n")})

md('''
<a href="https://colab.research.google.com/github/madsondeluna/curso-dinamica-molecular-202x/blob/main/MD_Martini.ipynb" target="_parent"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a>
''')

md('''
# Dinâmica Molecular Coarse-Grained com Martini 3

Este notebook simula dois sistemas de natureza diferente com o mesmo campo de força coarse-grained (CG):

| Parte | Sistema | Fenômeno | Tempo simulado |
|-------|---------|----------|----------------|
| A | Ubiquitina (1UBQ) em água + 0,15 M NaCl | Dinâmica de uma proteína enovelada | 0,5 ns (equilibração) + 15 ns |
| B | 128 DSPC + 1 536 W | Automontagem espontânea de uma bicamada lipídica | 15 ns |

A Parte A usa a mesma proteína do curso atomístico anterior (`legado/MD_GMX.ipynb`). A Parte B reproduz o tutorial oficial *Lipid bilayers, Part I: Self-assembly* do Martini 3.

Cada etapa termina com um **cartão de etapa**: visualização 3D interativa (arraste para girar, role para aproximar), a explicação do que aconteceu, os números do sistema naquele momento e a mesma etapa comparada em três resoluções:

| Resolução | Partícula | Exemplo de campo de força | Usado aqui |
|-----------|-----------|---------------------------|------------|
| All-atom (AA) | cada átomo, inclusive todos os H | CHARMM36m (`legado/MD_GMX.ipynb`), AMBER99SB-ILDN | referência |
| United-atom (UA) | átomos pesados; CH, CH₂ e CH₃ viram uma partícula, só os H polares ficam | GROMOS 54A7 | referência |
| Coarse-grained (CG) | bead de ~4 átomos pesados | Martini 3 | **caso desta aula** |
''')

md('''
| Autoria | |
|---|---|
| Autor | **Madson A. de Luna-Aragão**, MSc |
| Vínculos | Programa Interunidades de Pós-Graduação em Bioinformática, Universidade Federal de Minas Gerais (UFMG), Belo Horizonte, Brasil<br>Centre National de la Recherche Scientifique (CNRS), Institut de Biologie et Chimie des Protéines (IBCP), Lyon, França |
| E-mail | [madsondeluna@gmail.com](mailto:madsondeluna@gmail.com) |
| Site | [madsondeluna.com](https://madsondeluna.com) |
| GitHub | [github.com/madsondeluna](https://github.com/madsondeluna) |
| LinkedIn | [linkedin.com/in/madsonaragao](https://www.linkedin.com/in/madsonaragao/) |
| Repositório | [madsondeluna/curso-dinamica-molecular-202x](https://github.com/madsondeluna/curso-dinamica-molecular-202x) |
| Guia da aula | [README com todas as figuras](https://github.com/madsondeluna/curso-dinamica-molecular-202x#readme) |
''')

md('''
## Antes de começar: ambiente limpo e sessão ativa

Execute as três células abaixo **antes de qualquer outra coisa**. A primeira reinicia o ambiente sozinha se ele já tiver sido usado: a execução para com a mensagem *Ambiente reiniciado* (o Colab pode avisar que a sessão travou, é esperado). Nesse caso, rode **Ambiente de execução > Executar tudo** de novo. Os arquivos ficam em `/content/aula_martini/repN`. Só para simular do zero (`USE_PRECOMPUTED = False`) o notebook monta o Google Drive e grava em `MyDrive/aula_martini/repN`, para retomar do último checkpoint se a sessão cair.

**GPU.** Para simular do zero (`USE_PRECOMPUTED = False`), escolha **Ambiente de execução > Alterar o tipo de ambiente de execução > T4 GPU** antes de executar. O notebook detecta a GPU e a usa no `mdrun` quando o GROMACS instalado tem suporte; sem GPU, ou sem suporte, a simulação roda na CPU. Com a réplica pré-computada, a GPU não faz diferença.

| Parâmetro | Valor padrão | Efeito |
|-----------|--------------|--------|
| `REPLICA` | `1` | qual das três réplicas usar (1, 2 ou 3); cada uma tem sementes aleatórias diferentes |
| `USE_PRECOMPUTED` | `True` | baixa do GitHub a réplica já simulada: todas as etapas aparecem como "ja pronto" e os cartões e gráficos saem na hora. Com `False`, o notebook simula tudo no Colab (2 a 3 h) | Ao reabrir o notebook, cada etapa já concluída é detectada e reaproveitada; uma simulação interrompida continua do último checkpoint em vez de recomeçar.
''')

code('''
# Reinicio limpo: se esta celula ja rodou neste ambiente, reinicia o Python antes de recomecar.
# Variaveis e modulos de uma execucao anterior causam erros ao recomecar do inicio.
# Os arquivos em /content e no Drive continuam la; so a memoria do Python e limpa.
import os, sys, threading
if 'google.colab' in sys.modules and globals().get('_AULA_INICIADA'):
    # o erro abaixo interrompe o "Executar tudo"; 1 s depois o processo sai e o Colab sobe um ambiente novo
    threading.Timer(1.0, os._exit, args=(0,)).start()
    raise SystemExit('Ambiente reiniciado. Execute de novo: Ambiente de execucao > Executar tudo.')
_AULA_INICIADA = True
print('Ambiente limpo: pode seguir.')
''')

code('''
# Mantem a sessao ativa clicando automaticamente no botao de conexao
# Execute esta celula logo apos abrir o notebook
from IPython.display import Javascript, display
display(Javascript("""
function keepAlive() {
    const btn = document.querySelector('colab-toolbar-button#connect');
    if (btn) { btn.click(); console.log('Colab keep-alive ping'); }
}
setInterval(keepAlive, 60000);  // a cada 60 segundos
console.log('Keep-alive iniciado.');
"""))
''')

code('''
# Parametros da aula
REPLICA = 1              # 1, 2 ou 3: cada grupo usa uma replica diferente
USE_PRECOMPUTED = True   # True: baixa a replica ja simulada do GitHub; False: simula tudo do zero

# Diretorio de trabalho. O Google Drive so e montado para simular do zero:
# a corrida leva horas, e o Drive guarda checkpoints se a sessao cair.
# Com a replica pre-computada, tudo fica em /content e baixa de novo em segundos.
import os

BASE_DIR = '/content/aula_martini' if os.path.isdir('/content') else os.path.abspath('aula_martini')
if not USE_PRECOMPUTED:
    try:
        from google.colab import drive
        drive.mount('/content/drive')
        BASE_DIR = '/content/drive/MyDrive/aula_martini'
        print('Google Drive montado: a simulacao continua de onde parou se a sessao cair.')
    except Exception:
        print('Drive nao montado: arquivos nao persistem apos a sessao.')
WORK_DIR = os.path.join(BASE_DIR, f'rep{REPLICA}')
os.makedirs(WORK_DIR, exist_ok=True)
os.chdir(WORK_DIR)
print('Replica', REPLICA, '| diretorio de trabalho:', os.getcwd())
''')

md('''
# Instalação e Configuração do Ambiente

| Ferramenta | Função |
|------------|--------|
| GROMACS (apt) | Motor de simulação |
| `vermouth` (martinize2) | Converte a proteína atomística em modelo Martini 3 |
| `mdtraj` | DSSP para o martinize2 e análises em Python |
| `py3Dmol` | Visualização 3D |
| `viz/` deste repositório | Cartões de etapa e tema dos gráficos |

> **Ambiente:** todas as simulações rodam na CPU. O Colab gratuito (2 vCPUs) é suficiente.
''')

code('''
import os, sys, shutil, subprocess

def cached(files, label):
    """Retorna True se todos os arquivos existirem e nao estiverem vazios: etapa ja concluida."""
    missing = [f for f in files if not os.path.exists(f) or os.path.getsize(f) == 0]
    if not missing:
        print(f"  ja pronto: {label}")
        return True
    return False

# GROMACS (instalacao de sistema, necessaria a cada sessao do Colab)
if shutil.which('gmx') is None:
    os.system('apt-get update -q && apt-get install -y -q gromacs')
subprocess.run([sys.executable, '-m', 'pip', 'install', '-q',
                'vermouth', 'mdtraj', 'py3Dmol'], capture_output=True)
os.system("gmx --version 2>&1 | grep 'GROMACS version'")
os.system("martinize2 -V 2>&1 | tail -1")

# GPU: usada quando o Colab tem uma GPU e o GROMACS instalado tem suporte (CUDA ou OpenCL)
_gpu_hw = shutil.which('nvidia-smi') is not None and os.system('nvidia-smi -L > /dev/null 2>&1') == 0
_ver = subprocess.run(['gmx', '--version'], capture_output=True, text=True).stdout
_gpu_build = next((l.split(':', 1)[1].strip() for l in _ver.splitlines() if l.strip().startswith('GPU support')), 'disabled')
USE_GPU = _gpu_hw and _gpu_build.lower() != 'disabled'
print('GPU:', 'sim (' + _gpu_build + ')' if USE_GPU else
      ('detectada, mas o GROMACS instalado nao tem suporte: CPU' if _gpu_hw else 'nao detectada: CPU'))

# Uma unica simulacao por vez: 1 rank MPI, OpenMP usa todos os nucleos
MDRUN_CPU = 'gmx mdrun -ntmpi 1'
MDRUN = MDRUN_CPU + (' -nb gpu' if USE_GPU else '')
os.environ['GMX_MAXBACKUP'] = '-1'   # sobrescreve arquivos sem criar #backups#

import json, random, re

SEEDS_FILE = 'seeds.json'

def seeds():
    return json.load(open(SEEDS_FILE)) if os.path.exists(SEEDS_FILE) else {}

def record_seed(step, **values):
    """Guarda em seeds.json as sementes usadas em cada etapa (no diretorio da replica)."""
    s = seeds()
    s.setdefault(step, {}).update(values)
    json.dump(s, open(SEEDS_FILE, 'w'), indent=1)

def tool_seed(step):
    """Semente explicita para insert-molecules e genion: sorteada uma vez e reaproveitada."""
    s = seeds().get(step, {}).get('seed')
    if s is None:
        s = random.randint(1, 2**31 - 1)
        record_seed(step, seed=s)
    return s

def run_md(name, label, mdp, conf, top, extra=''):
    """grompp + mdrun com retomada.

    A etapa so conta como concluida quando <name>.gro existe: o mdrun escreve
    esse arquivo apenas no ultimo passo. Se a sessao cair no meio, a proxima
    execucao reaproveita o .tpr e continua do checkpoint <name>.cpt.
    As sementes sorteadas pelo grompp (gen_seed = -1, ld_seed = -1) so
    aparecem na saida dele; ficam registradas em seeds.json.
    """
    if cached([f'{name}.gro'], label):
        return
    if not os.path.exists(f'{name}.tpr'):
        os.system(f'gmx grompp -f {mdp} -c {conf} -p {top} {extra} -o {name}.tpr -quiet -maxwarn 0 '
                  f'> {name}_grompp.log 2>&1')
        out = open(f'{name}_grompp.log').read()
        if not os.path.exists(f'{name}.tpr'):
            print(out[-2000:])
            raise RuntimeError(f'grompp falhou em {label}')
        gen = re.search(r'Setting gen_seed to (-?\\d+)', out)
        ld = re.search(r'Setting the LD random seed to (-?\\d+)', out)
        record_seed(name, **({'gen_seed': int(gen.group(1))} if gen else {}),
                    **({'ld_seed': int(ld.group(1))} if ld else {}))
        print(f'  sementes de {label}:', seeds().get(name, {}))
    if os.path.exists(f'{name}.cpt'):
        print(f'  retomando {label} do checkpoint {name}.cpt')
    # -cpt 5: checkpoint a cada 5 min; -cpi sem arquivo existente comeca do zero
    status = os.system(f'{MDRUN} -deffnm {name} -cpi {name}.cpt -cpt 5 -quiet')
    if status != 0 and MDRUN != MDRUN_CPU and not os.path.exists(f'{name}.gro'):
        print(f'  {label}: a GPU falhou; repetindo na CPU')
        os.system(f'{MDRUN_CPU} -deffnm {name} -cpi {name}.cpt -cpt 5 -quiet')
    print(f'  {label}: ' + ('concluido' if os.path.exists(f'{name}.gro') else 'INCOMPLETO, execute a celula de novo'))

def gmx_cached(outputs, label, cmd):
    """Roda um comando gmx de analise apenas se as saidas ainda nao existirem."""
    if not cached(outputs, label):
        os.system(cmd + ' -quiet > /dev/null 2>&1')

def performance(log):
    """ns/dia registrado no fim de um .log do mdrun (None se a corrida nao terminou)."""
    if not os.path.exists(log):
        return None
    for line in open(log):
        if line.startswith('Performance:'):
            return float(line.split()[1])
    return None

def save_summary(**parts):
    """Acrescenta resultados desta replica a summary.json (usado na comparacao entre replicas)."""
    out = json.load(open('summary.json')) if os.path.exists('summary.json') else {}
    out.update(parts)
    out['replica'] = REPLICA
    json.dump(out, open('summary.json', 'w'))

def num(x, d=0):
    """Numero no formato brasileiro: espaco fino nos milhares, virgula decimal."""
    return f'{x:,.{d}f}'.replace(',', '\\u202f').replace('.', ',')

print('Ambiente pronto.')
''')

md('''
## Arquivos de entrada

| Arquivo | Origem |
|---------|--------|
| `minim.mdp`, `ions.mdp`, `npt_prot.mdp`, `md_prot.mdp`, `md_membrane.mdp` | pasta `mdp_martini/` deste repositório |
| `martini_v3.0.0.itp` (tipos de bead e matriz de Lennard-Jones) | [M3-Lipid-Parameters](https://github.com/Martini-Force-Field-Initiative/M3-Lipid-Parameters) |
| `martini_v3.0.0_ffbonded_v2.itp`, `martini_v3.0.0_phospholipids_PC_v2.itp` | parâmetros lipídicos v2 (2025), mesmo repositório |
| `martini_v3.0.0_solvents_v1.itp`, `martini_v3.0.0_ions_v1.itp` | água W e íons, mesmo repositório |
| `water.gro` | caixa de 400 beads W equilibrada, biblioteca cgmartini.nl |
| `viz/` | visualização das etapas (`simple_view.py`), das réplicas (`replica_views.py`) e gráficos |
''')

code('''
import os, sys
REPO = 'https://raw.githubusercontent.com/madsondeluna/curso-dinamica-molecular-202x/main'

def refresh(url, dest):
    """Baixa sempre a versao atual do GitHub; sem rede, mantem a copia que ja existir."""
    os.makedirs(os.path.dirname(dest) or '.', exist_ok=True)
    os.system(f'wget -q {url} -O {dest}.new')
    if os.path.exists(f'{dest}.new') and os.path.getsize(f'{dest}.new') > 0:
        os.replace(f'{dest}.new', dest)
    elif os.path.exists(f'{dest}.new'):
        os.remove(f'{dest}.new')
    return os.path.exists(dest)

# Parametros e codigo do visualizador acompanham a versao do notebook: baixados a cada execucao
mdp_files = ['minim.mdp', 'ions.mdp', 'npt_prot.mdp', 'md_prot.mdp', 'md_membrane.mdp', 'md_membrane_eq.mdp']
viz_files = ['viz/simple_view.py', 'viz/replica_plots.py', 'viz/replica_views.py',
             'viz/pure/__init__.py', 'viz/pure/palette.py', 'viz/pure/mpl.py']
ok = [refresh(f'{REPO}/mdp_martini/{f}', f) for f in mdp_files] + [refresh(f'{REPO}/{f}', f) for f in viz_files]
print('  arquivos .mdp e visualizador:', 'atualizados' if all(ok) else 'FALTANDO, verifique a conexao')

itp_files = ['martini_v3.0.0.itp', 'martini_v3.0.0_ffbonded_v2.itp',
             'martini_v3.0.0_phospholipids_PC_v2.itp',
             'martini_v3.0.0_solvents_v1.itp', 'martini_v3.0.0_ions_v1.itp']
if not cached(itp_files, 'campo de forca Martini 3'):
    BASE = 'https://raw.githubusercontent.com/Martini-Force-Field-Initiative/M3-Lipid-Parameters/main/ITPs'
    for f in itp_files:
        os.system(f'wget -q {BASE}/{f} -O {f}')

if not cached(['water.gro'], 'caixa de agua W'):
    os.system('wget -q https://cgmartini-library.s3.ca-central-1.amazonaws.com/1_Downloads/example_applications/solvent_systems/water.gro')

# Replica ja simulada: o pacote traz as saidas de todas as etapas
# a lista inclui os arquivos das versoes mais novas do pacote: uma copia antiga no Drive e baixada de novo
PACOTE = ['summary.json', 'md_prot_pbc.xtc', 'mem_whole.xtc']
if USE_PRECOMPUTED and not cached(PACOTE, f'replica {REPLICA} pre-computada'):
    os.system(f'wget -q {REPO}/results/rep{REPLICA}.tar.gz -O rep.tar.gz && tar xzf rep.tar.gz && rm rep.tar.gz')
    print('Replica baixada.' if os.path.exists('summary.json') else 'Replica nao encontrada no GitHub: o notebook vai simular.')

# Visualizador e tema dos graficos
sys.path.insert(0, os.path.abspath('viz'))
import re
import numpy as np
import matplotlib.pyplot as plt
from pure import palette
from replica_plots import use_paper_style, finish, COLORS
import simple_view
from simple_view import stage, card, style_beads, CHART, INK, BEADS, SS_COLORS
use_paper_style()                       # mesmo estilo das figuras do README
C1 = COLORS['protein']                  # proteina: laranja
CM = COLORS['membrane']                 # membrana: azul

def read_xvg(filepath, ncol=2):
    data = []
    with open(filepath) as f:
        for line in f:
            if not line.startswith(('#', '@')):
                parts = line.split()
                if len(parts) >= ncol:
                    try: data.append([float(x) for x in parts[:ncol]])
                    except ValueError: pass
    return np.array(data)

def xpm_to_distance_matrix(fname):
    """Le a matriz de um .xpm do GROMACS: cada cor tem o valor num comentario /* "0.0385" */."""
    lines = open(fname).read().splitlines()
    head = next(i for i, l in enumerate(lines) if re.match(r'^"\\s*\\d+\\s+\\d+\\s+\\d+\\s+\\d+', l))
    w, h, nc, cpp = map(int, lines[head].strip('",').split())
    value = {}
    for l in lines[head + 1: head + 1 + nc]:
        key = l[1:1 + cpp]
        m = re.search(r'/\\*\\s*"([^"]*)"\\s*\\*/', l)
        value[key] = float(m.group(1)) if m else float(len(value))
    rows = [l.split('"')[1] for l in lines[head + 1 + nc:] if l.startswith('"')][:h]
    arr = np.array([[value[r[c * cpp:(c + 1) * cpp]] for c in range(w)] for r in rows])
    return arr[::-1]          # a primeira linha do .xpm e o topo da imagem

def count_particles(gro):
    return int(open(gro).readlines()[1])

def count_res(gro, resname):
    return sum(1 for l in open(gro).readlines()[2:-1] if l[5:10].strip() == resname)

def box_edge(gro):
    return float(open(gro).readlines()[-1].split()[0])

print('Arquivos prontos.')
''')

md('''
# Fundamentos: do átomo ao bead

Antes de simular, o que a dinâmica molecular faz, o que um campo de força contém e o que muda quando trocamos de resolução. Todas as figuras deste capítulo estão em `imagens/fund_*.png` e são geradas por `tools/make_figures.py`.
''')

md('''
## 1. O ciclo da dinâmica molecular

A MD resolve as equações de Newton para todas as partículas, passo a passo. A cada passo $\\Delta t$, o programa calcula as forças a partir da energia potencial, atualiza as velocidades e depois as posições. O GROMACS usa o integrador leap-frog, que calcula velocidades na metade do passo.

![Ciclo da MD](https://raw.githubusercontent.com/madsondeluna/curso-dinamica-molecular-202x/main/imagens/fund_ciclo_md.png)
''')

md('''
## 2. O que há dentro de um campo de força

A energia potencial é a soma de termos simples. Os **ligados** mantêm a geometria da molécula (ligações, ângulos, diedros). Os **não ligados** agem entre pares de partículas: Lennard-Jones (repulsão de curto alcance e atração de van der Waals) e Coulomb (cargas). Um campo de força é o conjunto de parâmetros de todos esses termos, para cada tipo de partícula.

$$U = \\sum_{lig} \\tfrac{1}{2}k_b(r-r_0)^2 + \\sum_{ang} \\tfrac{1}{2}k_\\theta(\\theta-\\theta_0)^2 + \\sum_{died} k_\\phi[1+\\cos(n\\phi-\\delta)] + \\sum_{i \\lt j} 4\\varepsilon_{ij}\\left[\\left(\\tfrac{\\sigma_{ij}}{r_{ij}}\\right)^{12}-\\left(\\tfrac{\\sigma_{ij}}{r_{ij}}\\right)^{6}\\right] + \\sum_{i \\lt j} \\frac{q_iq_j}{4\\pi\\varepsilon_0\\varepsilon_r r_{ij}}$$

![Termos do campo de força](https://raw.githubusercontent.com/madsondeluna/curso-dinamica-molecular-202x/main/imagens/fund_campo_de_forca.png)
''')

md('''
## 3. Por que o passo de tempo é 2 fs (ou 20 fs)

O passo de tempo é limitado pelo movimento mais rápido do sistema. Em all-atom, é o estiramento de ligações com hidrogênio (período de ~11 fs); com essas ligações restritas (`constraints = h-bonds`), 2 fs é estável. No Martini não há hidrogênios e as ligações entre beads vibram devagar: 20 fs funcionam.

![Passo de tempo](https://raw.githubusercontent.com/madsondeluna/curso-dinamica-molecular-202x/main/imagens/fund_passo_de_tempo.png)
''')

md('''
## 4. Condições periódicas e cutoff

A caixa simulada se repete infinitamente em todas as direções, o que elimina superfícies artificiais. Cada partícula só interage com as vizinhas dentro de uma distância de corte $r_c$, que pode incluir cópias de partículas da caixa ao lado. Por isso a caixa precisa ser maior que $2r_c$ e a molécula não pode "ver" a própria imagem.

![PBC e cutoff](https://raw.githubusercontent.com/madsondeluna/curso-dinamica-molecular-202x/main/imagens/fund_pbc_cutoff.png)
''')

md('''
## 5. Quatro resoluções da mesma molécula

| Resolução | O que é representado | Exemplo |
|-----------|----------------------|---------|
| QM (mecânica quântica) | núcleos e elétrons; a ligação química emerge da equação de Schrödinger | DFT, Hartree-Fock |
| All-atom (AA) | cada átomo, inclusive todos os H; elétrons embutidos em cargas parciais | CHARMM36m, AMBER |
| United-atom (UA) | átomos pesados; CH, CH₂ e CH₃ alifáticos viram uma partícula | GROMOS 54A7 |
| Coarse-grained (CG) | beads de 2 a 4 átomos pesados | Martini 3 |

Só a QM descreve reações químicas (quebra e formação de ligações). As outras três usam ligações fixas e mudam apenas quantos detalhes guardam.

![Fenilalanina em quatro resoluções](https://raw.githubusercontent.com/madsondeluna/curso-dinamica-molecular-202x/main/imagens/fund_resolucoes_molecula.png)

O mapeamento Martini 3 depende do resíduo: glicina tem só o bead de backbone, aminoácidos aromáticos têm três ou quatro beads de cadeia lateral.

![Mapeamento de resíduos](https://raw.githubusercontent.com/madsondeluna/curso-dinamica-molecular-202x/main/imagens/fund_mapeamento_residuos.png)
''')

md('''
## 6. Um bead não é um átomo grande

Quatro moléculas de água viram um bead W. O bead é maior (σ = 0,47 nm contra 0,315 nm do oxigênio do TIP3P) e o poço de Lennard-Jones é mais fundo, porque precisa representar a atração média de quatro moléculas.

![Mapeamento 4:1 e Lennard-Jones](https://raw.githubusercontent.com/madsondeluna/curso-dinamica-molecular-202x/main/imagens/fund_mapeamento_lj.png)
''')

md('''
## 7. Quanto cada resolução alcança

Cada resolução cobre uma faixa de tamanhos e tempos. As faixas abaixo são ordens de grandeza típicas, não limites rígidos. Processos de interesse biológico, como enovelamento e fusão de membranas, ficam no limite ou além do all-atom e dentro do coarse-grained.

![Escalas](https://raw.githubusercontent.com/madsondeluna/curso-dinamica-molecular-202x/main/imagens/fund_escalas.png)

O custo por nanossegundo é o produto de três fatores: partículas, vizinhos dentro do cutoff e passos por ns. A figura usa o sistema da Parte A; os valores de all-atom e united-atom são contagens teóricas.

![Custo teórico](https://raw.githubusercontent.com/madsondeluna/curso-dinamica-molecular-202x/main/imagens/fund_custo_teorico.png)
''')

md('''
## 8. Paisagem de energia livre

A simulação caminha sobre uma paisagem de energia livre: passa a maior parte do tempo nos vales e cruza as barreiras raramente. Ao agrupar átomos em beads, o CG remove a rugosidade (barreiras pequenas de rotâmeros, ligações de hidrogênio, rearranjo de águas) e mantém os vales principais.

![Paisagem 2D](https://raw.githubusercontent.com/madsondeluna/curso-dinamica-molecular-202x/main/imagens/fund_paisagem_2d.png)

A QM, à esquerda, tem os mesmos vales do all-atom e uma rugosidade ainda mais fina, vinda dos elétrons (polarização, transferência de carga).

No enovelamento, a paisagem tem forma de funil: larga em cima (muitas conformações desenoveladas) e estreita no fundo (o estado nativo). A rede elástica do Martini troca o funil por um poço único: a proteína fica perto da estrutura inicial e não enovela nem desenovela.

![Funil de energia](https://raw.githubusercontent.com/madsondeluna/curso-dinamica-molecular-202x/main/imagens/fund_funil.png)

Minimizar e simular são coisas diferentes: a minimização só desce até o mínimo mais próximo; a dinâmica a temperatura finita explora a paisagem.

![Minimização e dinâmica](https://raw.githubusercontent.com/madsondeluna/curso-dinamica-molecular-202x/main/imagens/fund_minimizacao_vs_md.png)
''')

md('''
## 9. Tempos em resoluções diferentes não estão na mesma escala

Numa paisagem lisa, a mesma barreira é cruzada mais vezes no mesmo número de passos. Isso significa que **1 ns de Martini não equivale a 1 ns de all-atom**: o relógio da simulação CG anda mais rápido que o relógio físico. Para a água, Marrink *et al.* (2007) estimaram um fator de cerca de 4, mas o fator muda de processo para processo. Por isso, compare resoluções diferentes por eventos (formação da bicamada, contatos, estados visitados) e não por nanossegundos.

![Relógios diferentes](https://raw.githubusercontent.com/madsondeluna/curso-dinamica-molecular-202x/main/imagens/fund_relogios.png)

Com o mesmo ruído térmico nas quatro paisagens, as travessias da barreira em 1 milhão de passos caem de 52 no coarse-grained para 32 no united-atom, 18 no all-atom e 11 na QM: cada camada de rugosidade desacelera o relógio. Na QM, cada passo também custa muito mais, porque resolve os elétrons.
''')

md('''
## O que muda no Martini

O Martini agrupa, em média, **4 átomos pesados em 1 bead**. Quatro moléculas de água viram um único bead `W`. Com menos partículas e um potencial mais suave, o passo de integração sobe de 2 fs para **20 fs**.

| Parâmetro | All-atom (CHARMM36m) | United-atom (GROMOS 54A7) | Coarse-grained (Martini 3) |
|-----------|----------------------|---------------------------|----------------------------|
| Partícula | átomo | átomo pesado ou grupo CHₙ | bead de ~4 átomos pesados |
| Água | TIP3P, 3 sítios por molécula | SPC, 3 sítios por molécula | W, 1 bead = 4 moléculas |
| Passo de tempo | 2 fs | 2 fs | 20 fs (10 fs na equilibração) |
| Eletrostática | PME | PME ou reaction-field | reaction-field, `epsilon_r = 15` |
| Cutoff | 1,2 nm | 1,4 nm | 1,1 nm |
| Restrições | ligações com H | ligações com H | nenhuma |
| Estrutura secundária | emerge do campo de força | emerge do campo de força | fixada na topologia + rede elástica |

```
coulombtype   = reaction-field
rcoulomb      = 1.1      # cutoff eletrostatico (nm)
epsilon_r     = 15       # blindagem dieletrica implicita
epsilon_rf    = 0        # dieletrico infinito alem do cutoff
vdw-modifier  = Potential-shift-verlet
rvdw          = 1.1      # cutoff van der Waals (nm)
verlet-buffer-tolerance = 0.005
```

> **Tempo efetivo:** a superfície de energia do Martini é mais lisa que a atomística, e a difusão é mais rápida. Os tempos simulados não correspondem 1:1 a tempos reais; trate-os como ordens de grandeza.
''')

md('''
## Paisagem energética: por que o CG é mais rápido do que o número de partículas sugere

A **paisagem energética** é a energia livre do sistema em função das suas coordenadas. Uma simulação de MD caminha sobre ela: passa a maior parte do tempo nos vales (estados estáveis) e cruza as barreiras (transições) raramente, com probabilidade proporcional a $e^{-\\Delta G^{\\ddagger}/k_BT}$.

Ao agrupar átomos em beads, o modelo CG **integra** os graus de liberdade removidos: rotação de cada C-H, quebra e formação de cada ligação de hidrogênio, rearranjo de cada molécula de água. O que sobra é um potencial de força média, mais liso:

| Efeito | Consequência na simulação |
|--------|---------------------------|
| Desaparecem as barreiras pequenas e numerosas (rugosidade) | menos tempo preso em mínimos locais |
| Atrito entre partículas diminui | difusão mais rápida; Marrink *et al.* (2007) estimaram um fator de cerca de 4 a partir da difusão da água |
| Vales e barreiras grandes, definidos pela hidrofobicidade, permanecem | automontagem de membrana, partição e agregação continuam corretas |
| Entropia dos graus removidos fica embutida nos parâmetros | o balanço entalpia-entropia vale para a faixa de temperatura da parametrização |
| Na proteína, a rede elástica cria um vale único em torno da estrutura inicial | sem enovelamento nem troca de estrutura secundária |

O gráfico abaixo é um **esquema ilustrativo** (funções analíticas, não resultado de simulação) da mesma coordenada de reação nas quatro resoluções. A curva quântica (QM) é a mais rugosa: além de tudo o que o all-atom descreve, ela resolve os elétrons (polarização, formação e quebra de ligações), e cada passo de descida de resolução alisa a paisagem.
''')

code('''
# Esquema ilustrativo: paisagens 1D construidas com funcoes analiticas
x = np.linspace(-1.6, 1.6, 1200)
base = 4 * (x**2 - 1)**2                                   # dois vales e uma barreira
rough_aa = 0.9 * np.sin(23 * x) + 0.6 * np.sin(41 * x + 1.3)
rough_ua = 0.6 * np.sin(19 * x + 0.4) + 0.3 * np.sin(37 * x)
rough_qm = 0.7 * np.sin(67 * x + 0.7) + 0.5 * np.sin(89 * x + 2.1)   # detalhe eletronico: polarizacao, ligacoes
curves = {'Quântica (QM)': base + rough_aa + rough_qm + 2.6 * np.exp(-x**2 / 0.03),
          'All-atom': base + rough_aa + 1.2 * np.exp(-x**2 / 0.05),
          'United-atom': base + rough_ua + 0.8 * np.exp(-x**2 / 0.06),
          'Coarse-grained': 0.6 * base}
fig, ax = plt.subplots(figsize=(10, 3.8))
for (name, y), color in zip(curves.items(), palette.series(4)):
    ax.plot(x, y - y.min(), color=color, lw=1.6, label=name)
finish(ax, title='Mesma coordenada em quatro resoluções',
       xlabel='Coordenada de reação', ylabel='Energia livre (unidades demonstrativas)', legend=False)
ax.legend(frameon=False, loc='upper left', bbox_to_anchor=(1.01, 1))
ax.set_yticks([])
plt.show()
''')

# ---------------------------------------------------------------- Parte A
md('''
# Parte A: Ubiquitina em Água (1UBQ)

**Nesta parte você vai:** converter uma proteína real em modelo Martini, colocá-la numa caixa de água com sal, simular 15 ns e medir quanto cada região se move.

**Por que a ubiquitina.** É uma proteína pequena (76 resíduos) presente em todas as células eucarióticas. Ela funciona como uma etiqueta: enzimas ligam a Gly76 da cauda C-terminal a lisinas de outras proteínas. Cadeias de ubiquitina ligadas pela Lys48 marcam a proteína para ser degradada pelo proteassomo. Esse sistema de degradação é alvo terapêutico: o bortezomibe, inibidor do proteassomo, foi aprovado em 2003 para o tratamento do mieloma múltiplo (Bross *et al.*, 2004). Na análise, você vai ver que a cauda C-terminal é justamente a parte mais móvel da proteína.

Mesma proteína do curso atomístico (`legado/MD_GMX.ipynb`): 76 resíduos, carga líquida zero.
''')

code('''
import os
if not cached(['1UBQ_protein.pdb'], '1UBQ_protein.pdb'):
    os.system('wget -q https://files.rcsb.org/download/1UBQ.pdb')
    # Manter apenas atomos da proteina (remove aguas cristalograficas HETATM)
    with open('1UBQ.pdb') as f:
        prot = [l for l in f if l.startswith('ATOM')]
    with open('1UBQ_protein.pdb', 'w') as f:
        f.writelines(prot + ['END\\n'])
N_RES = len({l[22:26] for l in open('1UBQ_protein.pdb') if l.startswith('ATOM')})
N_HEAVY = sum(1 for l in open('1UBQ_protein.pdb') if l.startswith('ATOM'))

def stick_style(view):
    view.setStyle({}, {'stick': {'color': INK, 'radius': 0.15}})
    # cartoon pela estrutura secundaria que o 3Dmol calcula das ligacoes de hidrogenio do backbone
    view.addStyle({}, {'cartoon': {'colorscheme': {'prop': 'ss', 'map': {k: c for k, (_, c) in SS_COLORS.items()}},
                                   'opacity': 0.9}})

stage('A', 1, 'Estrutura cristalográfica', structure='1UBQ_protein.pdb', style=stick_style, box=False,
      text=['A estrutura 1UBQ vem de cristalografia de raios X. A difração não resolve hidrogênios, '
            'então o arquivo traz apenas os átomos pesados (C, N, O, S). As águas cristalográficas já foram removidas.',
            'É o mesmo ponto de partida nas três resoluções: o que muda é o que fazemos com estes átomos no passo seguinte.'],
      legend=[('Átomo pesado', 'ink'), 'SS'],
      facts={'resíduos': N_RES, 'átomos pesados': N_HEAVY, 'hidrogênios': 0, 'carga líquida': 0},
      compare={'aa': ('todos os H', '', 'O pdb2gmx adicionaria todos os hidrogênios: a proteína passa a ter cerca do dobro de partículas.'),
               'ua': ('H polares e aromáticos', '', 'Só os hidrogênios polares (N-H, O-H) e aromáticos ficam explícitos; os demais se fundem ao carbono.'),
               'cg': ('martinize2 + Martini 3', '', 'Não adiciona nada: agrupa os átomos pesados em beads.')})
''')

md('''
## Passo A1: Conversão para Martini 3 com martinize2

`martinize2` lê a estrutura atomística, agrupa os átomos de cada resíduo em beads segundo os mapeamentos do Martini 3 e escreve a topologia.

| Opção | Significado |
|-------|-------------|
| `-ff martini3001` | campo de força Martini 3.0.0 para proteínas |
| `-dssp` | calcula a estrutura secundária com o DSSP do `mdtraj` |
| `-elastic -ef 700 -el 0.5 -eu 0.9` | rede elástica: molas de 700 kJ/mol/nm² entre beads BB distantes 0,5 a 0,9 nm |
| `-cys auto` | detecta pontes dissulfeto (a ubiquitina não tem) |
| `-p backbone` | gera restrições de posição para os beads BB (`#ifdef POSRES`) |

**Por que a rede elástica?** No Martini, a estrutura secundária é imposta pela topologia e não emerge das interações. Sem molas entre os beads do backbone, a estrutura terciária se desfaz em poucos nanossegundos.

As comparações com all-atom e united-atom são **teóricas**: o número de partículas vem da composição química de cada resíduo, e o custo, de contagens de partículas, vizinhos e passos.

**Na aula:** como cada resíduo vira beads, e a mesma fenilalanina em quatro resoluções.

![Mapeamento de resíduos](https://raw.githubusercontent.com/madsondeluna/curso-dinamica-molecular-202x/main/imagens/fund_mapeamento_residuos.png)

![Quatro resoluções](https://raw.githubusercontent.com/madsondeluna/curso-dinamica-molecular-202x/main/imagens/fund_resolucoes_molecula.png)
''')

code('''
import os
if not cached(['ubq_cg.pdb', 'molecule_0.itp'], 'martinize2'):
    os.system('martinize2 -f 1UBQ_protein.pdb -x ubq_cg.pdb -o martinize.top '
              '-ff martini3001 -dssp -elastic -ef 700 -el 0.5 -eu 0.9 '
              '-cys auto -p backbone 2>&1 | grep -E "WARN|ERROR"')

# Contagem teorica de particulas em all-atom e united-atom, sem simular nada.
# H por residuo no interior da cadeia, pH 7 (Lys e Arg +1, Asp e Glu -1, His neutra)
H_ALL = dict(ALA=5, ARG=13, ASN=6, ASP=4, CYS=5, GLN=8, GLU=6, GLY=3, HIS=7, ILE=11, LEU=11,
             LYS=13, MET=9, PHE=9, PRO=7, SER=5, THR=7, TRP=10, TYR=9, VAL=9)
# H que um campo united-atom (GROMOS 54A7) mantem explicitos: polares (N-H, O-H, S-H) e aromaticos
H_UA = dict(ALA=1, ARG=6, ASN=3, ASP=1, CYS=2, GLN=3, GLU=1, GLY=1, HIS=4, ILE=1, LEU=1,
            LYS=4, MET=1, PHE=6, PRO=0, SER=2, THR=2, TRP=7, TYR=6, VAL=1)
residues = []
for l in open('1UBQ_protein.pdb'):
    if l.startswith('ATOM') and (not residues or residues[-1][0] != l[21:27]):
        residues.append((l[21:27], l[17:20]))
names = [r for _, r in residues]
N_AA = N_HEAVY + sum(H_ALL[r] for r in names) + 2   # N-terminal NH3+: 2 H a mais
N_UA = N_HEAVY + sum(H_UA[r] for r in names) + 2

SS = next((l[2:].strip() for l in open('molecule_0.itp')
           if l.startswith('; ') and len(l.strip()) > 20 and set(l[2:].strip()) <= set('HBEGITSCP')), '?')
N_CG = sum(1 for l in open('ubq_cg.pdb') if l.startswith('ATOM'))
N_RB = sum(1 for l in open('molecule_0.itp') if 'RUBBER_FC' in l and not l.startswith(('#', ';')))
print('Estrutura secundaria (DSSP):', SS)
print(f'All-atom (teorico): {N_AA} | united-atom (teorico): {N_UA} | Martini 3: {N_CG} | molas da rede elastica: {N_RB}')
''')

code('''
# Cartao 1: mapeamento atomistico -> CG, estruturas sobrepostas
def mapping_style(view):
    view.setStyle({'model': 0}, {'stick': {'color': INK, 'radius': 0.12}})
    view.setStyle({'model': 1, 'atom': 'BB'}, {'sphere': {'color': BEADS['BB'][1], 'radius': 1.4, 'opacity': 0.85}})
    view.setStyle({'model': 1, 'not': {'atom': 'BB'}}, {'sphere': {'color': BEADS['SC'][1], 'radius': 1.2, 'opacity': 0.85}})

import py3Dmol
v = py3Dmol.view(width=820, height=440)
v.addModel(open('1UBQ_protein.pdb').read(), 'pdb')
v.addModel(open('ubq_cg.pdb').read(), 'pdb')
mapping_style(v)
v.setBackgroundColor('#ffffff', 1)
v.zoomTo()

stage('A', 2, 'Mapeamento atomístico para beads', figure_html=v._make_html(),
      text=['Cada esfera é um bead Martini desenhado sobre a estrutura atomística original. '
            'O backbone de cada resíduo vira um bead BB; a cadeia lateral vira de zero (glicina) a cinco beads (triptofano).',
            f'A estrutura secundária calculada pelo DSSP fica gravada na topologia, e {N_RB} molas ligam beads BB '
            'próximos para manter a forma da proteína.'],
      legend=['BB', 'SC', ('Átomo pesado', 'ink')],
      facts={'átomos pesados': N_HEAVY, 'beads': N_CG, 'átomos por bead': num(N_HEAVY / N_CG, 1), 'molas elásticas': N_RB},
      compare={'aa': ('todos os H explícitos', f'{num(N_AA)} átomos', 'Resolve cada ligação de hidrogênio e a geometria de cada grupo. Contagem a pH 7.'),
               'ua': ('H polares e aromáticos', f'{num(N_UA)} partículas', 'CH, CH₂ e CH₃ alifáticos viram uma partícula cada (convenção GROMOS 54A7).'),
               'cg': ('Martini 3', f'{num(N_CG)} beads', f'{num(N_AA / N_CG, 1)} vezes menos partículas que o all-atom. Sem hidrogênios.')})
''')

md('''
## Passo A2: Topologia, caixa e solvatação

O `topol.top` inclui primeiro a matriz de interações do Martini 3, depois a proteína (`molecule_0.itp`), o solvente e os íons. A caixa cúbica deixa 1,2 nm entre a proteína e a borda. O `gmx solvate` preenche o volume com cópias de `water.gro`; `-radius 0.21` aproxima o raio de um bead Martini (o padrão de 0,105 nm é para átomos).

**Na aula:** a caixa periódica e por que um bead W não é um átomo grande.

![PBC e cutoff](https://raw.githubusercontent.com/madsondeluna/curso-dinamica-molecular-202x/main/imagens/fund_pbc_cutoff.png)

![Mapeamento 4:1 e Lennard-Jones](https://raw.githubusercontent.com/madsondeluna/curso-dinamica-molecular-202x/main/imagens/fund_mapeamento_lj.png)
''')

code('''
import os
TOP = """#include "martini_v3.0.0.itp"
#include "molecule_0.itp"
#include "martini_v3.0.0_solvents_v1.itp"
#include "martini_v3.0.0_ions_v1.itp"

[ system ]
Ubiquitin Martini 3

[ molecules ]
molecule_0    1
"""
if not cached(['ubq_solv.gro'], 'solvate'):
    # Reescrever topol.top garante que re-execucoes nao acumulem W/NA/CL
    with open('topol.top', 'w') as f:
        f.write(TOP)
    os.system('gmx editconf -f ubq_cg.pdb -o ubq_box.gro -c -d 1.2 -bt cubic -quiet')
    os.system('gmx solvate -cp ubq_box.gro -cs water.gro -radius 0.21 -o ubq_solv.gro -p topol.top -quiet')

N_W = count_res('ubq_solv.gro', 'W')
L = box_edge('ubq_solv.gro')
N_SOLV = count_particles('ubq_solv.gro')
stage('A', 3, 'Solvatação', structure='ubq_solv.gro',
      text=['O <code>gmx solvate</code> preenche a caixa com cópias de <code>water.gro</code>, uma caixa de água W já equilibrada.',
            f'Cada bead W representa quatro moléculas de água: os {num(N_W)} beads equivalem a {num(4 * N_W)} moléculas de H₂O.'],
      legend=['BB', 'SC', 'W'],
      facts={'beads de proteína': N_CG, 'beads W': num(N_W), 'aresta da caixa (nm)': num(L, 2), 'partículas no total': num(N_SOLV)},
      compare={'aa': ('TIP3P, 3 sítios', f'≈ {num(N_AA + 12 * N_W)} átomos', 'A mesma água em all-atom: 3 átomos por molécula. A água passa a ser mais de 95% do sistema.'),
               'ua': ('SPC, 3 sítios', f'≈ {num(N_UA + 12 * N_W)} partículas', 'A água continua com 3 sítios. O united-atom só economiza na proteína, e a água domina o custo.'),
               'cg': ('W, 1 bead = 4 H₂O', f'{num(N_SOLV)} beads', 'A redução maior está no solvente: 12 sítios atomísticos viram 1 bead.')})
''')

md('''
## Passo A3: Íons (0,15 M NaCl)

O `genion -conc` conta cada bead W como uma molécula de água e subestimaria a concentração por um fator 4. Por isso o número de pares é calculado a partir do volume da caixa:

$$N_{pares} = c \\cdot V \\cdot N_A = 0{,}15\\ \\text{mol/L} \\times V\\,[\\text{nm}^3] \\times 0{,}6022$$
''')

code('''
import os
V = box_edge('ubq_solv.gro') ** 3
N_PAIRS = round(0.15 * V * 0.6022)
if not cached(['ubq_ions.gro'], 'genion'):
    os.system('gmx grompp -f ions.mdp -c ubq_solv.gro -p topol.top -o ions.tpr -quiet')
    os.system(f'echo W | gmx genion -s ions.tpr -o ubq_ions.gro -p topol.top -seed {tool_seed("genion")} '
              f'-pname NA -nname CL -np {N_PAIRS} -nn {N_PAIRS} -neutral -quiet 2>&1 | tail -2')

stage('A', 4, 'Íons', structure='ubq_ions.gro',
      text=[f'O <code>genion</code> troca {2 * N_PAIRS} beads W por íons. A água aparece translúcida para os íons ficarem visíveis.',
            'Com <code>-conc 0.15</code>, o genion contaria 1 bead W como 1 molécula e poria um quarto dos íons. '
            'O número de pares vem do volume da caixa.'],
      legend=['BB', 'SC', 'NA', 'CL'],
      facts={'volume (nm³)': num(V, 1), 'pares NaCl': N_PAIRS, 'concentração (M)': '0,15', 'carga líquida': 0},
      compare={'aa': ('Na⁺ e Cl⁻ de 1 átomo', '', 'O genion -conc funciona direto: cada água é uma molécula. A camada de hidratação é feita de águas explícitas.'),
               'ua': ('Na⁺ e Cl⁻ de 1 átomo', '', 'Igual ao all-atom: íons e águas explícitos.'),
               'cg': ('bead TQ5 carregado', f'{N_PAIRS} + {N_PAIRS}', 'O íon é um bead pequeno. As águas ao redor não são resolvidas uma a uma.')})
''')

md('''
## Passo A4: Minimização de energia

Steepest descent até a força máxima ficar abaixo de 1000 kJ/mol/nm. A rede elástica e as cadeias laterais com geometria fixa impedem forças muito menores; o critério de 500 kJ/mol/nm do exercício atomístico não converge aqui.

![Minimização e dinâmica](https://raw.githubusercontent.com/madsondeluna/curso-dinamica-molecular-202x/main/imagens/fund_minimizacao_vs_md.png)
''')

code('''
import os
run_md('em', 'minimizacao', 'minim.mdp', 'ubq_ions.gro', 'topol.top')
gmx_cached(['em_potential.xvg'], 'energia da minimizacao',
           'echo Potential | gmx energy -f em.edr -o em_potential.xvg')

d = read_xvg('em_potential.xvg')
fig, ax = plt.subplots(figsize=(7, 3.5))
ax.plot(d[:, 0], d[:, 1], color=C1, lw=1.2)
finish(ax, title='Minimização de energia', subtitle='kJ/mol', xlabel='Passo', ylabel='Energia potencial')
plt.show()

''')

md('''
A energia caiu e estabilizou: as sobreposições entre beads foram removidas. Abaixo, o sistema minimizado, com a água translúcida.
''')

code('''
log = open('em.log').read()
fmax = float(re.search(r'Maximum force\\s+=\\s+([\\d.eE+-]+)', log).group(1))
stage('A', 5, 'Minimização de energia', structure='em.gro',
      text=['O steepest descent move cada partícula na direção da força até remover sobreposições entre beads. '
            'Não é dinâmica: não há tempo nem temperatura.',
            'A estrutura quase não muda; o que muda é a energia, que cai até um platô.'],
      legend=['BB', 'SC', 'NA', 'CL'],
      facts={'passos': len(d), 'energia final (kJ/mol)': num(d[-1, 1]), 'força máxima (kJ/mol/nm)': num(fmax), 'critério': '1000'},
      compare={'aa': ('emtol = 500', '', 'No curso atomístico o critério é 500 kJ/mol/nm, com ~30 mil átomos.'),
               'ua': ('emtol = 500 a 1000', '', 'Mesmo algoritmo, sistema do tamanho do all-atom por causa da água.'),
               'cg': ('emtol = 1000', f'{len(d)} passos', 'Menos partículas e forças mais suaves: segundos de CPU.')})
''')

md('''
## Passo A5: Equilibração NPT com restrição de posição (0,5 ns)

Os beads BB ficam presos às posições minimizadas (`-DPOSRES`, 1000 kJ/mol/nm²) enquanto água e íons relaxam. O passo de 10 fs é metade do valor de produção.

| Parâmetro | Valor |
|-----------|-------|
| Termostato | v-rescale, 310 K, `tau_t = 1 ps`, grupos Protein e non-Protein |
| Barostato | C-rescale isotrópico, 1 bar, `tau_p = 12 ps`, `compressibility = 3e-4 bar⁻¹` |

> **Tempo de execução estimado no Colab:** 2 a 4 minutos.
''')

code('''
import os
run_md('npt', 'NPT (0,5 ns)', 'npt_prot.mdp', 'em.gro', 'topol.top', extra='-r em.gro')
gmx_cached(['npt_energy.xvg'], 'temperatura e densidade',
           'printf "Temperature\\\\nDensity\\\\n" | gmx energy -f npt.edr -o npt_energy.xvg')

d = read_xvg('npt_energy.xvg', ncol=3)
fig, ax = plt.subplots(figsize=(8, 3.5))
ax.plot(d[:, 0], d[:, 1], color=C1, lw=0.8)
finish(ax, title='Temperatura durante o NPT', ylabel='Temperatura', subtitle='K', xlabel='Tempo (ps)')
plt.show()
''')

md('''
**O que se observa:** a temperatura oscila em torno de 310 K.

**Significado:** 310 K é a temperatura do corpo humano. O termostato v-rescale mantém a energia cinética média no alvo; as oscilações são as flutuações naturais de um sistema pequeno.
''')

code('''
fig, ax = plt.subplots(figsize=(8, 3.5))
ax.plot(d[:, 0], d[:, 2], color=C1, lw=0.8)
finish(ax, title='Densidade durante o NPT', ylabel='Densidade', subtitle='kg/m³', xlabel='Tempo (ps)')
plt.show()
print(f'Densidade media na segunda metade: {d[len(d)//2:, 2].mean():.0f} kg/m3')
''')

md('''
**O que se observa:** a densidade muda no começo e se estabiliza.

**Significado:** a caixa se ajustou à pressão de 1 bar. Sem essa etapa, a água ficaria comprimida ou rarefeita, e as interações da proteína com o solvente estariam erradas.
''')

md('''
Temperatura e densidade estabilizaram. Abaixo, o sistema ao fim da equilibração, com a água translúcida.
''')

code('''
stage('A', 6, 'Equilibração NPT', structure='npt.gro',
      text=['Temperatura e pressão chegam aos valores alvo enquanto o backbone fica preso. '
            'A caixa encolhe ou cresce um pouco até a densidade estabilizar.',
            'A densidade de água W Martini fica perto de 1000 kg/m³, como a água real.'],
      legend=['BB', 'SC', 'W', 'NA', 'CL'],
      facts={'temperatura média (K)': num(d[:, 1].mean(), 1), 'densidade (kg/m³)': num(d[len(d)//2:, 2].mean()),
             'aresta da caixa (nm)': num(box_edge('npt.gro'), 2), 'ns/dia': num(performance('npt.log') or 0)},
      compare={'aa': ('dt = 2 fs', '', 'NVT e depois NPT no curso atomístico. Ligações com H restritas pelo LINCS para permitir 2 fs.'),
               'ua': ('dt = 2 fs', '', 'Também restringe ligações com H; menos H do que no all-atom.'),
               'cg': ('dt = 10 fs aqui, 20 fs na produção', '', 'Sem hidrogênios, nenhuma vibração rápida limita o passo.')})
''')

md('''
## Passo A6: Simulação produtiva (15 ns)

Restrições removidas, passo de 20 fs: 750 000 passos.

> **Tempo de execução estimado no Colab:** 45 a 75 minutos. Com o Drive montado, a simulação é retomada do último checkpoint se a sessão cair.

**Na aula:** por que o Martini usa 20 fs e o all-atom 2 fs.

![Passo de tempo](https://raw.githubusercontent.com/madsondeluna/curso-dinamica-molecular-202x/main/imagens/fund_passo_de_tempo.png)
''')

code('''
import os
run_md('md_prot', 'MD produtiva (15 ns)', 'md_prot.mdp', 'npt.gro', 'topol.top', extra='-t npt.cpt')

# Indice com o grupo BB, reimagem e ajuste de rotacao/translacao
if not cached(['index_prot.ndx'], 'indice BB'):
    os.system('printf "a BB\\\\nq\\\\n" | gmx make_ndx -f md_prot.tpr -o index_prot.ndx -quiet > /dev/null 2>&1')
if not cached(['md_prot_pbc.xtc', 'md_prot_fit.xtc'], 'reimagem + ajuste'):
    os.system('printf "Protein\\\\nSystem\\\\n" | gmx trjconv -s md_prot.tpr -f md_prot.xtc -n index_prot.ndx '
              '-o md_prot_pbc.xtc -pbc mol -center -ur compact -quiet > /dev/null 2>&1')
    os.system('printf "BB\\\\nProtein\\\\n" | gmx trjconv -s md_prot.tpr -f md_prot_pbc.xtc -n index_prot.ndx '
              '-o md_prot_fit.xtc -fit rot+trans -quiet > /dev/null 2>&1')
if not cached(['md_prot_anim.pdb'], 'animacao proteina'):
    os.system('echo Protein | gmx trjconv -s md_prot.tpr -f md_prot_fit.xtc -n index_prot.ndx '
              '-o md_prot_anim.pdb -dt 250 -quiet > /dev/null 2>&1')

NSTEPS_PROT = int(re.search(r'nsteps\\s*=\\s*(\\d+)', open('md_prot.mdp').read()).group(1))
T_PROT = NSTEPS_PROT * 0.02 / 1000

PERF_CG = performance('md_prot.log')
from replica_views import export_protein_movie, single_view
pmovie, pmovie_water = export_protein_movie('.')       # filme completo do sistema: 0 a 15 ns, 31 quadros
stage('A', 7, 'Simulação produtiva',
      figure_html=lambda: single_view(pmovie, 'protein', animate=True, water_path=pmovie_water),
      text=[f'Filme completo dos {num(T_PROT)} ns, um quadro a cada 0,5 ns: proteína, íons e água translúcida, com a proteína centrada na caixa.',
            'As flutuações são pequenas: a rede elástica mantém o enovelamento. As regiões que se movem mais são as que têm menos molas.'],
      legend=['BB', 'SC', 'W', 'NA', 'CL'],
      facts={'tempo simulado (ns)': num(T_PROT), 'passo (fs)': 20, 'passos': num(NSTEPS_PROT), 'ns/dia': num(PERF_CG or 0)},
      compare={'aa': ('100 ps no curso atomístico', '', 'O curso atomístico simula 100 ps em 50 000 passos de 2 fs.'),
               'ua': ('dt = 2 fs', '', 'Custo por passo próximo do all-atom em sistemas solvatados: a água é a mesma.'),
               'cg': (f'{num(T_PROT)} ns em {num(NSTEPS_PROT)} passos', f'{num(PERF_CG or 0)} ns/dia', 'Compare com as estimativas teóricas mais abaixo.')})
''')

md('''
## Análise A: RMSD, RMSF, raio de giração e mapa de contatos

As mesmas análises do curso atomístico, calculadas sobre os beads BB (o equivalente CG do Cα). DSSP e ligações de hidrogênio não se aplicam: não há átomos de backbone nem hidrogênios, e a estrutura secundária está fixada na topologia.

**Na aula:** com a rede elástica, a proteína fica num poço único em torno da estrutura inicial (painel da direita).

![Funil de energia](https://raw.githubusercontent.com/madsondeluna/curso-dinamica-molecular-202x/main/imagens/fund_funil.png)
''')

code('''
import os
ARGS = '-s md_prot.tpr -f md_prot_fit.xtc -n index_prot.ndx'
gmx_cached(['rmsd_cg.xvg'],   'RMSD',     f'printf "BB\\\\nBB\\\\n" | gmx rms {ARGS} -o rmsd_cg.xvg')
gmx_cached(['rmsf_cg.xvg'],   'RMSF',     f'echo BB | gmx rmsf {ARGS} -o rmsf_cg.xvg -res')
gmx_cached(['gyrate_cg.xvg'], 'Rg',       f'echo Protein | gmx gyrate {ARGS} -o gyrate_cg.xvg')
gmx_cached(['mdmat_cg.xpm'],  'contatos', f'echo Protein | gmx mdmat {ARGS} -mean mdmat_cg.xpm')

rmsd, rmsf, rg = read_xvg('rmsd_cg.xvg'), read_xvg('rmsf_cg.xvg'), read_xvg('gyrate_cg.xvg')
fig, ax = plt.subplots(figsize=(8, 3.5))
ax.plot(rmsd[:, 0] / 1000, rmsd[:, 1], color=C1, lw=0.9)
finish(ax, title='RMSD dos beads de backbone', ylabel='RMSD', subtitle='nm', xlabel='Tempo (ns)')
plt.show()
''')

md('''
**O que se observa:** o RMSD sobe no primeiro nanossegundo e passa a oscilar em torno de 0,2 nm.

**Significado:** a proteína relaxa da forma do cristal para a forma em solução e depois só flutua. Com a rede elástica, um RMSD baixo é imposto pelo modelo: ele não prova que a proteína é estável, apenas que as molas a seguram.
''')

code('''
fig, ax = plt.subplots(figsize=(8, 3.5))
ax.plot(rmsf[:, 0], rmsf[:, 1], color=C1, lw=1.1)
finish(ax, title='RMSF por resíduo', ylabel='RMSF', subtitle='nm', xlabel='Resíduo')
plt.show()
''')

md('''
**O que se observa:** a maior parte da cadeia flutua menos de 0,1 nm; há picos perto do resíduo 9, nos resíduos 46-47 e na cauda C-terminal (72-76).

**Significado:** flexibilidade é função. A Gly76 da cauda C-terminal é o ponto em que a ubiquitina se liga covalentemente a lisinas de outras proteínas para marcá-las; a cauda precisa ser móvel para isso. A alça do resíduo 8 fica junto da Ile44 e da Val70, a região hidrofóbica que proteínas parceiras reconhecem.
''')

code('''
fig, ax = plt.subplots(figsize=(8, 3.5))
ax.plot(rg[:, 0] / 1000, rg[:, 1], color=C1, lw=0.9)
finish(ax, title='Raio de giração', ylabel='Rg', subtitle='nm', xlabel='Tempo (ns)')
plt.show()
''')

md('''
**O que se observa:** o raio de giração fica praticamente constante.

**Significado:** a proteína continua compacta, sem expandir nem desenovelar. Uma proteína que desenovela mostra o Rg crescendo.
''')

code('''
fig, ax = plt.subplots(figsize=(6.2, 5.4))
im = ax.imshow(xpm_to_distance_matrix('mdmat_cg.xpm'), cmap='viridis_r', origin='lower',
               extent=[1, N_RES, 1, N_RES], vmax=1.5)
finish(ax, title='Distância mínima média entre resíduos', xlabel='Resíduo', ylabel='Resíduo')
plt.colorbar(im, ax=ax, label='nm')
plt.show()
''')

md('''
**O que se observa:** a diagonal amarela são resíduos vizinhos na sequência; as manchas claras longe da diagonal são resíduos distantes na sequência mas próximos no espaço.

**Significado:** é a impressão digital do enovelamento. As faixas perpendiculares à diagonal são fitas β pareadas na folha; os blocos entre a hélice (resíduos 23-34) e a folha mostram a hélice apoiada sobre ela.
''')

md('''
Os gráficos resumem os 15 ns. Abaixo, a estrutura no último quadro e os números principais.
''')

code('''
top_res = rmsf[np.argsort(rmsf[:, 1])[-4:][::-1], 0].astype(int)
stage('A', 8, 'Análise', structure='md_prot.gro',
      text=['A rede elástica mantém o RMSD baixo e estável: ele não mede enovelamento, porque as molas impedem o desenovelamento.',
            'O RMSF continua informativo. A cauda C-terminal tem poucas molas e é a região mais móvel, como na simulação atomística.'],
      legend=['BB', 'SC', 'NA', 'CL'],
      facts={'RMSD médio (nm)': num(rmsd[len(rmsd)//2:, 1].mean(), 3), 'Rg médio (nm)': num(rg[:, 1].mean(), 3),
             'resíduos mais móveis': ', '.join(map(str, top_res)), 'RMSF máximo (nm)': num(rmsf[:, 1].max(), 2)},
      compare={'aa': ('todas as análises', '', 'RMSD, RMSF, Rg, contatos, DSSP ao longo do tempo e ligações de hidrogênio.'),
               'ua': ('quase todas', '', 'Ligações de hidrogênio continuam possíveis, porque os H polares são explícitos.'),
               'cg': ('RMSD, RMSF, Rg, contatos', '', 'Sem DSSP dinâmico nem ligações de hidrogênio: a estrutura secundária é parâmetro, não resultado.')})

save_summary(protein={'rmsd_t': rmsd[:, 0].round(1).tolist(), 'rmsd': rmsd[:, 1].round(4).tolist(),
                      'rg_t': rg[:, 0].round(1).tolist(), 'rg': rg[:, 1].round(4).tolist(),
                      'rmsf_res': rmsf[:, 0].astype(int).tolist(), 'rmsf': rmsf[:, 1].round(4).tolist()})
''')

md('''
## Estimativas teóricas de custo nas três resoluções

Os números saem do sistema CG desta réplica e de contagens. O custo de um passo de MD é dominado pelas interações entre pares dentro do cutoff, então cresce com

$$\\text{custo por ns} \\propto N \\times n_{viz} \\times \\frac{1}{\\Delta t}, \\qquad n_{viz} = \\rho \\cdot \\tfrac{4}{3}\\pi r_c^3$$

em que $N$ é o número de partículas, $\\rho$ a densidade de partículas, $r_c$ o cutoff e $\\Delta t$ o passo. All-atom e united-atom usam o mesmo volume da caixa CG, com cada bead W trocado por 4 moléculas de água de 3 sítios.

O fator é um **limite superior** do ganho do Martini: o custo real também tem PME, comunicação e partes que não escalam com $n_{viz}$, e sistemas CG pequenos não ocupam a CPU inteira. Ganhos medidos ficam abaixo deste número, mas na mesma direção.
''')

code('''
V_BOX = box_edge('md_prot.gro') ** 3
N_W_FINAL = count_res('md_prot.gro', 'W')
N_IONS = 2 * N_PAIRS
res3 = {
    # nome: (particulas, cutoff nm, passo fs)
    'All-atom':    (N_AA + 12 * N_W_FINAL + N_IONS, 1.2, 2),
    'United-atom': (N_UA + 12 * N_W_FINAL + N_IONS, 1.4, 2),
    'Martini 3':   (count_particles('md_prot.gro'), 1.1, 20),
}
cost = {}
rows = []
for name, (n, rc, dt) in res3.items():
    rho = n / V_BOX
    nviz = rho * 4 / 3 * np.pi * rc**3
    cost[name] = n * nviz / dt
    rows.append([name, num(n), num(rho, 1), num(nviz), num(dt), num(1e6 / dt)])
for r in rows:
    r.append(num(cost[r[0]] / cost['Martini 3']))

PERF_CG = performance('md_prot.log')
targets = [('15 ns', 15), ('100 ns', 100), ('1 µs', 1000), ('10 µs', 10000)]

def wall(ns, ns_day):
    h = ns / ns_day * 24
    if h < 1: return f'{num(h * 60)} min'
    if h < 48: return f'{num(h, 1)} h'
    if h < 24 * 365: return f'{num(h / 24)} dias'
    return f'{num(h / 24 / 365, 1)} anos'

from IPython.display import Markdown, display
head = ['Resolução', 'Partículas', 'ρ (nm⁻³)', 'Vizinhos no cutoff', 'Δt (fs)', 'Passos por ns', 'Custo por ns relativo ao Martini']
t1 = '| ' + ' | '.join(head) + ' |\\n|' + '---|' * len(head) + '\\n' + '\\n'.join('| ' + ' | '.join(r) + ' |' for r in rows)
head2 = ['Tempo simulado'] + list(res3) + ['Martini 3, tempo efetivo ÷ 4']
t2_rows = []
for label, ns in targets:
    row = [label]
    for name in res3:
        row.append(wall(ns * cost[name] / cost['Martini 3'], PERF_CG) if PERF_CG else '?')
    row.append(wall(ns / 4, PERF_CG) if PERF_CG else '?')
    t2_rows.append(row)
t2 = '| ' + ' | '.join(head2) + ' |\\n|' + '---|' * len(head2) + '\\n' + '\\n'.join('| ' + ' | '.join(r) + ' |' for r in t2_rows)
display(Markdown('**Contagens e fator de custo**\\n\\n' + t1 + '\\n\\n**Tempo de parede estimado** '
                 f'(Martini medido nesta réplica: {num(PERF_CG or 0)} ns/dia; all-atom e united-atom extrapolados pelo fator acima, limite superior)\\n\\n' + t2))
''')

md('''
**Como ler:** o united-atom tem quase o mesmo número de partículas do all-atom porque a água, que é a maior parte do sistema, continua com 3 sítios por molécula. O ganho do Martini vem de três fatores que se multiplicam: menos partículas (a água fica 12 vezes menor), menos vizinhos por partícula dentro do cutoff (a densidade de beads é cerca de 10 vezes menor) e passo 10 vezes maior. Sobre isso ainda vem o tempo efetivo: a paisagem mais lisa faz cada nanossegundo CG equivaler a vários nanossegundos de dinâmica atomística.
''')

# ---------------------------------------------------------------- Parte B
md('''
# Parte B: Automontagem de uma Bicamada Lipídica (DSPC)

**Nesta parte você vai:** espalhar 128 lipídios ao acaso numa caixa de água, deixar a simulação correr e ver a membrana se formar sozinha; depois medir sua espessura e a área de cada lipídio.

**Por que isso importa.** Toda célula é delimitada por uma bicamada lipídica, e ninguém a monta peça por peça: ela se forma porque as caudas apolares evitam a água. É o mesmo princípio que forma lipossomos, as vesículas usadas para entregar fármacos e as vacinas de mRNA em nanopartículas lipídicas.

Lipídios distribuídos ao acaso em água se organizam sozinhos em uma bicamada, guiados pelo efeito hidrofóbico. O sistema segue o tutorial oficial do Martini 3:

| Componente | Quantidade |
|------------|-----------|
| DSPC (fosfatidilcolina, 2 caudas de 18 carbonos) | 128 moléculas × 12 beads |
| W | 1 536 beads (12 por lipídio = 48 moléculas de água por lipídio) |
| Caixa inicial | 7,5 × 7,5 × 7,5 nm |
| Temperatura | 340 K, acima da transição gel-fluido do DSPC (~328 K) |

Na Parte A a estrutura estava dada e a topologia a mantinha. Aqui nada é imposto: a bicamada tem de surgir das interações entre beads.

Réplica 1 desta aula, em projeção lateral (cabeças polares em azul, caudas em cinza):

![Automontagem em quadros](https://raw.githubusercontent.com/madsondeluna/curso-dinamica-molecular-202x/main/imagens/fund_automontagem_quadros.png)
''')

code('''
import os
if not cached(['DSPC-em.gro'], 'estrutura de 1 DSPC'):
    # O tutorial usa a estrutura de DPPC: mesmo numero de beads e mesmos nomes do DSPC Martini 3
    os.system('wget -q https://cgmartini-library.s3.ca-central-1.amazonaws.com/1_Downloads/'
              'ff_parameters/martini2/lipidome/pc/dppc/DPPC-em.gro -O DSPC-em.gro')
    os.system("sed -i 's/DPPC/DSPC/' DSPC-em.gro")

stage('B', 1, 'Um lipídio DSPC', structure='DSPC-em.gro', box=False, height=380,
      text=['A cabeça polar tem colina (NC3, carga +1) e fosfato (PO4, carga −1). Dois beads de glicerol ligam a cabeça às caudas.',
            'Cada cauda de 18 carbonos vira 4 beads apolares. É essa diferença entre cabeça e caudas que move a automontagem.'],
      legend=['NC3', 'PO4', 'GL', 'C'],
      facts={'beads': 12, 'átomos pesados': 54, 'carga líquida': 0, 'caudas': '2 × C18:0'},
      compare={'aa': ('CHARMM36', '142 átomos', 'C₄₄H₈₈NO₈P com todos os 88 hidrogênios.'),
               'ua': ('GROMOS / Berger', '54 partículas', 'Um por átomo pesado; cada CH₂ e CH₃ das caudas vira uma partícula.'),
               'cg': ('Martini 3', '12 beads', '4,5 átomos pesados por bead em média.')})
''')

md('''
## Passo B1: Lipídios em posições aleatórias e solvatação

`gmx insert-molecules` coloca 128 cópias em posições e orientações aleatórias, rejeitando sobreposições. `gmx solvate -maxsol 1536` adiciona exatamente 1 536 beads W.
''')

code('''
import os
# Topologia fixa: 128 DSPC + 1536 W (reescrita sempre, e identica a cada execucao)
with open('dspc.top', 'w') as f:
    f.write("""#include "martini_v3.0.0.itp"
#include "martini_v3.0.0_ffbonded_v2.itp"
#include "martini_v3.0.0_phospholipids_PC_v2.itp"
#include "martini_v3.0.0_solvents_v1.itp"

[ system ]
DSPC bilayer self-assembly

[ molecules ]
DSPC   128
W      1536
""")
# Configuracao aleatoria: gerada uma unica vez e reaproveitada depois
if not cached(['mem_noW.gro'], 'lipidios aleatorios'):
    os.system(f'gmx insert-molecules -ci DSPC-em.gro -box 7.5 7.5 7.5 -nmol 128 -radius 0.21 -try 500 '
              f'-seed {tool_seed("insert-molecules")} -o mem_noW.gro -quiet 2>&1 | grep -i added')
if not cached(['mem_solv.gro'], 'lipidios + agua'):
    os.system('gmx solvate -cp mem_noW.gro -cs water.gro -o mem_solv.gro -maxsol 1536 -radius 0.21 -quiet')

stage('B', 2, 'Lipídios em posições aleatórias', structure='mem_noW.gro',
      text=['128 lipídios em posições e orientações sorteadas. Não há membrana, micela nem qualquer ordem: é o estado inicial.'],
      legend=['NC3', 'PO4', 'GL', 'C'],
      facts={'lipídios': 128, 'beads': num(count_particles('mem_noW.gro')), 'aresta da caixa (nm)': '7,5', 'tentativas por lipídio': 500})
''')

md('''
A mesma caixa com a água: os beads W preenchem os espaços entre os lipídios.
''')

code('''
stage('B', 3, 'Solvatação', structure='mem_solv.gro',
      text=['1 536 beads W, 12 por lipídio.'],
      legend=['NC3', 'PO4', 'GL', 'C', 'W'],
      facts={'beads W': num(1536), 'H₂O equivalentes': num(4 * 1536), 'água por lipídio': 48, 'partículas': num(count_particles('mem_solv.gro'))},
      compare={'aa': ('CHARMM36 + TIP3P', f'≈ {num(128 * 142 + 3 * 4 * 1536)} átomos', '18 176 átomos de lipídio mais 18 432 de água.'),
               'ua': ('GROMOS + SPC', f'≈ {num(128 * 54 + 3 * 4 * 1536)} partículas', 'As caudas encolhem muito, a água não.'),
               'cg': ('Martini 3 + W', f'{num(128 * 12 + 1536)} beads', 'O mesmo sistema com 12 vezes menos partículas que o all-atom.')})
''')

md('''
## Passo B2: Minimização e simulação de automontagem (15 ns)

O acoplamento de pressão é **isotrópico**: ainda não há bicamada, logo não há direção normal para tratar separadamente. Lipídios e água são acoplados ao termostato em grupos separados.

> **Duração:** o tutorial oficial usa 30 ns. Este notebook usa 15 ns; a seção de réplicas mostra quanto o tempo de formação varia entre sementes.
>
> **Tempo de execução estimado no Colab:** 45 a 75 minutos. Com o Drive montado, a simulação é retomada do último checkpoint se a sessão cair: basta executar a célula de novo.

**Na aula:** o relógio do Martini anda mais rápido que o físico, porque a paisagem é mais lisa.

![Relógios diferentes](https://raw.githubusercontent.com/madsondeluna/curso-dinamica-molecular-202x/main/imagens/fund_relogios.png)

Com o mesmo ruído térmico nas quatro paisagens, as travessias da barreira em 1 milhão de passos caem de 52 no coarse-grained para 32 no united-atom, 18 no all-atom e 11 na QM: cada camada de rugosidade desacelera o relógio. Na QM, cada passo também custa muito mais, porque resolve os elétrons.
''')

code('''
import os
run_md('mem_em', 'minimizacao membrana', 'minim.mdp', 'mem_solv.gro', 'dspc.top')
run_md('mem_md', 'automontagem (15 ns)', 'md_membrane.mdp', 'mem_em.gro', 'dspc.top')
if not cached(['mem_anim.pdb'], 'animacao membrana'):
    os.system('echo DSPC | gmx trjconv -s mem_md.tpr -f mem_md.xtc -o mem_anim.pdb -pbc mol -dt 500 -quiet > /dev/null 2>&1')
N_FRAMES = open('mem_anim.pdb').read().count('MODEL')
T_SIM = float(re.search(r'nsteps\\s*=\\s*(\\d+)', open('md_membrane.mdp').read()).group(1)) * 0.02 / 1000

PERF_MEM = performance('mem_md.log')
from replica_views import export_movie, single_view
movie, movie_water = export_movie('.')            # 0 a 15 ns, um quadro a cada 0,5 ns
stage('B', 5, 'Automontagem', figure_html=lambda: single_view(movie, 'membrane', animate=True, water_path=movie_water),
      text=['Filme completo de 0 a 15 ns, um quadro a cada 0,5 ns, com a água translúcida e a normal final da bicamada na vertical. A sequência típica é agregação em micelas pequenas, '
            'fusão em uma estrutura contínua e organização em duas camadas com as caudas voltadas para dentro.',
            'A bicamada pode surgir em qualquer plano da caixa e atravessar a borda periódica.'],
      legend=['NC3', 'PO4', 'GL', 'C', 'W'],
      facts={'tempo simulado (ns)': num(T_SIM), 'passo (fs)': 20, 'quadros': N_FRAMES, 'ns/dia': num(PERF_MEM or 0)},
      compare={'aa': ('CHARMM36', '', 'Possível, mas caro: dezenas de nanossegundos com ~27 mil átomos.'),
               'ua': ('GROMOS / Berger', '', 'A automontagem de DPPC em bicamada foi simulada em united-atom por Marrink et al. (2001).'),
               'cg': ('Martini 3', f'{num(PERF_MEM or 0)} ns/dia', 'Cabe numa aula: o tutorial oficial roda em minutos numa CPU de laboratório.')})
''')

md('''
## Análise B: Energia potencial e parâmetro de ordem

**Energia potencial:** cai em degraus a cada rearranjo (fusão de agregados, fechamento da bicamada) e estabiliza quando a estrutura final se forma.

**Parâmetro de ordem nemático $S$:** para cada lipídio, o vetor $\\hat{u}$ vai do bead PO4 ao ponto médio das extremidades das caudas (C4A, C4B). O tensor de ordem é

$$Q = \\left\\langle \\tfrac{3}{2}\\,\\hat{u}\\hat{u}^T - \\tfrac{1}{2} I \\right\\rangle$$

e $S$ é o seu maior autovalor. Com orientações aleatórias, $S$ fica perto de 0,2 (e não de 0: o maior autovalor de uma amostra finita de 128 vetores é sempre positivo); com todos os lipídios alinhados, $S \\to 1$. O autovetor correspondente é a normal da bicamada, encontrada sem supor em qual eixo ela se formou.

![Paisagem 2D](https://raw.githubusercontent.com/madsondeluna/curso-dinamica-molecular-202x/main/imagens/fund_paisagem_2d.png)

A QM, à esquerda, tem os mesmos vales do all-atom e uma rugosidade ainda mais fina, vinda dos elétrons (polarização, transferência de carga).
''')

code('''
import os
import mdtraj as md
gmx_cached(['mem_potential.xvg'], 'energia potencial',
           'echo Potential | gmx energy -f mem_md.edr -o mem_potential.xvg')
gmx_cached(['mem_whole.xtc'], 'trajetoria com moleculas inteiras',
           'echo System | gmx trjconv -s mem_md.tpr -f mem_md.xtc -o mem_whole.xtc -pbc mol')

traj = md.load('mem_whole.xtc', top='mem_em.gro')
top = traj.topology
po4 = top.select('name PO4')
if not cached(['mem_order.npz'], 'parametro de ordem'):
    c4a, c4b = top.select('name C4A'), top.select('name C4B')
    u = 0.5 * (traj.xyz[:, c4a] + traj.xyz[:, c4b]) - traj.xyz[:, po4]
    u /= np.linalg.norm(u, axis=2, keepdims=True)
    Q = 1.5 * np.einsum('fli,flj->fij', u, u) / u.shape[1] - 0.5 * np.eye(3)
    evals, evecs = np.linalg.eigh(Q)
    np.savez('mem_order.npz', t_ns=traj.time / 1000, S=evals[:, -1], normal=evecs[-1, :, -1])
o = np.load('mem_order.npz')
t_ns, S, normal = o['t_ns'], o['S'], o['normal']

pot = read_xvg('mem_potential.xvg')
fig, ax = plt.subplots(figsize=(8, 3.5))
ax.plot(pot[:, 0] / 1000, pot[:, 1], color=CM, lw=0.7)
ax.set_ylim(pot[pot[:, 0] > 200, 1].min() - 300, pot[pot[:, 0] > 200, 1].max() + 300)
finish(ax, title='Energia potencial durante a automontagem', ylabel='Energia', subtitle='kJ/mol', xlabel='Tempo (ns)')
plt.show()
''')

md('''
**O que se observa:** a energia cai rápido no começo e depois em degraus, até estabilizar.

**Significado:** cada degrau é um rearranjo: lipídios se agregando, agregados se fundindo, a bicamada se fechando. A energia cai porque as caudas apolares deixam de tocar a água e passam a tocar outras caudas. É o efeito hidrofóbico que monta as membranas das células.
''')

code('''
fig, ax = plt.subplots(figsize=(8, 3.5))
ax.plot(t_ns, S, color=CM, lw=1.1)
ax.set_ylim(0, 1)
finish(ax, title='Parâmetro de ordem S', ylabel='S', xlabel='Tempo (ns)')
plt.show()
''')

md('''
**O que se observa:** o S começa perto de 0,2 (orientações ao acaso) e, quando a bicamada se forma, sobe para cerca de 0,7.

**Significado:** os lipídios passam a apontar todos na mesma direção, perpendicular ao plano da membrana. O momento da subida é o tempo de formação da bicamada.
''')

code('''

above = np.where(S > 0.5)[0]
T_FORM = t_ns[above[0]] if len(above) else None
print(f'S inicial: {S[0]:.2f} | S final: {S[-1]:.2f}')
print('S > 0.5 a partir de', f'{T_FORM:.1f} ns' if T_FORM is not None else 'nunca (bicamada nao formada, ver nota abaixo)')
print(f'Normal da bicamada (ultimo frame): [{normal[0]:+.2f} {normal[1]:+.2f} {normal[2]:+.2f}]')
''')

md('''
## Análise B: que estrutura se formou? Percolação com condições periódicas

O parâmetro de ordem diz se os lipídios estão alinhados, mas não diz qual estrutura eles formaram. Numa caixa periódica, uma bicamada pode atravessar a borda em diagonal e um agregado pode parecer contínuo só porque encosta na própria imagem. A **dimensionalidade de percolação** resolve isso: o agregado de caudas é percorrido vizinho a vizinho, atravessando as bordas, e conta-se em quantas direções independentes ele se fecha sobre si mesmo através das imagens periódicas.

| Dimensão | Estrutura |
|----------|-----------|
| 0 | micela ou vesícula isolada |
| 1 | cilindro atravessando a caixa |
| 2 | bicamada (lamelar), em qualquer orientação |
| 3 | fase contínua, com conexões ou poros em todas as direções |

A bicamada conta como formada quando a dimensão é 2 **e** $S$ > 0,5.
''')

code('''
from collections import deque

def percolation(x, L, cut=0.6):
    \"\"\"Dimensionalidade do agregado de caudas respeitando as condicoes periodicas.\"\"\"
    d = x[:, None, :] - x[None, :, :]
    shift = np.round(d / L)                          # imagem periodica de cada par
    D = np.linalg.norm(d - shift * L, axis=-1)       # distancia minima entre imagens
    nb = [np.where((D[i] < cut) & (np.arange(len(x)) != i))[0] for i in range(len(x))]
    img = np.full((len(x), 3), np.nan); seen = np.zeros(len(x), bool); loops = []; n_agg = 0
    for s0 in range(len(x)):
        if seen[s0]: continue
        n_agg += 1; img[s0] = 0; seen[s0] = True; q = deque([s0])
        while q:
            i = q.popleft()
            for j in nb[i]:
                cand = img[i] + shift[i, j]          # em que copia da caixa j aparece, vista de i
                if not seen[j]:
                    seen[j] = True; img[j] = cand; q.append(j)
                elif np.any(cand != img[j]):
                    loops.append(cand - img[j])      # o agregado se fechou atraves da borda
    return (np.linalg.matrix_rank(np.array(loops)) if loops else 0), n_agg

STRUCT = {0: 'micela ou vesícula', 1: 'cilindro', 2: 'bicamada', 3: 'fase contínua'}
tails_idx = top.select('name =~ "C[1-4][AB]"')
if not cached(['mem_percolation.npz'], 'percolacao'):
    frames = np.linspace(0, traj.n_frames - 1, 16).astype(int)
    dims = []
    for f in frames:
        L = traj.unitcell_lengths[f]
        dims.append(percolation(traj.xyz[f][tails_idx] % L, L)[0])
    np.savez('mem_percolation.npz', t_ns=traj.time[frames] / 1000, dim=np.array(dims))
pc = np.load('mem_percolation.npz')
DIM = int(pc['dim'][-1])
FORMED = bool(DIM == 2 and S[-1] > 0.5)

fig, ax = plt.subplots(figsize=(8, 2.8))
ax.step(pc['t_ns'], pc['dim'], where='mid', color=CM, lw=1.6)
ax.set_yticks([0, 1, 2, 3]); ax.set_yticklabels([STRUCT[k] for k in range(4)]); ax.set_ylim(-0.4, 3.4)
finish(ax, title='Estrutura do agregado de caudas ao longo do tempo', xlabel='Tempo (ns)')
plt.show()
print(f'Estrutura final: {STRUCT[DIM]} (dimensao {DIM}) | S final {S[-1]:.2f} | bicamada formada: {"sim" if FORMED else "nao"}')
''')

md('''
## Passo B3: Equilibração semi-isotrópica da bicamada (5 ns)

Na automontagem a pressão é acoplada de forma isotrópica: a caixa cresce ou encolhe por igual nas três direções. Com uma bicamada formada, isso prende a área da membrana à forma cúbica da caixa. Na equilibração **semi-isotrópica**, o plano da membrana (xy) e a normal (z) se ajustam separadamente, e a área relaxa até o valor do campo de força.

O GROMACS assume a membrana no plano xy. Se a bicamada se formou com a normal em x ou y, o sistema é girado antes por uma permutação cíclica dos eixos, que é uma rotação própria da caixa cúbica.

> **Tempo de execução estimado no Colab:** 15 a 25 minutos. Só roda se a bicamada se formou.
''')

code('''
import os

def rotate_gro(src, dst, axis):
    \"\"\"Permutacao ciclica dos eixos de um .gro (posicoes, velocidades, caixa) que leva o eixo `axis` para z.\"\"\"
    order = {0: (1, 2, 0), 1: (2, 0, 1), 2: (0, 1, 2)}[axis]
    lines = open(src).read().splitlines()
    out = lines[:2]
    for l in lines[2:-1]:
        pos = [l[20:28], l[28:36], l[36:44]]
        vel = [l[44:52], l[52:60], l[60:68]] if len(l) >= 68 else None
        new = l[:20] + ''.join(pos[k] for k in order)
        if vel:
            new += ''.join(vel[k] for k in order)
        out.append(new)
    box = lines[-1].split()
    out.append('   ' + '   '.join(box[k] for k in order))
    open(dst, 'w').write('\\n'.join(out) + '\\n')

if FORMED:
    axis = int(np.argmax(np.abs(normal)))
    if not cached(['mem_rot.gro'], 'bicamada orientada no plano xy'):
        rotate_gro('mem_md.gro', 'mem_rot.gro', axis)
        print(f'Normal ao longo de {"xyz"[axis]}: sistema girado para a normal ficar em z.' if axis != 2 else 'Normal ja em z.')
    run_md('mem_eq', 'equilibracao semi-isotropica (5 ns)', 'md_membrane_eq.mdp', 'mem_rot.gro', 'dspc.top')
    gmx_cached(['mem_eq_whole.xtc'], 'trajetoria de equilibracao',
               'echo System | gmx trjconv -s mem_eq.tpr -f mem_eq.xtc -o mem_eq_whole.xtc -pbc mol')
    gmx_cached(['mem_eq_box.xvg'], 'dimensoes da caixa',
               'printf "Box-X\\\\nBox-Z\\\\n" | gmx energy -f mem_eq.edr -o mem_eq_box.xvg')
    box = read_xvg('mem_eq_box.xvg', ncol=3)
    fig, ax = plt.subplots(figsize=(8, 3.2))
    ax.plot(box[:, 0] / 1000, box[:, 1], color=CM, lw=1.2, label='x = y (plano da membrana)')
    ax.plot(box[:, 0] / 1000, box[:, 2], color=CM, lw=1.2, ls=(0, (6, 3)), label='z (normal)')
    finish(ax, title='Caixa durante a equilibração semi-isotrópica', subtitle='nm', xlabel='Tempo (ns)')
    plt.show()
else:
    print('Bicamada nao formada nesta replica: equilibracao semi-isotropica nao se aplica.')
''')

md('''
A caixa deixou de ser cúbica. Abaixo, a bicamada ao fim da equilibração semi-isotrópica, vista na orientação em que foi girada (normal em z).
''')

code('''
if FORMED:
    stage('B', 6, 'Equilibração semi-isotrópica', structure='mem_eq.gro',
          text=['A caixa deixou de ser cúbica: o plano da membrana e a normal se ajustam de forma independente.',
                'Espessura e área por lipídio são medidas nesta trajetória.'],
          legend=['NC3', 'PO4', 'GL', 'C', 'W'],
          facts={'plano xy final (nm)': num(box[-1, 1], 2), 'normal z final (nm)': num(box[-1, 2], 2),
                 'tempo (ns)': 5, 'acoplamento': 'semi-isotrópico'})
else:
    print('Bicamada nao formada nesta replica: equilibracao semi-isotropica nao se aplica.')
''')

md('''
## Análise B: Perfil de densidade, espessura e área por lipídio

As posições dos beads ao longo da normal (z) são centralizadas no meio da bicamada e convertidas em histogramas. A distância entre os picos de PO4 das duas monocamadas é a espessura $D_{PP}$. A área por lipídio é a área do plano xy dividida pelo número de lipídios de cada monocamada.
''')

code('''
if FORMED:
    teq = md.load('mem_eq_whole.xtc', top='mem_rot.gro')
    last = teq[teq.n_frames // 2:]                  # segunda metade da equilibracao
    tq = teq.topology
    po4e, tails = tq.select('name PO4'), tq.select('name =~ "C[1-4][AB]"')
    groups = {'NC3': tq.select('name NC3'), 'PO4': po4e, 'GL': tq.select('name =~ "GL[12]"'),
              'Caudas': tails, 'W': tq.select('resname W')}
    gcolor = {'NC3': BEADS['NC3'][1], 'PO4': BEADS['PO4'][1], 'GL': BEADS['GL'][1],
              'Caudas': '#8c8c8c', 'W': BEADS['W'][1]}            # mesmas cores dos cartoes
    Lz = last.unitcell_lengths[:, 2]
    z = last.xyz[:, :, 2]
    # Centro da bicamada por media circular das caudas (robusto a borda periodica)
    ang = 2 * np.pi * z[:, tails] / Lz[:, None]
    center = (np.arctan2(np.sin(ang).mean(1), np.cos(ang).mean(1)) / (2 * np.pi)) * Lz
    zc = (z - center[:, None] + Lz[:, None] / 2) % Lz[:, None] - Lz[:, None] / 2
    bins = np.linspace(-Lz.mean() / 2, Lz.mean() / 2, 81)
    mid = 0.5 * (bins[1:] + bins[:-1])
    fig, ax = plt.subplots(figsize=(8, 4))
    for name, idx in groups.items():
        h, _ = np.histogram(zc[:, idx].ravel(), bins=bins)
        ax.plot(mid, h / last.n_frames, label=name, color=gcolor[name], lw=1.6)
    finish(ax, title='Perfil de densidade ao longo da normal', subtitle='beads por fatia',
                    xlabel='Distância ao centro da bicamada (nm)')
    plt.show()
    h_po4, _ = np.histogram(zc[:, po4e].ravel(), bins=bins)
    D_PP = mid[mid > 0][np.argmax(h_po4[mid > 0])] - mid[mid < 0][np.argmax(h_po4[mid < 0])]
    upper = int((zc[-1, po4e] > 0).sum())
    lower = len(po4e) - upper
    APL = (last.unitcell_lengths[:, 0] * last.unitcell_lengths[:, 1]).mean() / (len(po4e) / 2)
    print(f'Espessura D_PP: {D_PP:.2f} nm | area por lipidio: {APL:.3f} nm2 | lipidios por monocamada: {upper} e {lower}')
else:
    D_PP = APL = None
    upper = lower = None
    print(f'Bicamada nao se formou nesta replica ({STRUCT[DIM]}): espessura e area nao se aplicam.')

''')

md('''
O perfil de densidade mostra as duas monocamadas. Abaixo, a estrutura final da bicamada usada nas medidas, com a água translúcida.
''')

code('''
stage('B', 7, 'Bicamada formada' if FORMED else 'Estado final: sem bicamada plana',
      structure='mem_eq.gro' if FORMED else 'mem_md.gro',
      text=(['As cabeças (NC3, PO4) ficam voltadas para a água e as caudas formam o núcleo apolar.',
             'Nenhuma destas propriedades foi imposta. Espessura e área por lipídio emergem dos parâmetros dos beads.']
            if FORMED else
            [f'Nesta réplica os lipídios formaram uma {STRUCT[DIM]} (dimensão de percolação {DIM}), e não uma bicamada: '
             f'as orientações não se alinham (S final = {num(S[-1], 2)}). O tutorial do Martini avisa que a formação no tempo simulado não é garantida.',
             'É um resultado legítimo, e é o motivo de rodar réplicas: compare com as outras duas na seção de réplicas.']),
      legend=['NC3', 'PO4', 'GL', 'C', 'W'],
      facts={'estrutura': STRUCT[DIM], 'S final': num(S[-1], 2),
             'formou em (ns)': num(T_FORM, 1) if FORMED and T_FORM is not None else 'não formou',
             'espessura D_PP (nm)': num(D_PP, 2) if FORMED else 'não se aplica',
             'área por lipídio (nm²)': num(APL, 3) if FORMED else 'não se aplica'},
      compare={'aa': ('CHARMM36', '', 'Referência para estrutura fina: ordem de cada C-H, hidratação das cabeças, ligações de hidrogênio.'),
               'ua': ('GROMOS / Berger', '', 'Ordem das caudas por segmento CH₂, a partir das posições dos átomos de carbono.'),
               'cg': ('Martini 3', '', 'Espessura, área, fases e automontagem. Não resolve a ordem de cada ligação C-H.')})

save_summary(membrane={'t_ns': t_ns.round(2).tolist(), 'S': S.round(3).tolist(),
                       'pot_t': (pot[::10, 0] / 1000).round(2).tolist(), 'pot': pot[::10, 1].round(1).tolist(),
                       't_form': None if T_FORM is None or not FORMED else round(float(T_FORM), 2),
                       'formed': bool(FORMED), 'dim': DIM,
                       'd_pp': round(float(D_PP), 3) if FORMED else None,
                       'apl': round(float(APL), 4) if FORMED else None,
                       'leaflets': [upper, lower] if FORMED else None})
''')

md('''
# Réplicas e sementes

Uma simulação é uma amostra de uma trajetória possível. Uma conclusão tirada de uma única corrida (a bicamada se forma aos 5,5 ns, o resíduo 9 é flexível) pode ser acaso daquela trajetória. Com **réplicas independentes**, cada resultado vira média ± desvio entre corridas, e um evento raro em uma só fica visível como raro.

Ter várias corridas não basta: elas precisam ser **independentes**, e quem decide isso são as sementes dos geradores de números aleatórios.

| Semente | Onde | O que sorteia | Se ficar fixa |
|---------|------|---------------|---------------|
| `gen_seed` | `.mdp`, com `gen_vel = yes` | velocidades iniciais (distribuição de Maxwell-Boltzmann) | todas as réplicas partem das mesmas velocidades |
| `ld_seed` | `.mdp` | ruído do termostato estocástico (v-rescale) e da dinâmica de Langevin | o mesmo ruído a cada passo |
| `-seed` do `insert-molecules` | linha de comando | posição e orientação de cada lipídio | todas as réplicas da Parte B partem da mesma caixa |
| `-seed` do `genion` | linha de comando | quais beads W viram íons | mesma distribuição inicial de íons |

Com a mesma `.tpr` e as mesmas sementes, duas corridas produzem a mesma trajetória, a menos de diferenças de arredondamento em paralelo, que só separam as trajetórias depois de algum tempo. São cópias, não réplicas. O erro clássico é copiar um `.mdp` com `gen_seed` fixo e rodar "três réplicas" que são a mesma corrida.

O inverso também é um problema: com `-1` o GROMACS sorteia uma semente nova, mas não a grava no `.tpr`. Sem registrar o valor, a corrida não pode ser reproduzida. Este notebook faz as duas coisas: sorteia (`gen_seed = -1`, `ld_seed = -1`, `-seed` aleatório nas ferramentas) e guarda cada valor usado em `seeds.json`.

**Onde a independência começa:** as velocidades são sorteadas na primeira etapa de dinâmica (aqui, o NPT da Parte A e a automontagem da Parte B). Etapas com `gen_vel = no` herdam as velocidades do checkpoint. Para uma réplica nova, refaça a partir dessa etapa, não só a produção.

**Como gerar uma réplica:** crie uma pasta por réplica (por exemplo `MyDrive/aula_martini_rep2`), troque `WORK_DIR` na segunda célula e execute o notebook. Cada pasta terá o próprio `seeds.json`. **Para reproduzir uma corrida exata:** copie os valores de `seeds.json` para `gen_seed` e `ld_seed` nos `.mdp` e para o `-seed` das ferramentas.
''')

code('''
# Sementes usadas nesta pasta de trabalho
import json
print(json.dumps(seeds(), indent=1))
''')

md('''
## Comparação entre as três réplicas

As três réplicas partem da mesma estrutura e diferem só nas sementes: velocidades iniciais, ruído do termostato, posição dos íons e, na Parte B, a posição de cada lipídio. A primeira figura mostra cada réplica (linha contínua, tracejada e pontilhada, na mesma cor); a segunda, a média com faixa de ± 1 desvio padrão. Na ubiquitina, a faixa cinza marca a α-hélice detectada pelo DSSP.

O que procurar: na proteína, as três curvas seguem a **mesma tendência** (o RMSD sobe e estabiliza, o RMSF tem os mesmos picos). Na membrana, o $S$ sai de ~0,2 e, quando a bicamada se forma, chega a ~0,7; o **momento** da formação muda de réplica para réplica, e uma réplica pode não chegar lá em 15 ns. Uma propriedade do sistema é a que se repete; um detalhe de uma corrida só é o que não se repete.

A célula usa o `summary.json` desta réplica e baixa do GitHub os resumos das outras.
''')

code('''
from replica_plots import load_summaries, plot_metric, plot_scalar
sums = load_summaries(BASE_DIR, REPO)
print('Replicas disponiveis:', sorted(sums))

names = {'t_form': ('Formacao da bicamada', 'ns'), 'd_pp': ('Espessura D_PP', 'nm'), 'apl': ('Area por lipidio', 'nm2')}
''')

md('''
### Proteína: RMSD
''')

code('''
plot_metric(sums, 'rmsd', 'replicas'); plt.show()
''')

md('''
**O que se observa:** as três réplicas sobem e oscilam na mesma faixa, entre 0,15 e 0,25 nm.

**Significado:** o comportamento da proteína é reprodutível. As diferenças entre as linhas são flutuação térmica, não comportamentos diferentes.
''')

code('''
plot_metric(sums, 'rmsd', 'mean'); plt.show()
''')

md('''
**O que se observa:** a faixa de desvio padrão é estreita.

**Significado:** quando as réplicas concordam, a média representa bem o sistema.

### Proteína: RMSF por resíduo
''')

code('''
plot_metric(sums, 'rmsf', 'replicas', ss=SS); plt.show()
''')

md('''
**O que se observa:** os picos caem nos mesmos resíduos nas três réplicas: perto do 9, em 46-47 e na cauda C-terminal. A faixa cinza é a α-hélice, uma das regiões mais rígidas.

**Significado:** a flexibilidade de cada região é propriedade da proteína, não de uma corrida.
''')

code('''
plot_metric(sums, 'rmsf', 'mean', ss=SS); plt.show()
''')

md('''
### Membrana: parâmetro de ordem
''')

code('''
plot_metric(sums, 'S', 'replicas'); plt.show()
''')

md('''
**O que se observa:** em duas réplicas, o S sobe de ~0,2 para ~0,7 entre 4 e 7 ns; na terceira, fica perto de 0,2 e só começa a subir no fim.

**Significado:** a formação da bicamada é um evento cinético. O destino final é o mesmo, mas o tempo para chegar lá depende da semente, e uma corrida pode não chegar no tempo simulado.
''')

code('''
plot_metric(sums, 'S', 'mean'); plt.show()
''')

md('''
**O que se observa:** depois de 6 ns a faixa de desvio fica larga, e a média para em ~0,5.

**Significado:** a média esconde que duas réplicas formaram a bicamada e uma não. Nenhuma réplica tem S = 0,5. Olhe sempre as réplicas antes de confiar na média.

### Membrana: energia potencial
''')

code('''
plot_metric(sums, 'pot', 'replicas'); plt.show()
''')

md('''
**O que se observa:** a energia cai nas três réplicas; a que não formou a bicamada termina um pouco acima das outras.

**Significado:** a bicamada é o estado de menor energia. A réplica sem bicamada parou num vale mais raso da paisagem energética.

### Membrana: valores finais
''')

code('''
fig, stat = plot_scalar(sums, 't_form'); plt.show()
''')

md('''
**Tempo de formação.** Com só duas réplicas que formaram a bicamada, o desvio padrão é uma estimativa grosseira: para um tempo confiável seriam necessárias muitas réplicas.
''')

code('''
fig, stat = plot_scalar(sums, 'd_pp'); plt.show()
''')

md('''
**Espessura.** Distância entre os planos de fosfato das duas monocamadas, da ordem de 4 nm, como nas membranas das células. As duas réplicas concordam bem.
''')

code('''
fig, stat = plot_scalar(sums, 'apl'); plt.show()
''')

md('''
**Área por lipídio.** Cerca de 0,60 nm²: o espaço que cada lipídio ocupa no plano da membrana. Quanto maior a área, mais fluida e permeável a membrana; colesterol e caudas saturadas diminuem a área, caudas insaturadas aumentam.
''')

code('''
# Resumo numerico da triplicata (so replicas em que a bicamada se formou)
for k, (label, unit) in names.items():
    v = np.array([x['membrane'][k] for x in sums.values() if x.get('membrane', {}).get(k) is not None])
    if len(v):
        print(f'{label}: {v.mean():.3f} +- {v.std(ddof=1) if len(v) > 1 else 0:.3f} {unit} (n = {len(v)})')
''')

md('''
## Réplicas lado a lado em 3D

Cada réplica tem o próprio visualizador, seguido do resultado dela. Cores e representação seguem a convenção do Martini no VMD: bastões nas ligações entre beads, esferas nas cabeças polares e água translúcida.
''')

code('''
import os
from replica_views import export_replica, single_view, SEM
from simple_view import _legend as card_legend
from IPython.display import Markdown, display

export_replica('.')                                   # arquivos de visualizacao desta replica
VIZ = os.path.join(BASE_DIR, 'summaries', 'viz')
os.makedirs(VIZ, exist_ok=True)

def viz_path(r, kind):
    ext = 'xyz' if kind.endswith('water') else 'pdb'          # agua em XYZ: arquivo menor
    own = os.path.join(BASE_DIR, f'rep{r}', f'viz_{kind}.{ext}')
    if os.path.exists(own) and os.path.getsize(own):
        return own
    remote = os.path.join(VIZ, f'rep{r}_{kind}.{ext}')
    if not os.path.exists(remote) or os.path.getsize(remote) == 0:
        os.system(f'wget -q {REPO}/results/viz/rep{r}_{kind}.{ext} -O {remote}')
    return remote

state = {r: ('bicamada' if s.get('membrane', {}).get('formed') else 'sem bicamada') for r, s in sums.items()}
LEG_PROT = [('Backbone (BB)', SEM['protein']), ('Cadeia lateral (SC)', SEM['protein_sc'])]
LEG_MEM = [('Colina NC3', SEM['head_cat']), ('Fosfato PO4', SEM['po4']), ('Glicerol GL1, GL2', SEM['glyc']),
           ('Caudas', SEM['tail']), ('Água W', SEM['water'])]
for r in (1, 2, 3):
    for kind in ('prot', 'prot_water', 'mem_final', 'mem_water', 'mem_anim', 'mem_anim_water'):
        viz_path(r, kind)
print('Arquivos de visualizacao das tres replicas prontos.')
''')

md('''
### Ubiquitina nas três réplicas

As três partem da mesma estrutura. Arraste para girar; a animação cobre os 15 ns. Observe que a forma geral se mantém: a rede elástica segura o enovelamento.
''')

code('''
for r in (1, 2, 3):
    m = sums.get(r, {}).get('protein', {})
    rmsd = np.array(m.get('rmsd', [0]))
    display(Markdown(f'**Réplica {r}**: RMSD médio na segunda metade {num(rmsd[len(rmsd)//2:].mean(), 3)} nm'))
    single_view(viz_path(r, 'prot'), 'protein', animate=True, water_path=viz_path(r, 'prot_water'))
card_legend(LEG_PROT)
''')

md('''
### Membrana: estado final das três réplicas

Bicamadas vistas de perfil, com a normal na vertical e a água translúcida. A mesma receita, com sementes diferentes: compare as três antes de ler o resultado de cada uma.
''')

code('''
for r in (1, 2, 3):
    m = sums.get(r, {}).get('membrane', {})
    if m.get('formed'):
        txt = f'bicamada formada aos {num(m["t_form"], 1)} ns; espessura {num(m["d_pp"], 2)} nm, área por lipídio {num(m["apl"], 3)} nm²'
    else:
        txt = f'sem bicamada em 15 ns (S final {num(m["S"][-1], 2) if m.get("S") else "?"})'
    display(Markdown(f'**Réplica {r}**: {txt}'))
    single_view(viz_path(r, 'mem_final'), 'membrane', water_path=viz_path(r, 'mem_water'))
card_legend(LEG_MEM)
''')

md('''
### Automontagem nas três réplicas

Um quadro por nanossegundo, com a água translúcida e a normal final na vertical. Em cada réplica, observe em que momento as cabeças se separam em duas camadas.
''')

code('''
for r in (1, 2, 3):
    display(Markdown(f'**Réplica {r}**: {state.get(r, "?")} ao fim de 15 ns'))
    single_view(viz_path(r, 'mem_anim'), 'membrane', animate=True, water_path=viz_path(r, 'mem_anim_water'))
card_legend(LEG_MEM)
''')

md('''
### Início e fim em projeção ortográfica

**Membrana.** Primeira coluna: o início, com os lipídios espalhados pela caixa. Segunda: o fim visto de lado. Terceira: o fim visto de cima. As réplicas 1 e 2 mostram duas camadas de cabeças (azul e vermelho) separadas pelo núcleo de caudas (cinza); a réplica 3 terminou numa estrutura conectada.

![Membrana: início e fim nas três réplicas](https://raw.githubusercontent.com/madsondeluna/curso-dinamica-molecular-202x/main/imagens/martini_vistas_membrana.png)

**Proteína.** As três réplicas partem da mesma estrutura e terminam com a mesma forma geral; o que muda é a posição da cauda C-terminal, a parte flexível. Primeira coluna: o início; segunda: o fim visto de frente; terceira: o fim girado 90°.

![Ubiquitina: início e fim nas três réplicas](https://raw.githubusercontent.com/madsondeluna/curso-dinamica-molecular-202x/main/imagens/martini_vistas_proteina.png)
''')

md('''
# Membranas reais: células normais, de câncer e bactérias

A bicamada desta aula tem um único lipídio. Membranas de verdade misturam dezenas, e a composição de cada lado é diferente. Essa diferença decide como moléculas e fármacos interagem com cada tipo de célula.

![Composição das membranas](https://raw.githubusercontent.com/madsondeluna/curso-dinamica-molecular-202x/main/imagens/fund_membranas_celulas.png)

| Membrana | Monocamada externa | Monocamada interna | Superfície externa |
|----------|--------------------|--------------------|--------------------|
| Célula humana normal | PC, esfingomielina, colesterol | PE, PS, PC, colesterol | neutra |
| Célula de câncer | PC, esfingomielina, colesterol e **PS exposta** | PE, PS, PC, colesterol | aniônica |
| Bactéria (membrana interna de *E. coli*) | PE, PG, cardiolipina | igual à externa | aniônica, sem colesterol |

**Célula normal.** A célula gasta energia para manter a PS, um lipídio com carga negativa, só na face interna. Por fora, a membrana é neutra e rica em colesterol (van Meer *et al.*, 2008).

**Célula de câncer.** Muitas células tumorais perdem essa assimetria e expõem PS na face externa. Células tumorais de melanoma e carcinoma mostraram de 3 a 7 vezes mais PS na face externa que queratinócitos normais (Utsugi *et al.*, 1991). A superfície fica negativa, e isso é explorado no desenvolvimento de peptídeos anticâncer e de moléculas que reconhecem PS.

**Bactéria.** A membrana interna de *E. coli* é formada principalmente por PE, com PG e cardiolipina, ambos negativos, e não tem colesterol (Sohlenkamp e Geiger, 2016). A composição varia entre espécies e condições de crescimento.

**Exemplo translacional: peptídeos antimicrobianos.** Peptídeos do sistema imune, como magainina e LL-37, são catiônicos (carga positiva). Eles são atraídos pelas superfícies negativas das bactérias e interagem pouco com a superfície neutra e rica em colesterol das células humanas. Essa é uma das bases da seletividade que se busca em novos antibióticos (Zasloff, 2002). Os mesmos princípios valem para peptídeos que atacam células de câncer.

**Como simular essas membranas.** Todos esses lipídios têm parâmetros no Martini 3 (repositório M3-Lipid-Parameters: PC, PE, PS, PG, SM, cardiolipina e colesterol). Trocar a bicamada desta aula por uma membrana bacteriana ou tumoral é mudar a composição de entrada; o restante do protocolo é o mesmo.
''')

md('''
### Domínios laterais: balsas lipídicas

A membrana também não é homogênea no plano. Esfingomielina, colesterol e lipídios de caudas saturadas se agrupam em domínios mais ordenados e mais espessos (fase líquido-ordenada, Lo), cercados por fosfolipídios insaturados em fase líquido-desordenada (Ld). Esses domínios, as balsas lipídicas, são dinâmicos, têm dimensões nanométricas e concentram proteínas de sinalização (Simons e Ikonen, 1997; Lingwood e Simons, 2010).

![Balsa lipídica](https://raw.githubusercontent.com/madsondeluna/curso-dinamica-molecular-202x/main/imagens/fund_membrana_dominios.png)

**Câncer.** Células tumorais costumam acumular colesterol na membrana. As balsas reúnem receptores e moléculas de sinalização ligadas a crescimento, adesão, migração e invasão, e alterar o colesterol das balsas é estudado como estratégia terapêutica (Vona *et al.*, 2021).

**No Martini.** Uma mistura de PC saturada, PC insaturada e colesterol se separa sozinha em fases Lo e Ld em simulações Martini, com propriedades próximas das medidas experimentais (Risselada e Marrink, 2008). É o mesmo tipo de automontagem da Parte B, com três componentes em vez de um.
''')

md('''
### Insaturação das caudas: inflamação e câncer

Uma ligação dupla *cis* dobra a cauda. Caudas saturadas, como as 18:0 do DSPC desta aula, ficam retas e empacotam densamente; caudas poli-insaturadas, como o araquidonato (20:4, quatro duplas ligações), dobram, ocupam mais área e deixam a membrana mais fluida. As duplas ligações também são o ponto em que a cauda reage com oxigênio.

![Caudas saturadas e poli-insaturadas](https://raw.githubusercontent.com/madsondeluna/curso-dinamica-molecular-202x/main/imagens/fund_insaturacao.png)

**Inflamação.** A fosfolipase A2 libera o araquidonato da posição *sn*-2 dos fosfolipídios. Ciclo-oxigenases e lipoxigenases o convertem em prostaglandinas e leucotrienos, mediadores da inflamação (Dennis e Norris, 2015). Anti-inflamatórios não esteroides, como o ibuprofeno, inibem as ciclo-oxigenases.

**Câncer.** Tumores com lipogênese ativa produzem os próprios ácidos graxos e têm fosfolipídios mais saturados. Com menos duplas ligações, a membrana sofre menos peroxidação, e as células ficam mais protegidas de radicais livres e de quimioterápicos como a doxorrubicina (Rysman *et al.*, 2010). No sentido oposto, PE com araquidonato ou adrenato oxidados dispara a ferroptose, uma forma de morte celular regulada que é alvo de pesquisa de novos fármacos (Kagan *et al.*, 2017).

**No Martini.** A insaturação aparece no tipo de bead da cauda: o DSPC (18:0) só tem beads C; o DOPC (18:1) tem um bead D no lugar da dupla ligação. Trocar DSPC por DOPC na Parte B é trocar uma cauda reta por uma cauda dobrada.
''')

md('''
### Proteínas de membrana: receptores e transportadores

De 20 a 30% dos genes da maioria dos genomas codificam proteínas de membrana (Krogh *et al.*, 2001). Elas atravessam o núcleo hidrofóbico da bicamada, de cerca de 3 nm, com hélices de resíduos hidrofóbicos.

![Proteínas de membrana](https://raw.githubusercontent.com/madsondeluna/curso-dinamica-molecular-202x/main/imagens/fund_proteinas_membrana.png)

**Receptores.** Os receptores acoplados à proteína G (GPCRs) têm sete hélices transmembrana; cerca de 34% dos fármacos aprovados pela FDA agem em GPCRs (Hauser *et al.*, 2017). Os receptores tirosina-quinase, como o EGFR, têm uma hélice por monômero e sinalizam quando dois monômeros se juntam; o EGFR é alvo de terapias contra câncer.

**Transportadores.** A glicoproteína P (P-gp, ABCB1) usa ATP para expulsar quimioterápicos da célula. Sua superexpressão é uma das causas de resistência a múltiplos fármacos em tumores (Gottesman *et al.*, 2002).

**No Martini.** Simulações de dez proteínas de membrana em uma membrana plasmática com mais de 60 tipos de lipídio mostraram que cada proteína enriquece ou depleta lipídios específicos ao seu redor (Corradi *et al.*, 2018). O martinize2 desta aula converte proteínas de membrana do mesmo jeito que converteu a ubiquitina; o passo extra é inseri-las na bicamada.
''')

md('''
### Peptídeos antimicrobianos em membranas bacterianas

Peptídeos antimicrobianos (AMPs) como magainina e LL-37 são curtos, catiônicos e anfipáticos: uma face concentra lisinas e argininas, a outra é hidrofóbica. A face positiva se liga ao PG e à cardiolipina da bactéria; a face hidrofóbica entra no núcleo da bicamada. Três modelos descrevem o que acontece depois (Brogden, 2005; Shai, 2002):

![Mecanismos dos AMPs](https://raw.githubusercontent.com/madsondeluna/curso-dinamica-molecular-202x/main/imagens/fund_amp_membrana.png)

- **Tapete**: os peptídeos cobrem a superfície e, acima de uma concentração limite, desfazem a bicamada como um detergente.
- **Poro toroidal**: os peptídeos entram na membrana e arrastam as cabeças dos lipídios, que passam a forrar o poro junto com eles.
- **Barril**: os peptídeos se alinham como as aduelas de um barril e formam sozinhos a parede de um canal.

Parte dos AMPs atravessa a membrana sem rompê-la e age em alvos internos, como a síntese de DNA, de proteínas ou da parede celular (Brogden, 2005). A seletividade vem da carga: a superfície bacteriana é aniônica, e a face externa das células humanas é neutra e rica em colesterol.

**No Martini.** A formação de um poro envolve vários peptídeos e centenas de lipídios, numa escala de tempo além do alcance usual do all-atom; por isso esses processos costumam ser estudados em coarse-grained, com o peptídeo convertido pelo martinize2 e a membrana bacteriana montada com PE, PG e cardiolipina.
''')

md('''
# Proteína × membrana: duas naturezas, um campo de força

| | Parte A: ubiquitina | Parte B: DSPC |
|---|---|---|
| O que mantém a estrutura | topologia: estrutura secundária fixa + rede elástica | nada imposto: só interações entre beads |
| O que a simulação mostra | flutuações em torno de um estado conhecido | surgimento de uma estrutura a partir do acaso |
| Observável principal | RMSD, RMSF | parâmetro de ordem, espessura, área por lipídio |
| Limitação do modelo | não descreve enovelamento nem mudanças de estrutura secundária | sem hidrogênios: sem ligações H explícitas entre cabeças |

O Martini 3 foi parametrizado principalmente a partir de energias livres de partição entre solventes. Por isso ele descreve bem processos governados por hidrofobicidade (automontagem, partição de moléculas em membranas, agregação) e precisa de restrições extras para manter estruturas definidas por ligações de hidrogênio, como as de uma proteína.

| Pergunta | Resolução indicada |
|----------|--------------------|
| Ligações de hidrogênio, mecanismo de uma enzima, ligante no sítio ativo | all-atom |
| Dinâmica de proteína em sistemas grandes, com custo menor | united-atom |
| Membranas, automontagem, agregação, sistemas de centenas de nanômetros | coarse-grained |
''')

md('''
## Referências

- Souza, P. C. T. *et al.* Martini 3: a general purpose force field for coarse-grained molecular dynamics. *Nat. Methods* **2021**, *18*, 382-388. https://doi.org/10.1038/s41592-021-01098-3
- Kroon, P. C. *et al.* Martinize2 and Vermouth provide a unified framework for molecular topology generation. *eLife* **2025**, *12*, RP90627. https://doi.org/10.7554/eLife.90627
- Pedersen, K. B. *et al.* The Martini 3 lipidome: expanded and refined parameters improve lipid phase behavior. *ACS Cent. Sci.* **2025**, *11*, 1598-1610. https://doi.org/10.1021/acscentsci.5c00755
- Marrink, S. J.; Risselada, H. J.; Yefimov, S.; Tieleman, D. P.; de Vries, A. H. The MARTINI force field: coarse grained model for biomolecular simulations. *J. Phys. Chem. B* **2007**, *111*, 7812-7824. https://doi.org/10.1021/jp071097f
- Marrink, S. J.; Lindahl, E.; Edholm, O.; Mark, A. E. Simulation of the spontaneous aggregation of phospholipids into bilayers. *J. Am. Chem. Soc.* **2001**, *123*, 8638-8639. https://doi.org/10.1021/ja0159618
- McGibbon, R. T. *et al.* MDTraj: a modern open library for the analysis of molecular dynamics trajectories. *Biophys. J.* **2015**, *109*, 1528-1532. https://doi.org/10.1016/j.bpj.2015.08.015
- van Meer, G.; Voelker, D. R.; Feigenson, G. W. Membrane lipids: where they are and how they behave. *Nat. Rev. Mol. Cell Biol.* **2008**, *9*, 112-124. https://doi.org/10.1038/nrm2330
- Utsugi, T.; Schroit, A. J.; Connor, J.; Bucana, C. D.; Fidler, I. J. Elevated expression of phosphatidylserine in the outer membrane leaflet of human tumor cells and recognition by activated human blood monocytes. *Cancer Res.* **1991**, *51*, 3062-3066. PMID 2032247
- Sohlenkamp, C.; Geiger, O. Bacterial membrane lipids: diversity in structures and pathways. *FEMS Microbiol. Rev.* **2016**, *40*, 133-159. https://doi.org/10.1093/femsre/fuv008
- Zasloff, M. Antimicrobial peptides of multicellular organisms. *Nature* **2002**, *415*, 389-395. https://doi.org/10.1038/415389a
- Simons, K.; Ikonen, E. Functional rafts in cell membranes. *Nature* **1997**, *387*, 569-572. https://doi.org/10.1038/42408
- Lingwood, D.; Simons, K. Lipid rafts as a membrane-organizing principle. *Science* **2010**, *327*, 46-50. https://doi.org/10.1126/science.1174621
- Vona, R.; Iessi, E.; Matarrese, P. Role of cholesterol and lipid rafts in cancer signaling: a promising therapeutic opportunity? *Front. Cell Dev. Biol.* **2021**, *9*, 622908. https://doi.org/10.3389/fcell.2021.622908
- Risselada, H. J.; Marrink, S. J. The molecular face of lipid rafts in model membranes. *Proc. Natl. Acad. Sci. USA* **2008**, *105*, 17367-17372. https://doi.org/10.1073/pnas.0807527105
- Dennis, E. A.; Norris, P. C. Eicosanoid storm in infection and inflammation. *Nat. Rev. Immunol.* **2015**, *15*, 511-523. https://doi.org/10.1038/nri3859
- Rysman, E. *et al.* De novo lipogenesis protects cancer cells from free radicals and chemotherapeutics by promoting membrane lipid saturation. *Cancer Res.* **2010**, *70*, 8117-8126. https://doi.org/10.1158/0008-5472.CAN-09-3871
- Kagan, V. E. *et al.* Oxidized arachidonic and adrenic PEs navigate cells to ferroptosis. *Nat. Chem. Biol.* **2017**, *13*, 81-90. https://doi.org/10.1038/nchembio.2238
- Krogh, A.; Larsson, B.; von Heijne, G.; Sonnhammer, E. L. Predicting transmembrane protein topology with a hidden Markov model: application to complete genomes. *J. Mol. Biol.* **2001**, *305*, 567-580. https://doi.org/10.1006/jmbi.2000.4315
- Hauser, A. S. *et al.* Trends in GPCR drug discovery: new agents, targets and indications. *Nat. Rev. Drug Discov.* **2017**, *16*, 829-842. https://doi.org/10.1038/nrd.2017.178
- Gottesman, M. M.; Fojo, T.; Bates, S. E. Multidrug resistance in cancer: role of ATP-dependent transporters. *Nat. Rev. Cancer* **2002**, *2*, 48-58. https://doi.org/10.1038/nrc706
- Corradi, V. *et al.* Lipid-protein interactions are unique fingerprints for membrane proteins. *ACS Cent. Sci.* **2018**, *4*, 709-717. https://doi.org/10.1021/acscentsci.8b00143
- Brogden, K. A. Antimicrobial peptides: pore formers or metabolic inhibitors in bacteria? *Nat. Rev. Microbiol.* **2005**, *3*, 238-250. https://doi.org/10.1038/nrmicro1098
- Shai, Y. Mode of action of membrane active antimicrobial peptides. *Biopolymers* **2002**, *66*, 236-248. https://doi.org/10.1002/bip.10260
- Bross, P. F. *et al.* Approval summary for bortezomib for injection in the treatment of multiple myeloma. *Clin. Cancer Res.* **2004**, *10*, 3954-3964. https://doi.org/10.1158/1078-0432.CCR-03-0781
- Tutorial Martini 3, *Lipid bilayers, Part I: Self-assembly*: https://cgmartini.nl/docs/tutorials/Martini3/LipidsI/
''')

nb = {"cells": cells,
      "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
                   "language_info": {"name": "python"}, "colab": {"provenance": []}},
      "nbformat": 4, "nbformat_minor": 5}
for i, c in enumerate(nb["cells"]):
    c["id"] = f"cell-{i:02d}"
    c["source"] = c["source"].splitlines(keepends=True)
json.dump(nb, open(sys.argv[1], "w"), indent=1, ensure_ascii=False)
print(f'{len(cells)} celulas -> {sys.argv[1]}')
