#!/usr/bin/env bash
# Roda as tres replicas do MD_Martini.ipynb em sequencia, localmente, com o
# proprio notebook (REPLICA=N, USE_PRECOMPUTED=False). Saidas em runs/aula_martini/repN.
#
# Uso: PY=/caminho/do/python KERNEL=nome THREADS=6 tools/run_replicas.sh [1 2 3]
#   PY      python com jupyter, nbconvert, vermouth, mdtraj, py3Dmol, matplotlib
#   KERNEL  kernel jupyter desse python
#   THREADS threads OpenMP por simulacao (padrao 6)
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PY="${PY:?defina PY}"; KERNEL="${KERNEL:?defina KERNEL}"; THREADS="${THREADS:-6}"
BIN="$(dirname "$PY")"
REPS=("${@:-1 2 3}"); [ $# -eq 0 ] && REPS=(1 2 3)
RUNS="$ROOT/runs"
mkdir -p "$RUNS"

"$PY" "$ROOT/tools/build_martini_notebook.py" "$ROOT/MD_Martini.ipynb"

for n in "${REPS[@]}"; do
  d="$RUNS/aula_martini/rep$n"
  mkdir -p "$d/viz"
  # o GitHub pode estar atras do repositorio local: entradas copiadas daqui
  cp "$ROOT"/mdp_martini/*.mdp "$d/"
  cp -R "$ROOT"/viz/stage_viewer.py "$ROOT"/viz/stage.css "$ROOT"/viz/replica_plots.py "$ROOT"/viz/pure "$d/viz/"
  sed -e "s/REPLICA = 1              #/REPLICA = $n              #/" \
      -e "s/USE_PRECOMPUTED = True   #/USE_PRECOMPUTED = False  #/" \
      "$ROOT/MD_Martini.ipynb" > "$RUNS/rep$n.ipynb"
  echo "rep$n inicio $(date +%H:%M)"
  ( cd "$RUNS" && PATH="$BIN:$PATH" OMP_NUM_THREADS="$THREADS" "$BIN/jupyter" nbconvert --to notebook \
      --execute "rep$n.ipynb" --output "rep${n}_executed.ipynb" \
      --ExecutePreprocessor.kernel_name="$KERNEL" --ExecutePreprocessor.timeout=10800 > "rep$n.log" 2>&1 ) \
    && echo "rep$n ok $(date +%H:%M)" || { echo "rep$n FALHOU, ver runs/rep$n.log"; exit 1; }
done
