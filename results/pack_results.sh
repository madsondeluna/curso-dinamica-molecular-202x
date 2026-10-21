#!/usr/bin/env bash
# Empacota as replicas simuladas localmente para o notebook baixar do GitHub.
# Os arquivos de texto sao anonimizados numa copia: caminhos locais, nome da
# maquina e usuario viram valores genericos.
# Uso: results/pack_results.sh <BASE_DIR com rep1, rep2, rep3>
set -euo pipefail

BASE="$(cd "${1:?uso: pack_results.sh <BASE_DIR>}" && pwd)"
OUT="$(cd "$(dirname "$0")" && pwd)"
HOST="$(hostname)"; HOST_S="$(hostname -s)"; USER_NAME="$(id -un)"
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT

for rep in "$BASE"/rep*; do
  n="${rep##*/rep}"
  [ -f "$rep/summary.json" ] || { echo "rep$n sem summary.json, pulando" >&2; continue; }
  dst="$TMP/rep$n"; mkdir -p "$dst"
  # Fica de fora o que o notebook baixa ou regenera em segundos
  rsync -a --exclude='viz' --exclude='*.mdp' --exclude='martini_v3.0.0*.itp' \
    --exclude='water.gro' --exclude='*.trr' --exclude='md_prot_pbc.xtc' \
    --exclude='rep.tar.gz' \
    --exclude='._*' --exclude='.DS_Store' --exclude='*.cpt' "$rep/" "$dst/"
  grep -rIlE "/private/tmp/|/Volumes/|/Users/|$HOST_S|$USER_NAME" "$dst" | while read -r f; do
    sed -E -i '' \
      -e "s#/private/tmp/[^ ]*/venv/bin/##g" \
      -e "s#/private/tmp/[^ ]*/rep$n#/content/aula_martini/rep$n#g" \
      -e "s#$BASE#/content/aula_martini#g" \
      -e "s#/Volumes/[^ ]*#/content#g" -e "s#/Users/[^ ]*#/content#g" \
      -e "s#$HOST#host#g" -e "s#$HOST_S#host#g" -e "s#$USER_NAME#user#g" "$f"
  done
  left=$(grep -rIlE "/private/tmp/|/Volumes/|/Users/|$HOST_S|$USER_NAME" "$dst" || true)
  [ -z "$left" ] || { echo "rep$n: ainda ha dados locais em: $left" >&2; exit 1; }
  tar -czf "$OUT/rep$n.tar.gz" -C "$dst" .
  cp "$dst/summary.json" "$OUT/summary_rep$n.json"
  echo "rep$n: $(du -h "$OUT/rep$n.tar.gz" | cut -f1)"
done
