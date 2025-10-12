#!/usr/bin/env bash
set -euo pipefail

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

python_result="$(
  python3 "${DIR}/Python/test_phylo.py" \
  | awk -F'\t' '$1=="python"{gsub(/^[[:space:]]+|[[:space:]]+$/,"",$2); print $2; exit}'
)"

codon_result="$(
  codon run --release "${DIR}/Codon/test_phylo_codon.py" \
  | awk -F'\t' '$1=="codon"{gsub(/^[[:space:]]+|[[:space:]]+$/,"",$2); print $2; exit}'
)"

echo -e "Language\tRuntime"

echo -e "---------------------------"

echo -e "Python\t\t${python_result}"

echo -e "Codon\t\t${codon_result}"
