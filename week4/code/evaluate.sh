#!/usr/bin/env bash
set -euo pipefail

DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$DIR"

printf "Method\t\t\tLanguage\tRuntime\n"

python3 -u "${DIR}/main.py "$@""
codon run -release "${DIR}/main_codon.py "$@""
