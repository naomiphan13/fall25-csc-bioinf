#!/usr/bin/env bash
set -euo pipefail

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Table header (as in your screenshot)
printf "Method\t\t\tLanguage\tRuntime\n"

# Unbuffered Python so lines flush immediately
python3 -u main.py "$@"
codon run -release main_codon.py "$@"