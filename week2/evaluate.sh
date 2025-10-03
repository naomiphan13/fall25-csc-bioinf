#!/usr/bin/env bash
set -euo pipefail

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo -e "Biocodon Tests on Python"
echo -e "-------------------------------------"
result=$(python3 "${DIR}/code/test_motifs.py")
n50=$(echo "$result"| awk '{print $2}')
echo -e "-------------------------------------"
echo -e "Biocodon Tests on Codon"
echo -e "-------------------------------------"
result=$(codon run -release "${DIR}/code/test.py")
n50=$(echo "$result"| awk '{print $2}')
