#!/bin/bash
#$ -N refresh_circles
#$ -q UI,TELOMERE2
#$ -pe smp 56
#$ -j y
#$ -cwd

# Refresh all population outputs to the current code, WITHOUT re-extracting
# telomeres from the raw reads (the extraction is unchanged, so the existing
# *_output_telomeres.fasta are already current). Steps:
#   1) re-analyze each telomeres FASTA  -> *_output_pattern.txt (+ _pattern_reads.txt)
#        new scoring, length-aware consolidation (vt=0.93), CIRCLE SUMMARY verdict
#   2) annotate variation marks         -> *_output_pattern_annotated.txt
#   3) classify family reads (98% id)   -> *_output_family_reads.txt

set -uo pipefail

SCRATCH="/nfsscratch/amalkova"
REPO="$SCRATCH/Find-dna-strings"
RESULTS="$SCRATCH/circle_population_test/population_results"

# --- conda (consensus has parasail+biopython) ---
if command -v conda >/dev/null 2>&1; then
    eval "$(conda shell.bash hook)"
else
    for c in "$HOME/miniconda3" "$HOME/anaconda3" "$SCRATCH/miniconda3" "$SCRATCH/anaconda3"; do
        [ -f "$c/etc/profile.d/conda.sh" ] && source "$c/etc/profile.d/conda.sh" && break
    done
fi
conda activate consensus || { echo "ERROR: could not activate conda env 'consensus'"; exit 1; }

# edlib is needed by the read classifier (pure-stdlib elsewhere); install if missing
python -c "import edlib" 2>/dev/null || pip install --quiet edlib

export NSLOTS="${NSLOTS:-56}"
cd "$REPO" || { echo "ERROR: cannot cd to $REPO"; exit 1; }

echo "=============================================================="
echo "Refreshing population outputs from existing telomeres"
echo "  repo    : $REPO"
echo "  results : $RESULTS"
echo "  workers : $NSLOTS"
echo "=============================================================="

echo; echo ">>> 1) regenerate pattern files (from telomeres, no re-extraction)"
python refresh_from_telomeres.py "$RESULTS"

echo; echo ">>> 2) annotate variation marks"
python annotate_pattern_variations.py "$RESULTS"/*_output_pattern.txt

echo; echo ">>> 3) classify family reads (circularized vs single-instance, 98% identity)"
python classify_family_reads.py "$RESULTS"/*_output_pattern.txt

echo; echo "All outputs refreshed in $RESULTS"
