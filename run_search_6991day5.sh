#!/bin/bash
#$ -N search_6991day5
#$ -q UI,TELOMERE2
#$ -pe smp 56
#$ -j y
#$ -cwd

# Search 6991_day5_with_selection for ONE specific 278 bp circle, using
# population search mode (-pm -psm perfect -p <pattern>).
#
# 'perfect' mode reports, per telomeric read, the exact tandem copy count of
# the pattern, plus a copy-count histogram and a FASTA of every matching read.
# (The pattern is in TG orientation, matching how reads are normalized.)

set -uo pipefail

# ---------------------------------------------------------------- paths ----
SCRATCH="/nfsscratch/amalkova"
REPO="$SCRATCH/Find-dna-strings"
INPUT_DIR="$SCRATCH/circle_population_test"
OUT_DIR="$INPUT_DIR/population_results"
TMP_DIR="$OUT_DIR/_fasta_tmp"

INPUT_FILE="6991_day5_with_selection"
OUT_BASE="6991_day5_278bp"        # output files derive from this name
SEARCH_MODE="perfect"              # perfect | alignment | template_switching

PATTERN="TGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGGTGTGTGGGTGTGGTGTGGTGTGGTGTGTGTGGGTGTGGGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGGTGTGGGTGTGTGGTGTGTGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGGGTGTGGTGTGGGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGGTGTGTGGGTGTGGGTGTGTGTGGGTGTGGGTGTGG"

# -------------------------------------------------------------- conda -------
if command -v conda >/dev/null 2>&1; then
    eval "$(conda shell.bash hook)"
else
    for c in "$HOME/miniconda3" "$HOME/anaconda3" \
             "$SCRATCH/miniconda3" "$SCRATCH/anaconda3"; do
        if [ -f "$c/etc/profile.d/conda.sh" ]; then
            source "$c/etc/profile.d/conda.sh"
            break
        fi
    done
fi
conda activate consensus || { echo "ERROR: could not activate conda env 'consensus'"; exit 1; }

# ----------------------------------------------------- environment check ----
echo "--- environment ----------------------------------------------"
echo "  CONDA_DEFAULT_ENV : ${CONDA_DEFAULT_ENV:-<none>}"
echo "  python            : $(command -v python)"
python -c "import sys; print('  sys.executable    :', sys.executable)"
python -c "import parasail; print('  parasail          :', parasail.__version__)" \
    || { echo "ERROR: parasail not importable in '${CONDA_DEFAULT_ENV:-<none>}'"; exit 1; }
echo "--------------------------------------------------------------"

# ------------------------------------------------------------- settings -----
NSLOTS="${NSLOTS:-56}"
export MPLBACKEND=Agg

cd "$REPO" || { echo "ERROR: cannot cd to $REPO"; exit 1; }
mkdir -p "$OUT_DIR" "$TMP_DIR"

echo "Pattern length: ${#PATTERN} bp   mode: $SEARCH_MODE"

in="$INPUT_DIR/$INPUT_FILE"
fasta="$TMP_DIR/$INPUT_FILE.fasta"

# Convert FASTQ -> FASTA if needed.
first=$(head -c1 "$in")
if [ "$first" = "@" ]; then
    echo "Converting FASTQ -> FASTA"
    awk 'NR%4==1{print ">"substr($0,2)} NR%4==2{print}' "$in" > "$fasta"
    created_tmp=1
elif [ "$first" = ">" ]; then
    fasta="$in"
    created_tmp=0
else
    echo "ERROR: $INPUT_FILE does not start with '@' or '>'"; exit 1
fi

echo "Running population search mode"
python main.py "$fasta" -pm -psm "$SEARCH_MODE" -p "$PATTERN" -w "$NSLOTS" \
    -o "$OUT_DIR/${OUT_BASE}.txt"
status=$?

[ "$created_tmp" -eq 1 ] && rm -f "$fasta"

if [ $status -eq 0 ]; then
    echo "Done. See:"
    echo "  $OUT_DIR/${OUT_BASE}_search.txt          (summary + histogram)"
    echo "  $OUT_DIR/${OUT_BASE}_search_matches.fasta (reads containing the pattern)"
else
    echo "FAILED (exit $status)"
fi
