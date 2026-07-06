#!/bin/bash
#$ -N find_loops_population
#$ -q UI,TELOMERE2
#$ -pe smp 56
#$ -j y
#$ -cwd

# Run Find-dna-strings in POPULATION mode on every read-dump file in
# circle_population_test/. Each file is processed one at a time; population
# mode parallelizes internally across all -pe slots (NSLOTS), so we do NOT
# fan out files in parallel (that would oversubscribe the cores).
#
# NOTE: the input files are FASTQ (they start with '@', SRA dumps) even though
# the .fasta extension was left off. Find-dna-strings reads FASTA, so each file
# is converted FASTQ->FASTA on the fly before analysis.

set -uo pipefail

# ---------------------------------------------------------------- paths ----
SCRATCH="/nfsscratch/amalkova"
REPO="$SCRATCH/Find-dna-strings"
INPUT_DIR="$SCRATCH/circle_population_test"
OUT_DIR="$INPUT_DIR/population_results"
TMP_DIR="$OUT_DIR/_fasta_tmp"

# -------------------------------------------------------------- conda -------
# Activate the consensus env. Falls back through common conda locations
# (including scratch) if `conda` isn't already on PATH.
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
# Fail fast (before converting any files) if we're not in the right env or a
# required package is missing.
echo "--- environment ----------------------------------------------"
echo "  CONDA_DEFAULT_ENV : ${CONDA_DEFAULT_ENV:-<none>}"
echo "  python            : $(command -v python)"
python -c "import sys; print('  sys.executable    :', sys.executable)"
if ! python -c "import parasail; print('  parasail          :', parasail.__version__)"; then
    echo "ERROR: 'parasail' is not importable in env '${CONDA_DEFAULT_ENV:-<none>}'."
    echo "       main.py imports it at startup (via AlignmentStrategy), so the run"
    echo "       cannot proceed. Install it into the env, e.g.:"
    echo "         conda activate consensus && pip install parasail"
    exit 1
fi
echo "--------------------------------------------------------------"

# ------------------------------------------------------------- settings -----
NSLOTS="${NSLOTS:-56}"        # set by SGE; default 56 for manual runs
export MPLBACKEND=Agg          # headless safety (pop mode skips graphs anyway)

cd "$REPO" || { echo "ERROR: cannot cd to $REPO"; exit 1; }
mkdir -p "$OUT_DIR" "$TMP_DIR"

# Select the 6991_* and 7172_* read dumps; skip the small note (.txt) files.
mapfile -t FILES < <(find "$INPUT_DIR" -maxdepth 1 -type f \
    \( -name '6991_*' -o -name '7172_*' \) \
    ! -name '*.txt' ! -name '*.sh' -printf '%f\n' | sort)

echo "==============================================================="
echo "Find-dna-strings population run"
echo "  repo      : $REPO"
echo "  inputs    : $INPUT_DIR"
echo "  outputs   : $OUT_DIR"
echo "  workers   : $NSLOTS"
echo "  files     : ${#FILES[@]}"
printf '    - %s\n' "${FILES[@]}"
echo "==============================================================="

for f in "${FILES[@]}"; do
    in="$INPUT_DIR/$f"
    echo
    echo ">>> [$f] $(date '+%Y-%m-%d %H:%M:%S')"

    # Detect format and convert FASTQ->FASTA if needed.
    first=$(head -c1 "$in")
    if [ "$first" = "@" ]; then
        fasta="$TMP_DIR/$f.fasta"
        echo "    converting FASTQ -> FASTA"
        awk 'NR%4==1{print ">"substr($0,2)} NR%4==2{print}' "$in" > "$fasta"
        created_tmp=1
    elif [ "$first" = ">" ]; then
        fasta="$in"          # already FASTA
        created_tmp=0
    else
        echo "    WARNING: $f does not start with '@' or '>' -- skipping"
        continue
    fi

    echo "    running population mode"
    python main.py "$fasta" -pm -w "$NSLOTS" \
        -o "$OUT_DIR/${f}_output.txt" \
        > "$OUT_DIR/${f}.log" 2>&1
    status=$?

    if [ $status -eq 0 ]; then
        echo "    OK -> $OUT_DIR/${f}_output_pattern.txt"
    else
        echo "    FAILED (exit $status) -- see $OUT_DIR/${f}.log"
    fi

    # Drop the temp FASTA to reclaim space (only if we created it).
    [ "$created_tmp" -eq 1 ] && rm -f "$fasta"
done

echo
echo "All done. Results in: $OUT_DIR"
