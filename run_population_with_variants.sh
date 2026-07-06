#!/bin/bash
#$ -N pop_variants
#$ -q UI,TELOMERE2
#$ -pe smp 56
#$ -j y
#$ -cwd

# Run population mode on every .fq read dump in circle_population_test/, using the
# updated finder that also reports the TOP 15 NON-VARIANT (distinct-family) patterns
# alongside the normal top-15-by-score table. Then run the variation annotator on
# each result so the normal top-15 table is also marked with which patterns are
# variations of a higher-ranked one.
#
# Per file you get, in circle_population_test/population_results/:
#   <name>_output_pattern.txt           normal top-15 + TOP 15 NON-VARIANT section
#   <name>_output_pattern_annotated.txt same, with inline "variation of Rank N" marks
#   <name>_output_telomeres.fasta       extracted telomeric reads
#   <name>.log                          per-file run log
#
# Files are processed one at a time; population mode uses all -pe slots internally.

set -uo pipefail

# ---------------------------------------------------------------- paths ----
SCRATCH="/nfsscratch/amalkova"
REPO="$SCRATCH/Find-dna-strings"
INPUT_DIR="$SCRATCH/circle_population_test"
OUT_DIR="$INPUT_DIR/population_results"

# pattern-finder options (these are the program defaults, set explicitly)
TOP_PATTERNS=15
VARIANT_THRESHOLD=0.90
CANDIDATE_POOL=200

# -------------------------------------------------------------- conda -------
if command -v conda >/dev/null 2>&1; then
    eval "$(conda shell.bash hook)"
else
    for c in "$HOME/miniconda3" "$HOME/anaconda3" \
             "$SCRATCH/miniconda3" "$SCRATCH/anaconda3"; do
        if [ -f "$c/etc/profile.d/conda.sh" ]; then
            source "$c/etc/profile.d/conda.sh"; break
        fi
    done
fi
conda activate consensus || { echo "ERROR: could not activate conda env 'consensus'"; exit 1; }

# ----------------------------------------------------- environment check ----
echo "--- environment ----------------------------------------------"
echo "  CONDA_DEFAULT_ENV : ${CONDA_DEFAULT_ENV:-<none>}"
echo "  python            : $(command -v python)"
python -c "import parasail; print('  parasail          :', parasail.__version__)" \
    || { echo "ERROR: parasail not importable in '${CONDA_DEFAULT_ENV:-<none>}'"; exit 1; }
echo "--------------------------------------------------------------"

# ------------------------------------------------------------- settings -----
NSLOTS="${NSLOTS:-56}"
export MPLBACKEND=Agg

cd "$REPO" || { echo "ERROR: cannot cd to $REPO"; exit 1; }
mkdir -p "$OUT_DIR"

shopt -s nullglob
FQ_FILES=( "$INPUT_DIR"/*.fq )
if [ "${#FQ_FILES[@]}" -eq 0 ]; then
    echo "No .fq files found in $INPUT_DIR"; exit 1
fi

echo "==============================================================="
echo "Population run (with non-variant patterns + annotation)"
echo "  inputs        : $INPUT_DIR/*.fq  (${#FQ_FILES[@]} files)"
echo "  outputs       : $OUT_DIR"
echo "  workers       : $NSLOTS"
echo "  top_patterns  : $TOP_PATTERNS   variant_threshold: $VARIANT_THRESHOLD   candidate_pool: $CANDIDATE_POOL"
echo "==============================================================="

for fq in "${FQ_FILES[@]}"; do
    base="$(basename "${fq%.fq}")"
    fa="$INPUT_DIR/$base.fa"
    out="$OUT_DIR/${base}_output.txt"
    pattern_file="$OUT_DIR/${base}_output_pattern.txt"

    echo
    echo ">>> [$base] $(date '+%Y-%m-%d %H:%M:%S')"

    # Ensure a FASTA exists (reuse if already converted, else make it once).
    if [ ! -s "$fa" ]; then
        echo "    converting FASTQ -> FASTA ($base.fa)"
        if ! awk 'NR%4==1{print ">"substr($0,2)} NR%4==2{print}' "$fq" > "$fa"; then
            echo "    FAILED conversion -- skipping"; rm -f "$fa"; continue
        fi
    else
        echo "    using existing $base.fa"
    fi

    echo "    running population mode"
    python main.py "$fa" -pm -w "$NSLOTS" \
        -tn "$TOP_PATTERNS" -vt "$VARIANT_THRESHOLD" -cp "$CANDIDATE_POOL" \
        -o "$out" > "$OUT_DIR/${base}.log" 2>&1
    status=$?

    if [ $status -ne 0 ]; then
        echo "    FAILED population mode (exit $status) -- see $OUT_DIR/${base}.log"
        continue
    fi

    echo "    annotating top-15 with variation marks"
    python annotate_pattern_variations.py "$pattern_file" \
        >> "$OUT_DIR/${base}.log" 2>&1 \
        && echo "    OK -> ${base}_output_pattern.txt (+ _annotated.txt)" \
        || echo "    WARNING: annotation step failed -- see $OUT_DIR/${base}.log"
done

echo
echo "All done. Results in: $OUT_DIR"
