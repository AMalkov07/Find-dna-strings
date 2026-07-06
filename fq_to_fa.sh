#!/bin/bash
#$ -N fq_to_fa
#$ -q UI,TELOMERE2
#$ -pe smp 56
#$ -j y
#$ -cwd

# Create a .fa (FASTA) version of every .fq (FASTQ) file in circle_population_test/.
# Pure awk conversion (no conda/tools needed). Files are processed in parallel.
#
# For each <name>.fq -> <name>.fa in the same directory:
#   - keeps the read id (text after '@' on the header line)
#   - drops the '+' separator and quality line
# Assumes standard 4-line FASTQ records with single-line sequences (SRA dumps).

set -uo pipefail

INPUT_DIR="/nfsscratch/amalkova/circle_population_test"
NSLOTS="${NSLOTS:-56}"

cd "$INPUT_DIR" || { echo "ERROR: cannot cd to $INPUT_DIR"; exit 1; }

shopt -s nullglob
FQ_FILES=( *.fq )
if [ "${#FQ_FILES[@]}" -eq 0 ]; then
    echo "No .fq files found in $INPUT_DIR"; exit 1
fi

echo "==============================================================="
echo "FASTQ -> FASTA conversion"
echo "  dir          : $INPUT_DIR"
echo "  files        : ${#FQ_FILES[@]}"
echo "  max parallel : $NSLOTS"
echo "==============================================================="

max_jobs="$NSLOTS"
for in in "${FQ_FILES[@]}"; do
    out="${in%.fq}.fa"
    if [ -s "$out" ]; then
        echo "skip (already exists): $out"
        continue
    fi
    (
        echo "converting $in -> $out"
        if awk 'NR%4==1{print ">"substr($0,2)} NR%4==2{print}' "$in" > "$out"; then
            echo "done: $out ($(wc -l < "$out") lines)"
        else
            echo "FAILED: $in"; rm -f "$out"
        fi
    ) &
    # throttle to max_jobs concurrent conversions
    while [ "$(jobs -r | wc -l)" -ge "$max_jobs" ]; do wait -n; done
done
wait

echo
echo "All conversions complete. FASTA files written alongside the .fq files in:"
echo "  $INPUT_DIR"
ls -1 "$INPUT_DIR"/*.fa
