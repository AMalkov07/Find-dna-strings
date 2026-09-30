# population_followup

Scripts for following circles found by population mode (`main.py -pm`) across a set of samples,
e.g. a timecourse or a panel of survivors. They reuse the program's own reader, population search,
template-switching analysis and graph code; nothing here changes how circles are detected.

Run from the Find-dna-strings root with it on `PYTHONPATH` (the scripts import `utils`, `data_io`,
`analysis`), e.g. `PYTHONPATH=. python population_followup/<script>.py ...`.
Population mode is memory-hungry on long-telomere (Type II) samples (150–230 GB seen), so run the
heavy steps on a cluster node, not a laptop.

## Workflow

A strain folder `W/` holds one folder per sample and a `samples.txt` (one sample name per line,
in the order you want them plotted).

1. **Population mode on each sample** (program itself):
   ```
   python main.py W/S/S.fasta -pm -mtl 200 -w N -o W/S/S_output.txt
   python annotate_pattern_variations.py W/S/S_output_pattern.txt
   python classify_family_reads.py      W/S/S_output_pattern.txt
   ```
2. **Collect every flagged candidate** (families >= 90 bp marked `<- circle`) from all samples:
   ```
   python population_followup/build_units.py W W/samples.txt W/followup_all
   ```
   -> `units.tsv` (name, sample, family, sequence) and `units_info.tsv` (adds type: TG telomeric vs
   simple repeat, and the program's read count).
3. **Search every sample for every candidate** (exact matches, no length gate), one pass per sample:
   ```
   python population_followup/search_units.py W/S/S.fasta W/followup_all S N
   ```
   -> `<S>_<unit>_perfect_search.txt` (counts + copy histogram) and `_perfect_search_matches.fasta`.
4. **Plots**
   - Heatmap of all candidates x all samples (groups units that match >= 97% as one circle):
     `python population_followup/plot_all_candidates.py W/followup_all W/samples.txt STRAIN`
     -> `STRAIN_all<N>_circles_heatmap.png`, `STRAIN_all<N>_counts.tsv`
   - Timecourse of chosen units (log scale, read counts on bars):
     `python population_followup/plot_units_timecourse.py W/followup_all W/samples.txt STRAIN OUT_PREFIX UNIT=LABEL ...`
   - Survivor panel, one bar per survivor for its own main circle (98%-identity family counts):
     `python population_followup/plot_survivor_main_circles.py W/followup_all W/samples.txt STRAIN`
   - Every read carrying one unit, reference-mode template-switching style, one row per read:
     ```
     python population_followup/plot_reads_template_switching.py --units-tsv W/followup_all/units.tsv \
         --unit NAME --out OUT_PREFIX day6=W/followup_all/day6_NAME_perfect_search_matches.fasta \
         IT134=W/IT134/IT134_output_family_reads.txt:2
     ```

## Things to know when reading the results

- **Two kinds of per-read count.** Step 1's `_output_family_reads.txt` lists reads carrying each
  family at >= 98% identity (families >= 100 bp only, that sample's telomeres >= 200 bp only).
  Step 3's search counts **exact** copies in all telomeres. Exact matching undercounts long units
  (one nanopore error breaks a copy), so a unit can have thousands of family reads and ~0 exact hits.
- **Simple-repeat "circles".** `CTTT`, `AGGG`, `GGT`, `GTT` repeats pass the 75% G+T telomere test
  and show up as perfect tandem arrays of any length. `build_units.py` labels them
  (`simple repeat (...)`) and the heatmap draws them in red; they are not telomeric circles.
- **Candidates need >= 2 reads.** Population mode only reports a pattern present >= 2 times in
  >= 2 reads (`min_reads=2` is fixed in `find_best_population_patterns`), so a circle seen in a
  single read is never a candidate.
- **Rows grouped as one circle** (>= 97% rotation-aware containment) are usually dimers or variants
  of the same unit; a ~95% match to a circle is at the level of ordinary telomere background.
