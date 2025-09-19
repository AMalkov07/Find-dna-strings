from collections import defaultdict
from typing import List, Dict, Tuple, Optional, Set
import parasail
import re

from utils.data_structures import SeedExtendSettings, SeedExtendCluster, ImperfectAlignmentEvent, Config
from parasail.bindings_v2 import Result

class SeedAndExtend:
    def __init__(self, mutagenic_zone: str, pattern: str, seed_extend_settings: SeedExtendSettings, config: Config):
        self.mutagenic_zone = mutagenic_zone
        self.pattern = pattern
        self.seed_extend_settings = seed_extend_settings
        self.config = config

    def _build_kmer_index(self) -> Dict[str, List[int]]:
        mutagenic_zone = self.mutagenic_zone
        n_mutagenic_zone = len(mutagenic_zone)
        k = self.seed_extend_settings.k
        idx = defaultdict(list)

        for i in range(n_mutagenic_zone - k + 1):
            kmer = mutagenic_zone[i:i+k]
            idx[kmer].append(i)
        return idx

    def _find_seed_hits(self, kmer_indexes: Dict[str, List[int]]) -> List[Tuple[int, int]]:
        hits: List[Tuple[int, int]] = []
        mutagenic_zone = self.mutagenic_zone
        k = self.seed_extend_settings.k
        
        for qpos in range(len(mutagenic_zone) - k + 1):
            kmer = mutagenic_zone[qpos:qpos + k]
            for rpos in kmer_indexes.get(kmer, ()):
                hits.append((qpos, rpos))
        return hits

    def _cluster_hits_by_offset(self, hits: List[Tuple[int, int]]) -> List[SeedExtendCluster]:
        offset_tolerance = self.seed_extend_settings.offset_tolerance
        min_seeds = self.seed_extend_settings.min_seed

        offsets = [r - q for (q, r) in hits]
        # group offsets by rounding/nearby bucket
        # Greedy clustering by offset value
        offsets_sorted = sorted(set(offsets))
        clusters = []
        for off in offsets_sorted:
            placed = False
            for cl in clusters:
                #cluster_mean = statistics.median(cl["offsets"]) #<<<<<<<<note: figure out if the mean is actually the metric that should be used to measure distance from
                if abs(cl["offsets"][0] - off) <= offset_tolerance: 
                    cl['offsets'].append(off)
                    placed = True
                    break
            if not placed:
                clusters.append({'offsets': [off]})
        # For each cluster, gather supporting hits
        clusters_out: List[SeedExtendCluster] = []
        for cl in clusters:
            off_vals = set(cl['offsets'])
            supporting = [(q, r) for (q, r) in hits if (r - q) in off_vals]
            if len(supporting) < min_seeds: #min_seeds determines the minimum number of kmer's that have to be close enough to eachother for us to consider the cluster. To my understanding, if the kmer size is 4 and min_seed is 2, then a single substring match of 5 bp would qualify as a valid cluster
                continue
            q_positions = [q for q, _ in supporting]
            r_positions = [r for _, r in supporting]
            clusters_out.append(SeedExtendCluster(
                offset=int(sum(off_vals) / len(off_vals)),
                hits=supporting,
                qmin=min(q_positions),
                qmax=max(q_positions),
                rmin=min(r_positions),
                rmax=max(r_positions),
                n_seeds=len(supporting)
            ))
        return clusters_out

    def _parse_cigar_string(self, cigar_str: str) -> List[Tuple[str, int]]:
        ops = []
        cigar_op_re = re.compile(r'(\d+)([MIDNSHP=X])')
        for m in cigar_op_re.finditer(cigar_str):
            l = int(m.group(1))
            op = m.group(2)
            ops.append((op, l))
        return ops

    def _analyze_from_cigar(self, cigar_str: str, sub_mutagenic_zone: str, mutagenic_zone_start: int) -> ImperfectAlignmentEvent:
        #Fix: get rid of below 2 variables and rename stuff in this function
        pattern = self.pattern
        sub_mutagenic_zone = sub_mutagenic_zone

        ops: List[Tuple[str, int]] = self._parse_cigar_string(cigar_str)
        if len(ops) == 0:
            raise ValueError("Empty Cigar")
        # count consumption
        sub_mutagenic_zone_consumed = sum(l for (op, l) in ops if op in ('M', 'D', '=', 'X'))
        pattern_consumed = sum(l for (op, l) in ops if op in ('M', 'I', '=', 'X'))
        # soft clip lengths (if any)
        left_soft = ops[0][1] if ops[0][0] == 'S' else 0
        right_soft = ops[-1][1] if ops[-1][0] == 'S' else 0
        pattern_total = pattern_consumed + left_soft + right_soft

        beg_pattern = left_soft
        end_pattern = beg_pattern + pattern_consumed + right_soft - 1

        # If user expects entire pattern aligned (semi-global), validate this:
        if pattern_total != len(pattern):
            # Not fatal, but warn via exception text so caller notices.
            # (We could also choose to silently continue; you asked for clarity.)
            raise ValueError(
                f"CIGAR-consumed pattern length ({pattern_total}) != len(pattern) ({len(sub_mutagenic_zone)}). "
                "For a true semi-global alignment the entire pattern should be consumed. "
                "Check that you called an sg_trace_* function (no soft-clips) or that your "
                "CIGAR matches expectations."
            )

        # Compute beginning sub_mutagenic_zone coordinate if missing, try inclusive/exclusive interpretations
        if ops[0][0] == 'D':
            beg_sub_mutagenic_zone = ops[0][1]
        else:
            beg_sub_mutagenic_zone = 0

        # Now do the actual mapping walk; check bounds as we go
        sub_mutagenic_zone_pos = int(beg_sub_mutagenic_zone)
        pattern_pos = int(beg_pattern)

        matches = 0
        mismatches: List[Tuple[int, str, str]]= []
        insertions: List[Tuple[int, str]] = []  # (sub_mutagenic_zone_pos_where_insertion_happens, pattern_pos_start, length, seq)
        deletions: List[Tuple[int, str]] = []   # (sub_mutagenic_zone_pos_start, length, seq)
        match_positions_sub_mutagenic_zone = []  # all ref indices that are exact matches

        if ops[0][0] == 'D':
            ops = ops[1:]

        op_save = None
        L_save = None

        for op, L in ops:
            if pattern_pos >= len(pattern) or sub_mutagenic_zone_pos >= len(sub_mutagenic_zone) or pattern[pattern_pos] == 'N':
                if pattern_pos < len(pattern) and pattern[pattern_pos] == 'N':
                    pattern_pos = len(pattern)
                op_save = op
                L_save = L
                break
            if op in ('M', '=', 'X'):
                for i in range(L):
                    #r_idx = sub_mutagenic_zone_pos + i - beg_ref
                    #q_idx = pattern_pos + i - beg_pattern
                    r_idx = sub_mutagenic_zone_pos + i 
                    q_idx = pattern_pos + i
                    #r_idx = i + L
                    #q_idx = i + L
                    if r_idx < 0 or r_idx >= len(sub_mutagenic_zone) or q_idx < 0 or q_idx >= len(pattern):
                        raise IndexError(f"Index out of range during CIGAR walk: r_idx={r_idx}, q_idx={q_idx}")
                    rbase = sub_mutagenic_zone[r_idx]
                    qbase = pattern[q_idx]
                    if rbase == qbase:
                        matches += 1
                        match_positions_sub_mutagenic_zone.append(r_idx)
                    else:
                        #mismatches.append((r_idx, q_idx, rbase, qbase))
                        mismatches.append((q_idx+1, qbase, rbase)) # the relative r_idx should be r_idx - true_beg_sub_mutagenic_zone
                        #mismatches.append((r_idx-true_beg_sub_mutagenic_zone, q_idx-true_beg_ref, rbase, qbase))
                sub_mutagenic_zone_pos += L
                pattern_pos += L
            elif op == 'I': #treat I's like deletions not insertions
                # convention: insertion at current sub_mutagenic_zone_pos (between ref_pos-1 and ref_pos)
                for i in range(L):
                    q_idx = pattern_pos + i
                    q_base = pattern[q_idx]
                    deletions.append((q_idx+1, q_base))
                pattern_pos += L
            elif op == 'D': # treat D's like insertions not deletions
                for i in range(L):
                    r_idx = sub_mutagenic_zone_pos + i
                    r_base = sub_mutagenic_zone[r_idx]
                    insertions.append((pattern_pos, r_base))
                sub_mutagenic_zone_pos += L
                #pattern_pos += 1
            elif op == 'S':
                # soft clip: consumes pattern only (these were accounted for earlier)
                pattern_pos += L
            elif op == 'H':
                # hard clip: consumes neither
                pass
            elif op == 'N':
                # skipped region in sub_mutagenic_zoneerence (like intron) - consumes ref
                sub_mutagenic_zone_pos += L
            else:
                raise ValueError(f"Unhandled CIGAR op '{op}'")

        #incase sub_mutagenic_zoneerence is smaller then pattern and there are a bunch of deletions at the end:
        if pattern_pos < len(pattern) and op_save == 'I' and L_save:
            for i in range(L_save):
                q_idx = pattern_pos + i
                q_base = pattern[q_idx]
                deletions.append((q_idx+1, q_base))

        n_mistakes = len(insertions) + len(deletions) + len(mismatches)
        score = (len(pattern) - n_mistakes) / len(pattern)

        end_sub_mutagenic_zone = sub_mutagenic_zone_pos - 1

        output = ImperfectAlignmentEvent(
            full_telomer_start_index=None,
            mutagenic_zone_start_index=mutagenic_zone_start + beg_sub_mutagenic_zone,
            mutagenic_zone_end_index=mutagenic_zone_start + end_sub_mutagenic_zone,
            insertion_events=insertions,
            deletion_events=deletions,
            mismatch_events=mismatches,
            score = score 
        )

        return output

    def _simple_alignment(self) -> List[ImperfectAlignmentEvent]:
        alignments: List[ImperfectAlignmentEvent] = []
        max_mistakes: int = self.config.maximum_alignment_mutations
        n_pattern = len(self.pattern)
        n_mutagenic_zone = len(self.mutagenic_zone)

        window_size: int = min(n_pattern + n_pattern // max_mistakes, n_mutagenic_zone)

        region: str = self.mutagenic_zone
        n_region: int = len(region)

        match_score: int = self.seed_extend_settings.alignment_settings.match
        mismatch_score: int = self.seed_extend_settings.alignment_settings.mismatch
        gap_open_score: int = self.seed_extend_settings.alignment_settings.gap_open
        gap_extend_score: int = self.seed_extend_settings.alignment_settings.gap_extend
        matrix = parasail.matrix_create("ACGT", match_score, mismatch_score)

        slide_end: int = n_region - min(n_pattern, n_region) + 1
        for start in range(0, slide_end, 1):
            mutagenic_zone_start: int = start 
            sub_mutagenic_zone: str = region[start: start + window_size]

            parasail_alignment: Result = parasail.sg_trace_scan(self.pattern, sub_mutagenic_zone, gap_open_score, gap_extend_score, matrix)

            try:
                #Fix: make sure the below gives the correct cigar strings
                cigar_str: str = str(parasail_alignment.cigar.decode)
            except:
                raise ValueError("your version of parasail does not give access to the cigar string with the sg_trace_scan function. Please use different version of parasail")

            cigar_analysis: ImperfectAlignmentEvent = self._analyze_from_cigar(cigar_str, sub_mutagenic_zone, mutagenic_zone_start)

            n_mistakes = len(cigar_analysis.insertion_events) + len(cigar_analysis.deletion_events) + len(cigar_analysis.mismatch_events)

            if n_mistakes <= max_mistakes:
                alignments.append(cigar_analysis)

        return alignments

        

    def _extend_cluster(self, cluster: SeedExtendCluster) -> List[ImperfectAlignmentEvent]:
        alignments: List[ImperfectAlignmentEvent] = []
        max_mistakes = self.config.maximum_alignment_mutations
        # max_events = len(query) // 25
        n_pattern = len(self.pattern)
        n_mutagenic_zone = len(self.mutagenic_zone)

        window_size = min(n_pattern + n_pattern // max_mistakes, n_mutagenic_zone)

        cluster_start = max(0, cluster.rmin - cluster.qmin - self.seed_extend_settings.flank)
        cluster_end = min(n_mutagenic_zone, cluster_start + n_pattern + self.seed_extend_settings.flank)
        region = self.mutagenic_zone[cluster_start:cluster_end]
        n_region = len(region)

        match_score = self.seed_extend_settings.alignment_settings.match
        mismatch_score = self.seed_extend_settings.alignment_settings.mismatch
        gap_open_score = self.seed_extend_settings.alignment_settings.gap_open
        gap_extend_score = self.seed_extend_settings.alignment_settings.gap_extend
        matrix = parasail.matrix_create("ACGT", match_score, mismatch_score)

        slide_end = n_region - min(n_pattern, n_region) + 1
        for start in range(0, slide_end, 1):
            mutagenic_zone_start = cluster_start + start
            sub_mutagenic_zone = region[start: start + window_size]

            parasail_alignment: Result = parasail.sg_trace_scan(self.pattern, sub_mutagenic_zone, gap_open_score, gap_extend_score, matrix)

            try:
                #Fix: make sure the below gives the correct cigar strings
                cigar_str = str(parasail_alignment.cigar.decode)
            except:
                raise ValueError("your version of parasail does not give access to the cigar string with the sg_trace_scan function. Please use different version of parasail")

            cigar_analysis: ImperfectAlignmentEvent = self._analyze_from_cigar(cigar_str, sub_mutagenic_zone, mutagenic_zone_start)

            n_mistakes = len(cigar_analysis.insertion_events) + len(cigar_analysis.deletion_events) + len(cigar_analysis.mismatch_events)

            if n_mistakes <= max_mistakes:
                alignments.append(cigar_analysis)


        return alignments

    def _is_subset(self, aln_a: ImperfectAlignmentEvent, aln_b: ImperfectAlignmentEvent) -> bool:
        start_diff = (aln_a.mutagenic_zone_start_index) - (aln_b.mutagenic_zone_start_index)
        if abs(start_diff) > 10:
            return False

        deletions_a: Set[int] = {x[0] for x in aln_a.deletion_events}
        insertions_a: set[int] = {x[0] for x in aln_a.insertion_events}
        mismatches_a: Set[int] = {x[0] for x in aln_a.mismatch_events}

        for i in aln_b.insertion_events:
            if i[0] not in insertions_a:
                return False

        for d in aln_b.deletion_events:
            if d[0] not in deletions_a:
                return False

        for m in aln_b.mismatch_events:
            if m[0] not in mismatches_a:
                return False
    
        return True


    def _filter_redundant_alignments_by_error(self, alignments: List[ImperfectAlignmentEvent]) -> List[ImperfectAlignmentEvent]:

        kept: List[ImperfectAlignmentEvent] = []
        added_start_pos: Set = set()

        alignments = sorted(alignments, key=lambda x: x.score, reverse=True)

        for i, aln_a in enumerate(alignments):
            if aln_a.mutagenic_zone_start_index in added_start_pos:
                continue
            redundant = False
            #errors_a = aln_a["insertions"] + aln_a["deletions"] + aln_a["mismatches"]

            for aln_b in kept:
                #errors_b = aln_b["insertions"] + aln_b["deletions"] + aln_b["mismatches"]

                # If A’s errors are covered by B’s errors, then A is redundant (worse).
                if self._is_subset(aln_a, aln_b):
                    redundant = True
                    break
        
            if not redundant:
                kept.append(aln_a)
                added_start_pos.add(aln_a.mutagenic_zone_start_index)

        return kept

    def _best_non_overlapping_alignments(self, alignments: List[ImperfectAlignmentEvent]) -> List[ImperfectAlignmentEvent]:
        if not alignments:
            return []

        # sort by end (classic interval scheduling order)
        aligns = sorted(alignments, key=lambda x: x.mutagenic_zone_start_index)
        n = len(aligns)

        # dp[i] stores the best solution that ENDS at alignment i
        # tuple layout:
        #   (count, gaps_including_start, score, first_start, path_list)
        dp: List[Optional[Tuple[int, int, int, int, List[ImperfectAlignmentEvent]]]] = [None] * n

        #for i, (si, ei, sco_i, *rest_i) in enumerate(aligns):
        for i in range(n):
            curr_align = aligns[i]
            si = curr_align.mutagenic_zone_start_index
            ei = curr_align.mutagenic_zone_end_index
            sco_i = curr_align.score
            start_gap = 0 if si == 0 else 1 # adds an extra gap if the start index doesn't equal 0
            best = (1, start_gap, sco_i, si, [aligns[i]])

            # try to extend from any non-overlapping j < i
            for j in range(i):
                #sj, ej, sco_j, *rest_j = aligns[j]
                curr_cmp_align = aligns[j]
                sj = curr_cmp_align.mutagenic_zone_start_index
                ej = curr_cmp_align.mutagenic_zone_end_index
                sco_j = curr_cmp_align.score
                if ej < si:  # non-overlapping (abutting allowed)
                    cnt, gaps_w_start, sc, first_start, path = dp[j]

                    # add 0 if abut, else +1 gap if there's separation
                    gap_add = 0 if ej+1 == si else 1

                    cand = (
                        cnt + 1,
                        gaps_w_start + gap_add,
                        sc + sco_i,
                        first_start,
                        path + [aligns[i]],
                    )

                    # lexicographic preference: count ↑, gaps ↓, score ↑
                    if (
                        cand[0] > best[0] or
                        (cand[0] == best[0] and cand[1] < best[1]) or
                        (cand[0] == best[0] and cand[1] == best[1] and cand[2] > best[2])
                    ):
                        best = cand

            dp[i] = best

        # Final selection across endpoints i:
        # add the END-gap (+1 if the chain doesn't reach reference end).
        best_overall = None
        best_key = None

        for i, state in enumerate(dp):
            cnt, gaps_w_start, sc, first_start, path = state
            last_end = path[-1].mutagenic_zone_end_index  # end of the last chosen alignment
            end_gap = 0 if last_end+1 == len(self.mutagenic_zone) or self.seed_extend_settings.last_zone else 1
            total_gaps = gaps_w_start + end_gap

            key: Tuple[int, int, int] = (cnt, -total_gaps, sc)  # count ↑, gaps ↓, score ↑
            if best_overall is None or key > best_key:
                best_overall = {
                    "count": cnt,
                    "gaps": total_gaps,
                    "score": sc,
                    "path": path
                }
                best_key = key

        return best_overall['path']
        

    def execute(self) -> List[ImperfectAlignmentEvent]:
        skip_seeding: Optional[bool] = self.config.skip_seeding
        if skip_seeding and not skip_seeding:
            kmer_indexes: Dict[str, List[int]] = self._build_kmer_index()
            hits: List[Tuple[int, int]] = self._find_seed_hits(kmer_indexes)
            if not hits:
                return []

            clusters: List[SeedExtendCluster] = self._cluster_hits_by_offset(hits)

            if not clusters:
                return []

            results: List[ImperfectAlignmentEvent] = []
            all_used_regions = set() #Fix: add the use of the all_used_regions
            for cl in clusters:
                #Fix, add the use of the all_used_regions
                #alns, used_regions = self._extend_cluster(cl)
                alns: List[ImperfectAlignmentEvent] = self._extend_cluster(cl)
                if alns:
                    results += alns

        else:
            results: List[ImperfectAlignmentEvent] = self._simple_alignment()

        results = self._filter_redundant_alignments_by_error(results)
        final_result: List[ImperfectAlignmentEvent] = self._best_non_overlapping_alignments(results)

        return final_result

