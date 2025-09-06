"""
seed_and_extend.py

Seed-hit + extend pipeline for semiglobal alignments using parasail.

Requirements:
    pip install parasail

Author: ChatGPT (adapt for your environment)
"""

from collections import defaultdict, Counter
import parasail
from memory_profiler import profile
import statistics

from bisect import bisect_right

import parasail_functions
import parse_csv

# ---------- Masking code ----------
#import bisect

#def add_mask_interval(mask_intervals, start, end):
    #"""
    #Add interval [start,end] into mask_intervals (merge overlaps). mask_intervals
    #is kept sorted and non-overlapping.
    #"""
    #if start > end:
        #return
    ## find insertion point
    #i = bisect.bisect_left(mask_intervals, (start, end))
    ## merge with left neighbor if overlapping
    #if i > 0 and mask_intervals[i-1][1] >= start - 1:
        #i -= 1
        #start = min(start, mask_intervals[i][0])
        #end = max(end, mask_intervals[i][1])
        #mask_intervals.pop(i)
    ## merge with following intervals
    #while i < len(mask_intervals) and mask_intervals[i][0] <= end + 1:
        #start = min(start, mask_intervals[i][0])
        #end = max(end, mask_intervals[i][1])
        #mask_intervals.pop(i)
    #mask_intervals.insert(i, (start, end))


#def iterative_mask_and_search(query,
                              #reference,
                              #k=11,
                              #flank=60,
                              #match=2,
                              #mismatch=-2,
                              #gap_open=5,
                              #gap_extend=1,
                              #min_identity=0.90,
                              #offset_tolerance=8,
                              #min_seeds=2,
                              #score_threshold=None,
                              #max_iterations=10):
    #"""
    #Dynamic iterative mask-and-search:
      #- Build k-mer index once
      #- Produce seed hits
      #- While seeds remain:
          #- filter hits by mask
          #- cluster
          #- extend clusters into candidates
          #- select best candidate
          #- accept candidate -> add to mask_intervals
          #- remove masked hits and continue
    #Returns list of accepted alignments
    #"""

    ## 1) Prep
    #query = query.upper()
    #reference = reference.upper()
    #qlen = len(query)

    #index = build_kmer_index(reference, k)
    #seed_hits = find_seed_hits(query, index, k)  # list of (qpos,rpos)

    ## Sort seed_hits for deterministic behavior (optional)
    #seed_hits.sort()

    #mask_intervals = []  # list of (start,end), merged, sorted
    #accepted = []

    #iters = 0
    #while seed_hits and iters < max_iterations:
        #iters += 1

        ## 2) Apply soft mask: remove hits whose rpos in any mask interval
        #seed_hits = apply_soft_mask_to_hits(seed_hits, mask_intervals)
        #if not seed_hits:
            #break

        ## 3) Cluster the current seed hits
        #clusters = cluster_hits_by_offset(seed_hits,
                                          #offset_tolerance=offset_tolerance,
                                          #min_seeds=min_seeds)
        #if not clusters:
            #break

        ## 4) Extend all clusters to get candidates (you can also extend only top-N clusters)
        #candidates = []
        #for cl in clusters:
            #aln = extend_cluster(query, reference, cl,
                                 #flank=flank,
                                 #match=match,
                                 #mismatch=mismatch,
                                 #gap_open=gap_open,
                                 #gap_extend=gap_extend,
                                 #min_identity=min_identity)
            #if aln is not None:
                ## attach cluster info and unify ref coords if present
                #aln['cluster'] = cl
                ## ensure ref_start/end exist and are ints
                #s = aln.get('ref_start')
                #e = aln.get('ref_end')
                #if s is None or e is None:
                    #continue
                #aln['ref_start'] = int(s)
                #aln['ref_end'] = int(e)
                #candidates.append(aln)

        #if not candidates:
            ## nothing passed extension filters
            #break

        ## 5) Pick best candidate (by score, break ties by identity)
        #candidates.sort(key=lambda c: (c['score'], c['identity']), reverse=True)
        #best = candidates[0]

        ## optional score threshold
        #if score_threshold is not None and best['score'] < score_threshold:
            #break

        ## 6) Accept it
        #accepted.append(best)

        ## 7) Add its ref interval to mask_intervals (merge)
        #add_mask_interval(mask_intervals, best['ref_start'], best['ref_end'])

        ## 8) Remove seed hits that fall in the new mask (will be done at loop start)
        #seed_hits = apply_soft_mask_to_hits(seed_hits, mask_intervals)

    #return accepted



# ---------- Utilities to select non-overlapping alignments ----------

#from bisect import bisect_right

def overlaps(a_start, a_end, b_start, b_end, allow_partial_overlap=0):
    """
    True if intervals overlap by more than allow_partial_overlap bases.
    allow_partial_overlap can be 0 (no overlap allowed) or >0 to permit small overlaps.
    """
    # overlap length = min(a_end,b_end) - max(a_start,b_start) + 1
    ov = min(a_end, b_end) - max(a_start, b_start) + 1
    return ov > allow_partial_overlap

def greedy_select(candidates, key='score', allow_partial_overlap=0):
    """
    Greedy selection: sort by key desc and pick non-overlapping intervals.
    allow_partial_overlap=int: number of bases allowed to overlap (slack).
    Returns list of selected candidates (in descending key order).
    """
    # Sort by key desc; break ties by longer alignment then higher identity if present
    def sort_key(c):
        return (c.get(key, 0), c.get('ref_end') - c.get('ref_start') + 1, c.get('identity', 0))
    candidates_sorted = sorted(candidates, key=sort_key, reverse=True)

    selected = []
    used_intervals = []
    for c in candidates_sorted:
        s, e = c['ref_start'], c['ref_end']
        # check overlap with every already selected
        bad = False
        for us, ue in used_intervals:
            if overlaps(s, e, us, ue, allow_partial_overlap=allow_partial_overlap):
                bad = True
                break
        if not bad:
            selected.append(c)
            used_intervals.append((s, e))
    return selected


# Weighted Interval Scheduling (optimal by score)
def weighted_interval_scheduling(candidates, weight_key='score', allow_partial_overlap=0):
    """
    Optimal selection by weight using DP:
    - candidates: list of dicts with 'ref_start', 'ref_end', and weight_key
    - allow_partial_overlap: allowed overlap in bases (0 for strict non-overlap)
    Returns list of selected candidates (optimal by sum of weights).
    Complexity: O(n log n)
    """

    # 1) Sort candidates by end position
    items = sorted(candidates, key=lambda c: c['ref_end'])
    n = len(items)
    if n == 0:
        return []

    # 2) Precompute `p[j]`: index of the rightmost interval i < j that doesn't overlap j
    # We will use strict non-overlap subject to allow_partial_overlap slack.
    # Create an array of end positions for binary search:
    ends = [it['ref_end'] for it in items]

    # Helper to find the last index i < j that doesn't overlap j
    def find_last_non_overlapping(j):
        # find rightmost i with items[i].ref_end < items[j].ref_start - allow_partial_overlap
        target = items[j]['ref_start'] - allow_partial_overlap - 1
        # bisect_right returns insertion point; subtract 1 to get index <= target
        i = bisect_right(ends, target) - 1
        return i  # may be -1 if none

    p = [find_last_non_overlapping(j) for j in range(n)]

    # 3) DP table: M[j] = max weight considering intervals[0..j]
    M = [0] * n
    # For reconstruction:
    take = [False] * n

    def weight(j):
        return items[j].get(weight_key, 0)

    for j in range(n):
        incl = weight(j) + (M[p[j]] if p[j] != -1 else 0)
        excl = M[j-1] if j > 0 else 0
        if incl > excl:
            M[j] = incl
            take[j] = True
        else:
            M[j] = excl
            take[j] = False

    # 4) Reconstruct solution
    selected = []
    j = n - 1
    while j >= 0:
        if take[j]:
            selected.append(items[j])
            j = p[j]
        else:
            j -= 1

    selected.reverse()  # increasing ref_end order
    return selected


# ---------- Mask utilities (optional) ----------

def mask_reference_intervals(reference, intervals):
    """
    Hard mask: replace bases in each interval (inclusive) with 'N'.
    intervals: list of (start,end) pairs (0-based inclusive).
    Returns a new masked reference string.
    """
    if not intervals:
        return reference
    ref_list = list(reference)
    for (s, e) in intervals:
        if s < 0: s = 0
        if e >= len(ref_list): e = len(ref_list) - 1
        for i in range(s, e + 1):
            ref_list[i] = 'N'
    return ''.join(ref_list)

def apply_soft_mask_to_hits(seed_hits, mask_intervals):
    """
    Given list of seed hits (qpos,rpos) keep only those where rpos not in any mask interval.
    mask_intervals: list of (s,e)
    Returns filtered seed_hits list.
    """
    if not mask_intervals:
        return seed_hits
    # sort mask_intervals for efficient testing
    mask_intervals_sorted = sorted(mask_intervals, key=lambda x: x[0])
    out = []
    import bisect
    starts = [mi[0] for mi in mask_intervals_sorted]
    for qpos, rpos in seed_hits:
        # find the mask whose start is <= rpos (bisect_right-1)
        idx = bisect.bisect_right(starts, rpos) - 1
        if idx >= 0:
            s, e = mask_intervals_sorted[idx]
            if rpos >= s and rpos <= e:
                continue
        out.append((qpos, rpos))
    return out


def build_kmer_index(reference: str, k: int):
    """Return dict: kmer -> list of positions in reference."""
    idx = defaultdict(list)
    L = len(reference)
    for i in range(L - k + 1):
        kmer = reference[i:i + k]
        idx[kmer].append(i)
    return idx


def find_seed_hits(query: str, index: dict, k: int):
    """
    Return list of seed hits as (qpos, rpos).
    Only exact k-mer seeds are used here (fast and simple).
    """
    hits = []
    for qpos in range(len(query) - k + 1):
        kmer = query[qpos:qpos + k]
        for rpos in index.get(kmer, ()):
            hits.append((qpos, rpos))
    return hits


def cluster_hits_by_offset(hits, offset_tolerance=8, min_seeds=2):
    """
    Convert (qpos,rpos) hits to clusters by offset = rpos - qpos.
    Returns list of clusters, each cluster is dict:
      {
        'offset': representative_offset,
        'hits': list of (qpos,rpos),
        'qmin','qmax','rmin','rmax'
      }
    offset_tolerance controls how close offsets must be to cluster together.
    """
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
    clusters_out = []
    for cl in clusters:
        off_vals = set(cl['offsets'])
        supporting = [(q, r) for (q, r) in hits if (r - q) in off_vals]
        if len(supporting) < min_seeds: #min_seeds determines the minimum number of kmer's that have to be close enough to eachother for us to consider the cluster. To my understanding, if the kmer size is 4 and min_seed is 2, then a single substring match of 5 bp would qualify as a valid cluster
            continue
        q_positions = [q for q, _ in supporting]
        r_positions = [r for _, r in supporting]
        clusters_out.append({
            'offset': int(sum(off_vals) / len(off_vals)),
            'hits': supporting,
            'qmin': min(q_positions),
            'qmax': max(q_positions),
            'rmin': min(r_positions),
            'rmax': max(r_positions),
            'n_seeds': len(supporting),
        })
    return clusters_out

def filter_redundant_alignments(alignments, overlap_thresh=0.8, score_diff_thresh=0.05):
    """
    Remove highly redundant alignments, keeping only meaningful alternates.
    
    Args:
        alignments (list of dict): each with 'window_range', 'score', 'identity'.
        overlap_thresh (float): fraction overlap to be considered redundant.
        score_diff_thresh (float): relative score difference below which alignments
                                   are considered redundant.
    """
    kept = []
    
    for aln in sorted(alignments, key=lambda x: x["score"], reverse=True):
        redundant = False
        aln_start, aln_end = aln["window_range"]
        
        for kept_aln in kept:
            k_start, k_end = kept_aln["window_range"]

            # Compute overlap
            overlap = max(0, min(aln_end, k_end) - max(aln_start, k_start))
            overlap_frac = overlap / min(aln_end - aln_start, k_end - k_start)

            # Compare scores/identity
            score_diff = abs(aln["score"] - kept_aln["score"]) / max(aln["score"], kept_aln["score"])

            if overlap_frac >= overlap_thresh and score_diff <= score_diff_thresh:
                redundant = True
                break
        
        if not redundant:
            kept.append(aln)

    return kept


def extend_cluster(query,
                   reference,
                   cluster,
                   flank=60,
                   window_size = 150,
                   step_size = 10,
                   match=2,
                   mismatch=-3,
                   gap_open=5,
                   gap_extend=2,
                   min_identity=0.88,
                   used_regions=None):
    """
    Extract a ref window around the cluster and run a semi-global alignment (extension).
    Returns a dict with alignment info if it passes filters; else None.

    Notes:
     - flank: number of reference bases added on both sides to allow indels and contextual alignment.
     - min_identity: fraction of query bases that must be exact matches (matches / query_len).
    """
    alignments = []
    if used_regions is None:
        used_regions = set() #bit inefficient and pointless
    
    max_mistakes = len(query) // 12
    max_events = len(query) // 25
    #max_mistakes = len(query) // 10
    #max_mistakes = len(query) // 8

    qlen = len(query)
    n_ref = len(reference)
    n_query = len(query)

    window_size = min(n_query + n_query//12, n_ref)

    sliding_size = min(n_query + n_query//12//2, n_ref)

    
    cluster_start = max(0, cluster['rmin'] - cluster['qmin'] - flank)
    cluster_end = min(len(reference), cluster_start + qlen + flank)
    region = reference[cluster_start:cluster_end]
    # scoring matrix for DNA (simple)
    matrix = parasail.matrix_create("ACGT", match, mismatch)

    #start alignmnents
    #for r_offset in range(0, max_mistakes, step_size):
        ##sub_ref_start = cluster_start + start # this is the actual value in the full reference string
        #sub_ref_start = cluster_start
        ##if sub_ref_start in used_regions:
            ##continue
        ##else:
            ##used_regions.add(sub_ref_start)
        #sub_ref_end = sub_ref_start + window_size

        #if n_query < n_ref:
            #break


        ## Perform semiglobal alignment
        #sub_ref = region[cluster_start:cluster_start + r_offset + n_query]

        #if len(sub_ref) > n_query:
            #pad_length = len(sub_ref) - n_query
            #query_padded = query + "N" * pad_length
        #else:
            #query_padded = query

        #query_padded = query
        
        
        ##res = parasail.nw_trace_scan(query_padded, sub_ref, gap_open, gap_extend, matrix)
        #res = parasail.sg_trace_scan(query_padded, sub_ref, gap_open, gap_extend, matrix)
        ##result = parasail.sg_dx_trace(query_seq, sub_ref, gap_open, gap_extend, matrix)

        #cigar_str = parasail_functions.try_get_cigar_string(res)
        #if cigar_str is None:
            #raise RuntimeError("Could not obtain CIGAR string from parasail result. Ensure you used a *_trace_* semi-global function.")


        #analysis = parasail_functions.analyze_alignment_from_cigar(
            #query=query_padded,
            #ref=sub_ref,
            #cigar_str=cigar_str,
        #)

        
        #n_mistakes = len(analysis['insertions']) + len(analysis['deletions']) + len(analysis['mismatches'])
        #score = (n_query - n_mistakes) / n_query

        ##identity = analysis['matches'] / len(query)

        #if n_mistakes <= max_mistakes:
        ##if True:

            #aln = {
                ##'score': res.score, #check if score exists
                ##'score': (n_query-n_mistakes)/n_query, #check if score exists
                #'score': res.score,
                #'matches': analysis['matches'],
                #'aligned_length': analysis['end_ref'] - analysis['beg_ref'] + 1,
                ##'identity': round(identity, 4),
                #'ref_start': analysis['beg_ref'],
                #'ref_end': analysis['end_ref'],
                #'cigar': cigar_str,
                #'n_seeds': cluster.get('n_seeds', 0),
                #'mismatches': analysis['mismatches'],
                #'insertions': analysis['insertions'],
                #'deletions': analysis['deletions'],
                #'true_aligned_window_start': sub_ref_start,
                #'absolute_ref_start': sub_ref_start + analysis['beg_ref'],
                #'absolute_ref_end': sub_ref_start + analysis['end_ref']
            #}

            #alignments.append(aln)


    #slide window
    #for start in range(0, len(region) - sliding_size+1, step_size):
    for start in range(0, len(region) - min(n_query, len(region)) + 1, step_size):
        sub_ref_start = cluster_start + start # this is the actual value in the full reference string
        #if sub_ref_start in used_regions:
            #continue
        #else:
            #used_regions.add(sub_ref_start)
        sub_ref_end = sub_ref_start + window_size


        # Perform semiglobal alignment
        sub_ref = region[start:start + window_size]

        if len(sub_ref) > n_query:
            pad_length = len(sub_ref) - n_query
            query_padded = query + "N" * pad_length
        else:
            query_padded = query

        query_padded = query
        
        
        #res = parasail.nw_trace_scan(query_padded, sub_ref, gap_open, gap_extend, matrix)
        res = parasail.sg_trace_scan(query_padded, sub_ref, gap_open, gap_extend, matrix)
        #result = parasail.sg_dx_trace(query_seq, sub_ref, gap_open, gap_extend, matrix)

        cigar_str = parasail_functions.try_get_cigar_string(res)
        if cigar_str is None:
            raise RuntimeError("Could not obtain CIGAR string from parasail result. Ensure you used a *_trace_* semi-global function.")


        analysis = parasail_functions.analyze_alignment_from_cigar(
            query=query_padded,
            ref=sub_ref,
            cigar_str=cigar_str,
        )

        insertions_events, _, _= parse_csv.count_events(analysis['insertions'], "insertions")
        deletion_events, _, _ = parse_csv.count_events(analysis['deletions'], "deletions")
        mismatch_events = len(analysis['mismatches'])
        n_events = insertions_events + deletion_events + mismatch_events

        
        n_mistakes = len(analysis['insertions']) + len(analysis['deletions']) + len(analysis['mismatches'])
        score = (n_query - n_mistakes) / n_query

        #identity = analysis['matches'] / len(query)

        if n_mistakes <= max_mistakes or n_events < max_events:
        #if True:

            aln = {
                #'score': res.score, #check if score exists
                #'score': (n_query-n_mistakes)/n_query, #check if score exists
                'score': res.score,
                'matches': analysis['matches'],
                'aligned_length': analysis['end_ref'] - analysis['beg_ref'] + 1,
                #'identity': round(identity, 4),
                'ref_start': analysis['beg_ref'],
                'ref_end': analysis['end_ref'],
                'cigar': cigar_str,
                'n_seeds': cluster.get('n_seeds', 0),
                'mismatches': analysis['mismatches'],
                'insertions': analysis['insertions'],
                'deletions': analysis['deletions'],
                'true_aligned_window_start': sub_ref_start,
                'absolute_ref_start': sub_ref_start + analysis['beg_ref'],
                'absolute_ref_end': sub_ref_start + analysis['end_ref']
            }

            alignments.append(aln)


    return alignments, used_regions

    # If Parasail produced a CIGAR (trace variant), try to decode and compute query-consumed length.
    # Some parasail versions give .cigar.decode attribute; we'll try to be flexible.
    try:
        if hasattr(aln, 'cigar') and aln.cigar is not None:
            # decode to SAM CIGAR string if available
            # Many parasail wrappers have aln.cigar.decode() or aln.cigar.decode
            if callable(getattr(aln.cigar, 'decode', None)):
                cigar_str = aln.cigar.decode()
            else:
                cigar_str = getattr(aln.cigar, 'decode', None) or str(aln.cigar)
            # If CIGAR object contains begin positions:
            # Some parasail CIGARs include .ref_begin and .query_begin
            try:
                a = aln.cigar
                # use attributes if provided
                aln_query_start = getattr(a, 'query_begin', None)
                aln_ref_start = getattr(a, 'ref_begin', None)
                aln_query_end = getattr(a, 'query_end', None)
                aln_ref_end = getattr(a, 'ref_end', None)
            except Exception:
                pass
    except Exception:
        # ignore cigar extraction errors — not all builds have trace enabled
        cigar_str = None

    # Basic sanity checks to ensure entire query was aligned.
    # Prefer to use explicit query-begin/end if available (safe). Otherwise, fallback.
    # If query begin/end are not present, we can attempt a weaker check:
    query_consumed = None
    try:
        if aln_query_start is not None and aln_query_end is not None:
            query_consumed = (aln_query_end - aln_query_start + 1)
        else:
            # some sg_stats results expose `.length` (alignment length) and `.matches`
            # but `.length` is alignment column count (includes gaps), so can't directly
            # infer query consumption. We'll assume full query if `.length >= qlen` and
            # `.matches` is present. This is a conservative heuristic; if you need exact
            # behavior parse CIGAR by using sg_trace variant instead.
            if hasattr(aln, 'length'):
                if aln.length >= qlen:
                    query_consumed = qlen
    except Exception:
        pass

    if query_consumed is None:
        # fallback: be conservative — require matches be at least 90% and alignment length
        # not too short.
        if not hasattr(aln, 'matches'):
            return None
        identity = aln.matches / qlen
        if identity < min_identity or (len(analysis["insertions"]) + len(analysis['deletions']) + len(analysis['mismatches'])) > (len(query) // 12):
            return None
        # Accept but mark that we couldn't confirm full query consumption exactly.
        # Compute approximate ref coordinates: the alignment was done against ref_subseq starting at window_start.
        # The best alignment end position can sometimes be taken from aln.end_ref if present.
        aln_ref_start = window_start
        aln_ref_end = window_start + (getattr(aln, 'end_ref', 0) if hasattr(aln, 'end_ref') else len(ref_subseq)) - 1
        return {
            'score': aln.score,
            'matches': getattr(aln, 'matches', None),
            'aligned_length': getattr(aln, 'length', None),
            'identity': round(identity, 4),
            'ref_start': analysis['beg_ref'],
            'ref_end': aln_ref_end,
            'cigar': cigar_str,
            'n_seeds': cluster.get('n_seeds', 0),
            'mismatches': analysis['mismatches'],
            'insertions': analysis['insertions'],
            'deletions': analysis['deletions'],
            'analysis_end': analysis['end_ref']
        }

    # If we get here, we have query_consumed (an integer)
    if query_consumed != qlen:
        # query not fully consumed (some of query was clipped) -> reject
        return None

    # compute identity and check threshold
    matches = getattr(aln, 'matches', None)
    if matches is None:
        # cannot compute identity -> reject
        return None
    identity = matches / qlen
    if identity < min_identity:
        return None
    
    if len(analysis["insertions"]) + len(analysis['deletions']) + len(analysis['mismatches']) >  len(query) // 12:
        return None

    # compute ref coordinate mapping: if we have aln_ref_start/ end (relative to ref_subseq),
    # convert to reference coordinates
    if aln_ref_start is None:
        # try to get end_ref attribute from result object (some variants provide it)
        ref_end_attr = getattr(aln, 'end_ref', None)
        ref_begin_attr = getattr(aln, 'begin_ref', None)
        if ref_begin_attr is not None:
            aln_ref_start = window_start + ref_begin_attr
        elif ref_end_attr is not None and hasattr(aln, 'length'):
            # approximate start
            aln_ref_end = window_start + ref_end_attr
            aln_ref_start = max(0, aln_ref_end - aln.length + 1)
        else:
            aln_ref_start = window_start

    if aln_ref_end is None:
        # try end_ref
        end_ref = getattr(aln, 'end_ref', None)
        if end_ref is not None:
            aln_ref_end = window_start + end_ref
        else:
            # best-effort: compute from start + alignment length
            aln_ref_end = min(len(reference) - 1, aln_ref_start + getattr(aln, 'length', qlen) - 1)

    return {
        'score': aln.score,
        'matches': matches,
        'aligned_length': getattr(aln, 'length', None),
        'identity': round(identity, 4),
        'ref_start': analysis['beg_ref'],
        'ref_end': aln_ref_end,
        #'ref_end': analysis['end_ref'],
        'cigar': cigar_str,
        'n_seeds': cluster.get('n_seeds', 0),
        'mismatches': analysis['mismatches'],
        'insertions': analysis['insertions'],
        'deletions': analysis['deletions'],
        'analysis_end': analysis['end_ref']
    }

# consider putting this function inside of the filter_redundant_alignments_by_error function
def is_subset_with_tolerance(aln_a, aln_b, tol=0, start_diff_tolerance = 10):
    """
    Return True if all errors in errors_a are also in errors_b within tolerance.
    """

    start_diff = (aln_a['absolute_ref_start']) - (aln_b['absolute_ref_start'])
    if abs(start_diff) > 10:
        return False

    deletions_a = {x[0] for x in aln_a['deletions']} # makes a set
    insertions_a = {x[0] for x in aln_a['insertions']} # makes a set
    mismatches_a = {x[0] for x in aln_a['mismatches']} # makes a set

    for i in aln_b['insertions']:
        #if i[0] + start_diff not in insertions_a:
        if i[0] not in insertions_a:
            return False

    for d in aln_b['deletions']:
        #if d[0] + start_diff not in deletions_a:
        if d[0] not in deletions_a:
            return False

    for m in aln_b['mismatches']:
        #if m[0] + start_diff not in mismatches_a:
        if m[0] not in mismatches_a:
            return False
    
    return True


def filter_redundant_alignments_by_errors(alignments, tol=0):
    """
    Remove redundant alignments where one alignment's mistakes are a subset of another's.
    
    Args:
        alignments (list of dict): each dict has:
            - "score"
            - "insertions": list[int]
            - "deletions": list[int]
            - "mismatches": list[int]
        tol (int): position tolerance for matching mistakes
    
    Returns:
        list of dict
    """
    kept = []
    added_start_pos = set()
    alignments = sorted(alignments, key=lambda x: x["score"], reverse=True)

    for i, aln_a in enumerate(alignments):
        if aln_a['absolute_ref_start'] in added_start_pos:
            continue
        redundant = False
        #errors_a = aln_a["insertions"] + aln_a["deletions"] + aln_a["mismatches"]

        for aln_b in kept:
            #errors_b = aln_b["insertions"] + aln_b["deletions"] + aln_b["mismatches"]

            # If A’s errors are covered by B’s errors, then A is redundant (worse).
            if is_subset_with_tolerance(aln_a, aln_b, tol=tol):
                redundant = True
                break
        
        if not redundant:
            kept.append(aln_a)
            added_start_pos.add(aln_a['absolute_ref_start'])


    return kept

def best_nonoverlapping_alignments(alignments, reference_length, last_elem, reference_start=0):
    """
    Select the best non-overlapping subset of alignments with priorities:
      1) maximize count
      2) minimize number of gaps:
           - gap between consecutive alignments if end_j < start_i
           - +1 if first chosen start > reference_start
           - +1 if last  chosen end   < reference_length
      3) maximize total score

    Args
    ----
    alignments : list[tuple]
        Each item is (start, end, score, *rest). Coordinates are 0-based half-open [start, end).
    reference_length : int
        End coordinate of the reference (i.e., last valid index + 1).
    reference_start : int, default 0
        Start coordinate of the reference.

    Returns
    -------
    result : dict
        {
          "count": int,
          "gaps": int,             # includes start/end gaps per your definition
          "score": float,
          "path": list[tuple],     # chosen alignments in increasing order of end
        }
    """
    if not alignments:
        # With no alignments, by your definition there are two gaps: start and end
        return []

    # sort by end (classic interval scheduling order)
    aligns = sorted(alignments, key=lambda x: x['absolute_ref_end'])
    n = len(aligns)

    # dp[i] stores the best solution that ENDS at alignment i
    # tuple layout:
    #   (count, gaps_including_start, score, first_start, path_list)
    dp = [None] * n

    #for i, (si, ei, sco_i, *rest_i) in enumerate(aligns):
    for i in range(n):
        curr_align = aligns[i]
        si = curr_align['absolute_ref_start']
        ei = curr_align['absolute_ref_end']
        sco_i = curr_align['score']
        # start a chain with only alignment i
        start_gap = 0 if si == reference_start else 1
        best = (1, start_gap, sco_i, si, [aligns[i]])

        # try to extend from any non-overlapping j < i
        for j in range(i):
            #sj, ej, sco_j, *rest_j = aligns[j]
            curr_cmp_align = aligns[j]
            sj = curr_cmp_align['absolute_ref_start']
            ej = curr_cmp_align['absolute_ref_end']
            sco_j = curr_cmp_align['score']
            #if ej <= si:  # non-overlapping (abutting allowed)
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
        last_end = path[-1]['absolute_ref_end']  # end of the last chosen alignment
        end_gap = 0 if last_end+1 == reference_length or last_elem else 1
        total_gaps = gaps_w_start + end_gap

        key = (cnt, -total_gaps, sc)  # count ↑, gaps ↓, score ↑
        if best_overall is None or key > best_key:
            best_overall = {
                "count": cnt,
                "gaps": total_gaps,
                "score": sc,
                "path": path
            }
            best_key = key

    return best_overall['path']


def seed_and_extend_pipeline(query,
                             reference,
                             last_elem=False,
                             k=11,
                             flank=60,
                             match=2,
                             mismatch=-2,
                             gap_open=5,
                             gap_extend=1,
                             min_identity=0.90,
                             offset_tolerance=6,
                             min_seeds=2):
    """
    Full pipeline:
      - build index
      - find seeds
      - cluster seeds
      - extend clusters and filter by identity & full query coverage
    Returns list of alignment dicts.
    """
    # normalize (uppercase)
    #print("_____________________________________________________________")
    query = query.upper()
    reference = reference.upper()

    index = build_kmer_index(reference, k)
    hits = find_seed_hits(query, index, k)
    if not hits:
        return []

    clusters = cluster_hits_by_offset(hits,
                                      offset_tolerance=offset_tolerance,
                                      min_seeds=min_seeds)


    if not clusters:
        return []

    results = []
    all_used_regions = set()
    for cl in clusters:
        alns, used_regions = extend_cluster(query,
                             reference,
                             cl,
                             flank=flank,
                             window_size=150,
                             step_size=1, #modify step size
                             match=match,
                             mismatch=mismatch,
                             gap_open=gap_open,
                             gap_extend=gap_extend,
                             min_identity=min_identity,
                             used_regions=all_used_regions)

        #if used_regions:
            #for ur in used_regions:
                #all_used_regions.append(ur)

        if alns:
            for aln in alns:
                # attach cluster info
                aln['cluster_offset'] = cl['offset']
                aln['cluster_qmin'] = cl['qmin']
                aln['cluster_qmax'] = cl['qmax']
                aln['cluster_rmin'] = cl['rmin']
                aln['cluster_rmax'] = cl['rmax']
                results.append(aln)


    results = sorted(results, key=lambda x: x["score"], reverse=True)
    #print(f"alns: {results}")
    # sort by identity then score
    #results.sort(key=lambda x: (x['identity'], x['score']), reverse=True)
    #note: should probably integerate filtering into extension function
    results = filter_redundant_alignments_by_errors(results, 0)
    final_results = best_nonoverlapping_alignments(results, len(reference), last_elem, 0)
    #selected_set = weighted_interval_scheduling(results)


    return final_results
    
    return selected_set


# Example usage:
if __name__ == "__main__":
    query = "TGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTG"
    #query = "TGTGGTGTGGGTGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGGG"
    # small toy reference - replace with your actual ~5k sequence
    reference = "GTGTGTGTGGTGTGTGGGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGTGTGTGTGTGGGTGTGTGGGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTG"
    #reference = "GTGGTGTGGTGTGTGGGTGTGTGGGTGTGGGTGTGGTGTGGATGTGGTGTGGGTGTGGATGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGGTGTGGTGTGTGTGGGTGTGGTGTGTGTGGGTGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGGTGTGTGGTGTGTGGGTGTGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGGTGTGGTGTGTGGGGTGTGGTGTGGGTGTGTGTGGGTGTGGTGTGTGGGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGGGTGTGGTGTGGGTGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGTGTGTGGGTGTGGTGTGTGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGTGTGGTGTGGGTGTGGGTGTGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGGGTGTGGTGTGGGTGTGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGGGTGTGGTGTGGGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGTGTGGTGTGTGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGGGTGTGGTGTGTGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGGGTGTGGTGTGTGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGTGTGTGGGTGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGTGTGGTGTGGGTGTGGTGTGTGGTGTGTGGGTGGTGTGTGGGTGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGTGTGGGGTGTGGTGTGGGTGTGTGTGGGTGTGGTGTGTGGGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGGGTGTGTGGGGTGTGGTGTGGGTGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGGGTGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGGGTGTGGTGTGTGGTGTGTGGTGTGGGTGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGGTGTGGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGGGTGTGGTGTGGGTGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGGGTGTGGTGTGGGTGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGGGTGTGGT"

    res = seed_and_extend_pipeline(query, reference,
                                   k=15,
                                   flank=60,
                                   match=2,
                                   mismatch=-2,
                                   gap_open=10,
                                   gap_extend=3,
                                   min_identity=0.90,
                                   offset_tolerance=6,
                                   min_seeds=2)

    #res = iterative_mask_and_search(query, reference, k=11, flank=80,
                                        #min_identity=0.90, offset_tolerance=8,
                                        #min_seeds=2, score_threshold=None)

    print("Found alignments:", len(res))
    for a in res:
        print(a)
