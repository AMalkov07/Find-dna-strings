#!/usr/bin/env python3
"""
parasail_find_match_sg.py

Find the best semi-global alignment (entire query must be aligned) of a short query
against a long reference using parasail. Use a *_trace_* semi-global aligner so we get
a CIGAR and can enumerate matches, mismatches, insertions, and deletions and their
locations on the reference.

Positions printed are 0-based by default; pass --one-based to print 1-based positions.

Requires: parasail (pip install parasail)
"""

import re
import argparse
import sys
from cigar import Cigar

try:
    import parasail
except Exception as e:
    raise RuntimeError("parasail library not found. Install with `pip install parasail`.") from e


CIGAR_OP_RE = re.compile(r'(\d+)([MIDNSHP=X])')

# ---------------------------
# Utility / parsing functions
# ---------------------------

def parse_cigar_string(cigar_str):
    """Return list of (op, length) pairs from a CIGAR string like '8M2I3D'."""
    ops = []
    for m in CIGAR_OP_RE.finditer(cigar_str):
        l = int(m.group(1))
        op = m.group(2)
        ops.append((op, l))
    #print(f"ops: {ops}")
    return ops


def try_get_cigar_string(result):
    """
    Try several ways to obtain a CIGAR string from a parasail result.

    Many parasail Python bindings expose methods differently; this function tries:
      1) parasail.cigar_decode(result.cigar)
      2) result.cigar.decode or str(result.cigar)
      3) result.traceback() or result.get_traceback() returning aligned strings
         and builds a CIGAR from the aligned strings.

    Returns a CIGAR string (like '5M1I4M') or None if not found.
    """
    # 1) parasail helper
    try:
        if hasattr(parasail, 'cigar_decode') and hasattr(result, 'cigar') and result.cigar is not None:
            s = parasail.cigar_decode(result.cigar)
            if isinstance(s, str) and s:
                return s
    except Exception:
        pass

    # 2) cigar object decode or str()
    try:
        cigar_obj = getattr(result, 'cigar', None)
        if cigar_obj is not None:
            if hasattr(cigar_obj, 'decode'):
                dec = cigar_obj.decode
                if callable(dec):
                    return dec()
                else:
                    return str(dec)
            # fallback to str(...) if that includes M/I/D ops
            s = str(cigar_obj)
            if s and any(ch in s for ch in 'MID'):
                return s
    except Exception:
        pass

    # 3) try getting aligned strings and rebuild CIGAR
    try:
        tb_func = None
        if hasattr(result, 'traceback') and callable(getattr(result, 'traceback')):
            tb_func = result.traceback
        elif hasattr(result, 'get_traceback') and callable(getattr(result, 'get_traceback')):
            tb_func = result.get_traceback

        if tb_func is not None:
            tb = tb_func()
            # common return shape: (ref_aln, comp_aln, query_aln) or (ref_aln, query_aln)
            if isinstance(tb, (list, tuple)) and len(tb) >= 2:
                ref_aln = tb[0]
                query_aln = tb[-1]
                # build cigar
                cigar = []
                i = 0
                while i < len(ref_aln):
                    if ref_aln[i] != '-' and query_aln[i] != '-':
                        j = i
                        while j < len(ref_aln) and ref_aln[j] != '-' and query_aln[j] != '-':
                            j += 1
                        cigar.append((j - i, 'M'))
                        i = j
                    elif ref_aln[i] == '-' and query_aln[i] != '-':
                        j = i
                        while j < len(ref_aln) and ref_aln[j] == '-':
                            j += 1
                        cigar.append((j - i, 'I'))
                        i = j
                    else:
                        j = i
                        while j < len(ref_aln) and query_aln[j] == '-':
                            j += 1
                        cigar.append((j - i, 'D'))
                        i = j
                return ''.join(f"{l}{op}" for l, op in cigar)
    except Exception:
        pass

    return None


# --------------------------------
# Analyze alignment from a CIGAR
# --------------------------------
def analyze_alignment_from_cigar(query, ref, cigar_str,
                                 beg_ref=None, beg_query=None,
                                 end_ref=None, end_query=None,
                                 clip_ref_edges=True):
    """
    Robust CIGAR -> coordinate analyzer for semi-global alignments.

    Parameters
    ----------
    query, ref : str
        Original sequences.
    cigar_str : str
        Textual CIGAR (e.g. "5M1I4M").
    beg_ref, end_ref, beg_query, end_query : int or None
        If available, coordinates from parasail result. end_ref/beg_ref are
        assumed to be 0-based indexes (but function tries both inclusive/exclusive).
    clip_ref_edges : bool
        If True, DROP or trim deletions/insertions/mismatches that lie entirely
        outside the span defined by the first and last exact *match*. That
        implements your "clip edges of the reference alignment" request.

    Returns
    -------
    dict with keys:
      'beg_ref','end_ref','beg_query','end_query',
      'matches','mismatches','insertions','deletions','cigar'
    - mismatches: list of (ref_idx, query_idx, ref_base, query_base)
    - insertions: list of (ref_pos, query_pos, length, seq)  (ref_pos is where insertion occurs)
    - deletions: list of (ref_pos_start, length, seq)
    """

    ops = parse_cigar_string(cigar_str)
    if len(ops) == 0:
        raise ValueError("Empty CIGAR")

    # count consumption
    ref_consumed = sum(l for (op, l) in ops if op in ('M', 'D', '=', 'X'))
    query_consumed = sum(l for (op, l) in ops if op in ('M', 'I', '=', 'X'))
    # soft clip lengths (if any)
    left_soft = ops[0][1] if ops[0][0] == 'S' else 0
    right_soft = ops[-1][1] if ops[-1][0] == 'S' else 0
    query_total = query_consumed + left_soft + right_soft

    # infer beg_query / end_query if not provided
    if beg_query is None:
        # prefer explicit beg_query if user passed; else infer from leading S
        beg_query = left_soft
    if end_query is None:
        end_query = beg_query + query_consumed + right_soft - 1

    # If user expects entire query aligned (semi-global), validate this:
    if query_total != len(query):
        # Not fatal, but warn via exception text so caller notices.
        # (We could also choose to silently continue; you asked for clarity.)
        raise ValueError(
            f"CIGAR-consumed query length ({query_total}) != len(query) ({len(query)}). "
            "For a true semi-global alignment the entire query should be consumed. "
            "Check that you called an sg_trace_* function (no soft-clips) or that your "
            "CIGAR matches expectations."
        )

    # Compute beginning ref coordinate if missing, try inclusive/exclusive interpretations
    if ops[0][0] == 'D':
        beg_ref = ops[0][1]
    if beg_ref is None:
        beg_ref = 0

    # Now do the actual mapping walk; check bounds as we go
    ref_pos = int(beg_ref)
    query_pos = int(beg_query)

    matches = 0
    mismatches = []
    insertions = []  # (ref_pos_where_insertion_happens, query_pos_start, length, seq)
    deletions = []   # (ref_pos_start, length, seq)
    match_positions_ref = []  # all ref indices that are exact matches

    if ops[0][0] == 'D':
        ops = ops[1:]

    op_save = None
    L_save = None

    for op, L in ops:
        if query_pos >= len(query) or ref_pos >= len(ref):
            op_save = op
            L_save = L
            break
        if op in ('M', '=', 'X'):
            for i in range(L):
                #r_idx = ref_pos + i - beg_ref
                #q_idx = query_pos + i - beg_query
                r_idx = ref_pos + i 
                q_idx = query_pos + i
                #r_idx = i + L
                #q_idx = i + L
                if r_idx < 0 or r_idx >= len(ref) or q_idx < 0 or q_idx >= len(query):
                    raise IndexError(f"Index out of range during CIGAR walk: r_idx={r_idx}, q_idx={q_idx}")
                rbase = ref[r_idx]
                qbase = query[q_idx]
                if rbase == qbase:
                    matches += 1
                    match_positions_ref.append(r_idx)
                else:
                    #mismatches.append((r_idx, q_idx, rbase, qbase))
                    mismatches.append((q_idx+1, qbase, rbase)) # the relative r_idx should be r_idx - true_beg_ref
                    #mismatches.append((r_idx-true_beg_ref, q_idx-true_beg_ref, rbase, qbase))
            ref_pos += L
            query_pos += L
        elif op == 'I': #treat I's like deletions not insertions
            # insertion relative to reference (consumes query only)
            # convention: insertion at current ref_pos (between ref_pos-1 and ref_pos)
            for i in range(L):
                q_idx = query_pos + i
                q_base = query[q_idx]
                deletions.append((q_idx+1, q_base))
            #insertions.append((query_pos, ref_pos, L, seq))
            #insertions.append((ref_pos-true_beg_ref, query_pos-true_beg_ref, L, seq))
            query_pos += L
            #ref_pos += L
        elif op == 'D': # treat D's like insertions not deletions
            for i in range(L):
                r_idx = ref_pos + i
                r_base = ref[r_idx]
                insertions.append((query_pos, r_base))
            ref_pos += L
            #query_pos += 1
        elif op == 'S':
            # soft clip: consumes query only (these were accounted for earlier)
            query_pos += L
        elif op == 'H':
            # hard clip: consumes neither
            pass
        elif op == 'N':
            # skipped region in reference (like intron) - consumes ref
            ref_pos += L
        else:
            raise ValueError(f"Unhandled CIGAR op '{op}'")

    #incase reference is smaller then query and there are a bunch of deletions at the end:
    if query_pos < len(query) and op_save == 'I' and L_save:
        for i in range(L_save):
            q_idx = query_pos + i
            q_base = query[q_idx]
            deletions.append((q_idx+1, q_base))
        

    #computed_end_ref = ref_pos - 1
    #computed_end_query = query_pos - 1

    ## If caller provided end_ref / end_query, prefer them (but computed ones are returned if caller omitted)
    #if end_ref is None:
        #end_ref = computed_end_ref
    #if end_query is None:
        #end_query = computed_end_query

    ## If clipping of ref edges requested, compute first/last exact match and filter events outside that span
    #if clip_ref_edges and match_positions_ref:
        #first_match_ref = min(match_positions_ref)
        #last_match_ref = max(match_positions_ref)

        ## filter mismatches to those inside [first_match_ref, last_match_ref]
        #mismatches = [mm for mm in mismatches if first_match_ref <= mm[0] <= last_match_ref]

        ## filter/clip insertions: insertion at ref_pos is kept only if ref_pos is within [first,last]
        #insertions = [ins for ins in insertions if first_match_ref <= ins[0] <= last_match_ref]

        ## trim deletions that partially overlap the match span
        #trimmed_deletions = []
        #for (rpos, L, seq) in deletions:
            #del_start = rpos
            #del_end = rpos + L - 1
            #overlap_start = max(del_start, first_match_ref)
            #overlap_end = min(del_end, last_match_ref)
            #if overlap_start <= overlap_end:
                ## compute trimmed indices relative to original deletion
                #trim_left = overlap_start - del_start
                #trim_right = overlap_end - del_start  # inclusive
                #trimmed_seq = seq[trim_left:trim_right+1]
                ##trimmed_deletions.append([overlap_start, len(trimmed_seq), trimmed_seq])
                #trimmed_deletions.append([overlap_start-true_beg_ref, len(trimmed_seq), trimmed_seq])
            ## else — deletion entirely outside match interval -> drop
        #deletions = trimmed_deletions
    #else:
        ## if there were no exact matches and clip_ref_edges requested,
        ## we keep everything (alternatively we could choose to drop edge-deletions,
        ## but it's ambiguous when there are no exact matches).
        #pass

    end_ref = ref_pos - 1

    return {
        #'beg_ref': int(beg_ref),
        'beg_ref': int(beg_ref),
        'end_ref': int(end_ref),
        'beg_query': int(beg_query),
        'end_query': int(end_query),
        'matches': int(matches),
        'mismatches': mismatches,
        'insertions': insertions,
        'deletions': deletions,
        'cigar': cigar_str,
    }

##def analyze_alignment_from_cigar(query, ref, cigar_str,
                                 ##beg_ref=None, beg_query=None,
                                 ##end_ref=None, end_query=None,
                                 ##clip_ref_edges=True):
    ##"""
    ##Robust CIGAR -> coordinate analyzer for semi-global alignments.

    ##Parameters
    ##----------
    ##query, ref : str
        ##Original sequences.
    ##cigar_str : str
        ##Textual CIGAR (e.g. "5M1I4M").
    ##beg_ref, end_ref, beg_query, end_query : int or None
        ##If available, coordinates from parasail result. end_ref/beg_ref are
        ##assumed to be 0-based indexes (but function tries both inclusive/exclusive).
    ##clip_ref_edges : bool
        ##If True, DROP or trim deletions/insertions/mismatches that lie entirely
        ##outside the span defined by the first and last exact *match*. That
        ##implements your "clip edges of the reference alignment" request.

    ##Returns
    ##-------
    ##dict with keys:
      ##'beg_ref','end_ref','beg_query','end_query',
      ##'matches','mismatches','insertions','deletions','cigar'
    ##- mismatches: list of (ref_idx, query_idx, ref_base, query_base)
    ##- insertions: list of (ref_pos, query_pos, length, seq)  (ref_pos is where insertion occurs)
    ##- deletions: list of (ref_pos_start, length, seq)
    ##"""

    ##ops = parse_cigar_string(cigar_str)
    ##if len(ops) == 0:
        ##raise ValueError("Empty CIGAR")

    ### count consumption
    ##ref_consumed = sum(l for (op, l) in ops if op in ('M', 'D', '=', 'X'))
    ##query_consumed = sum(l for (op, l) in ops if op in ('M', 'I', '=', 'X'))
    ### soft clip lengths (if any)
    ##left_soft = ops[0][1] if ops[0][0] == 'S' else 0
    ##right_soft = ops[-1][1] if ops[-1][0] == 'S' else 0
    ##query_total = query_consumed + left_soft + right_soft

    ### infer beg_query / end_query if not provided
    ##if beg_query is None:
        ### prefer explicit beg_query if user passed; else infer from leading S
        ##beg_query = left_soft
    ##if end_query is None:
        ##end_query = beg_query + query_consumed + right_soft - 1

    ### If user expects entire query aligned (semi-global), validate this:
    ##if query_total != len(query):
        ### Not fatal, but warn via exception text so caller notices.
        ### (We could also choose to silently continue; you asked for clarity.)
        ##raise ValueError(
            ##f"CIGAR-consumed query length ({query_total}) != len(query) ({len(query)}). "
            ##"For a true semi-global alignment the entire query should be consumed. "
            ##"Check that you called an sg_trace_* function (no soft-clips) or that your "
            ##"CIGAR matches expectations."
        ##)

    ### Compute beginning ref coordinate if missing, try inclusive/exclusive interpretations
    ##if beg_ref is None:
        ###if end_ref is None:
            ###raise ValueError("Need beg_ref OR end_ref from parasail result to compute reference coordinates.")
        #### Candidate 1: end_ref is inclusive index of last consumed ref position (common)
        ###cand1 = int(end_ref) - ref_consumed + 1
        #### Candidate 2: end_ref may be exclusive (one-past-end)
        ###cand2 = int(end_ref) - ref_consumed
        ###chosen = None
        ###for cand in (cand1, cand2):
            ###if cand is None:
                ###continue
            ###if cand >= 0 and cand + ref_consumed - 1 < len(ref):
                ###chosen = int(cand)
                ###break
        ###if chosen is None:
            ###raise ValueError(
                ###"Unable to infer beg_ref from end_ref + CIGAR. "
                ###f"end_ref={end_ref}, ref_consumed={ref_consumed}, ref_len={len(ref)}. "
                ###"Check parasail's beg/end coordinate semantics for your binding."
            ###)
        ##beg_ref = 0
    ##true_beg_ref = 0
    ##if ops[0][0] == 'D':
        ##true_beg_ref = ops[0][1]
        

    ### Now do the actual mapping walk; check bounds as we go
    ##ref_pos = int(beg_ref)
    ##query_pos = int(beg_query)

    ##matches = 0
    ##mismatches = []
    ##insertions = []  # (ref_pos_where_insertion_happens, query_pos_start, length, seq)
    ##deletions = []   # (ref_pos_start, length, seq)
    ##match_positions_ref = []  # all ref indices that are exact matches

    ##for op, L in ops:
        ##if op in ('M', '=', 'X'):
            ##for i in range(L):
                ###r_idx = ref_pos + i - beg_ref
                ###q_idx = query_pos + i - beg_query
                ##r_idx = ref_pos + i 
                ##q_idx = query_pos + i
                ###r_idx = i + L
                ###q_idx = i + L
                ##if r_idx < 0 or r_idx >= len(ref) or q_idx < 0 or q_idx >= len(query):
                    ##raise IndexError(f"Index out of range during CIGAR walk: r_idx={r_idx}, q_idx={q_idx}")
                ##rbase = ref[r_idx]
                ##qbase = query[q_idx]
                ##if rbase == qbase:
                    ##matches += 1
                    ##match_positions_ref.append(r_idx)
                ##else:
                    ###mismatches.append((r_idx, q_idx, rbase, qbase))
                    ##mismatches.append((q_idx, r_idx, qbase, rbase)) # the relative r_idx should be r_idx - true_beg_ref
                    ###mismatches.append((r_idx-true_beg_ref, q_idx-true_beg_ref, rbase, qbase))
            ##ref_pos += L
            ##query_pos += L
        ##elif op == 'I':
            ### insertion relative to reference (consumes query only)
            ### convention: insertion at current ref_pos (between ref_pos-1 and ref_pos)
            ##seq = query[query_pos:query_pos+L]
            ##insertions.append((ref_pos, query_pos, L, seq))
            ###insertions.append((query_pos, ref_pos, L, seq))
            ###insertions.append((ref_pos-true_beg_ref, query_pos-true_beg_ref, L, seq))
            ##query_pos += L
        ##elif op == 'D':
            ### deletion relative to reference (consumes ref only)
            ##seq = ref[ref_pos:ref_pos+L]
            ##deletions.append((ref_pos, L, seq))
            ###deletions.append((ref_pos-true_beg_ref, L, seq))
            ##ref_pos += L
        ##elif op == 'S':
            ### soft clip: consumes query only (these were accounted for earlier)
            ##query_pos += L
        ##elif op == 'H':
            ### hard clip: consumes neither
            ##pass
        ##elif op == 'N':
            ### skipped region in reference (like intron) - consumes ref
            ##ref_pos += L
        ##else:
            ##raise ValueError(f"Unhandled CIGAR op '{op}'")

    ##computed_end_ref = ref_pos - 1
    ##computed_end_query = query_pos - 1

    ### If caller provided end_ref / end_query, prefer them (but computed ones are returned if caller omitted)
    ##if end_ref is None:
        ##end_ref = computed_end_ref
    ##if end_query is None:
        ##end_query = computed_end_query

    ### If clipping of ref edges requested, compute first/last exact match and filter events outside that span
    ##if clip_ref_edges and match_positions_ref:
        ##first_match_ref = min(match_positions_ref)
        ##last_match_ref = max(match_positions_ref)

        ### filter mismatches to those inside [first_match_ref, last_match_ref]
        ##mismatches = [mm for mm in mismatches if first_match_ref <= mm[0] <= last_match_ref]

        ### filter/clip insertions: insertion at ref_pos is kept only if ref_pos is within [first,last]
        ##insertions = [ins for ins in insertions if first_match_ref <= ins[0] <= last_match_ref]

        ### trim deletions that partially overlap the match span
        ##trimmed_deletions = []
        ##for (rpos, L, seq) in deletions:
            ##del_start = rpos
            ##del_end = rpos + L - 1
            ##overlap_start = max(del_start, first_match_ref)
            ##overlap_end = min(del_end, last_match_ref)
            ##if overlap_start <= overlap_end:
                ### compute trimmed indices relative to original deletion
                ##trim_left = overlap_start - del_start
                ##trim_right = overlap_end - del_start  # inclusive
                ##trimmed_seq = seq[trim_left:trim_right+1]
                ###trimmed_deletions.append([overlap_start, len(trimmed_seq), trimmed_seq])
                ##trimmed_deletions.append([overlap_start-true_beg_ref, len(trimmed_seq), trimmed_seq])
            ### else — deletion entirely outside match interval -> drop
        ##deletions = trimmed_deletions
    ##else:
        ### if there were no exact matches and clip_ref_edges requested,
        ### we keep everything (alternatively we could choose to drop edge-deletions,
        ### but it's ambiguous when there are no exact matches).
        ##pass

    ##return {
        ###'beg_ref': int(beg_ref),
        ##'beg_ref': int(true_beg_ref),
        ##'end_ref': int(end_ref),
        ##'beg_query': int(beg_query),
        ##'end_query': int(end_query),
        ##'matches': int(matches),
        ##'mismatches': mismatches,
        ##'insertions': insertions,
        ##'deletions': deletions,
        ##'cigar': cigar_str,
    ##}


#def analyze_alignment_from_cigar(query, ref, cigar_str, beg_ref=None, beg_query=None, end_ref=None, end_query=None):
    #"""
    #Walk the CIGAR to detect matches, mismatches, insertions, deletions and their reference locations.

    #Inputs:
      #- query, ref: strings
      #- cigar_str: textual CIGAR covering the aligned region
      #- beg_ref/beg_query OR end_ref/end_query must be provided (parasail may provide beg_* or end_*)
    #Outputs: dictionary with beg_ref, end_ref, beg_query, end_query, matches, mismatches, insertions, deletions, cigar
    #"""
    #ops = parse_cigar_string(cigar_str)

    ## count how many reference/query positions the CIGAR consumes
    ##ref_consumed = sum(l for (op, l) in ops if op in ('M', 'D', '=', 'X'))
    ##query_consumed = sum(l for (op, l) in ops if op in ('M', 'I', '=', 'X'))

    #if beg_ref is None:
        #if end_ref is None:
            #raise ValueError("Need beg_ref OR end_ref to compute coordinates.")

        #beg_ref = 0
        ##if ops[0][0] == 'D':
            ##beg_ref = ops[0][1]
        ##else:
            ##beg_ref = 0
        ##beg_ref = int(end_ref) - ref_consumed + 1

    #if beg_query is None:
        #if end_query is None:
            #raise ValueError("Need beg_query OR end_query to compute coordinates.")
        ##beg_query = int(end_query) - query_consumed + 1
        #beg_query = 0

    #ref_pos = beg_ref
    #query_pos = beg_query

    #matches = 0
    #mismatches = []
    #insertions = []  # (ref_pos_where_insertion_happens, query_pos_start, length, seq)
    #deletions = []   # (ref_pos_start, length, seq)

    #for op, L in ops:
        #if op in ('M', '=', 'X'):
            #for i in range(L):
                #r_idx = ref_pos + i
                #q_idx = query_pos + i
                #rbase = ref[r_idx]
                #qbase = query[q_idx]
                #if rbase == qbase:
                    #matches += 1
                #else:
                    #mismatches.append((r_idx, q_idx, rbase, qbase))
            #ref_pos += L
            #query_pos += L
        #elif op == 'I':
            #insertions.append((ref_pos, query_pos, L, query[query_pos:query_pos+L]))
            #query_pos += L
        #elif op == 'D':
            #deletions.append((ref_pos, L, ref[ref_pos:ref_pos+L]))
            #ref_pos += L
        #elif op == 'S':
            ## soft clip: consumes query only (not typical for semi-global internal alignment)
            #query_pos += L
        #elif op == 'H':
            ## hard clip: consumes neither (ignore)
            #pass
        #elif op == 'N':
            ## skipped region in reference (consume ref)
            #ref_pos += L
        #else:
            #raise ValueError(f"Unhandled CIGAR op '{op}'")

    #return {
        #'beg_ref': beg_ref,
        #'end_ref': (ref_pos - 1),
        #'beg_query': beg_query,
        #'end_query': (query_pos - 1),
        #'matches': matches,
        #'mismatches': mismatches,
        #'insertions': insertions,
        #'deletions': deletions,
        #'cigar': cigar_str,
    #}


# --------------------------------
# Main alignment + reporting
# --------------------------------

def find_best_sg_match_and_report(query, ref, gap_open=10, gap_extend=1, match=2, mismatch=-2):
    """
    Run a parasail semi-global (sg) trace alignment to align the entire query into the reference.

    Returns an analysis dict:
      - beg_ref, end_ref, beg_query, end_query
      - matches, mismatches (list), insertions (list), deletions (list)
      - cigar, parasail_score, and the raw parasail result (for debugging)
    """
    # Try to create a small DNA matrix; fall back to blosum62 if not available.
    matrix = None
    try:
        if hasattr(parasail, 'matrix_create'):
            matrix = parasail.matrix_create("ACGT", match, mismatch)
    except Exception:
        matrix = None
    if matrix is None:
        matrix = parasail.blosum62

    # Pick a semi-global trace function. prefer 8-bit striped for speed; try fallbacks gracefully.
    res = None
    tried_names = [
        'sg_trace_scan', 'sg_trace_striped_8', 'sg_trace_scan_8', 'sg_trace_striped_16',
        'sg_trace_scan_16', 'sg_trace_striped_32'
    ]
    for name in tried_names:
        if hasattr(parasail, name):
            res = getattr(parasail, name)(query, ref, gap_open, gap_extend, matrix)
            break
    if res is None:
        # last resort: try generic sg_trace (some builds may export without suffix)
        for name in ('sg_trace', 'sg_trace_striped', 'sg_trace_scan'):
            if hasattr(parasail, name):
                # try to call with common signatures; some may require profile objects instead
                try:
                    res = getattr(parasail, name)(query, ref, gap_open, gap_extend, matrix)
                    break
                except Exception:
                    continue
    if res is None:
        raise RuntimeError("Could not find any sg_trace_* function in your parasail binding. "
                           "Inspect `dir(parasail)` to see available functions.")

    # extract coordinates and stats if available
    end_ref = getattr(res, 'end_ref', None)
    end_query = getattr(res, 'end_query', None)
    beg_ref = getattr(res, 'beg_ref', None)
    beg_query = getattr(res, 'beg_query', None)
    matches_count = getattr(res, 'matches', None)

    #print(f"beg_ref: {beg_ref}")
    #print(f"beg_query: {beg_query}")
    #print(f"end_ref: {end_ref}")
    #print(f"end_query: {end_query}")
    #print(f"matches_count: {matches_count}")

    cigar_str = try_get_cigar_string(res)
    if cigar_str is None:
        raise RuntimeError("Could not obtain CIGAR string from parasail result. Ensure you used a *_trace_* semi-global function.")

    analysis = analyze_alignment_from_cigar(
        query=query,
        ref=ref,
        cigar_str=cigar_str,
        beg_ref=beg_ref,
        beg_query=beg_query,
        end_ref=end_ref,
        end_query=end_query
    )

    if matches_count is not None:
        analysis['matches_parasail'] = int(matches_count)

    analysis['parasail_score'] = getattr(res, 'score', None)
    analysis['parasail_result_obj'] = res

    return analysis


def print_report(analysis, one_based=False):
    """Pretty-print the analysis dictionary produced by find_best_sg_match_and_report."""
    offset = 1 if one_based else 0
    print("Alignment CIGAR:", analysis['cigar'])
    #c = Cigar(analysis['cigar'])
    #print(f"cigar human readable: {list(c.items())}")
    print(f"Reference start (0-based): {analysis['beg_ref']}  end (0-based): {analysis['end_ref']}")
    if one_based:
        print(f"Reference start (1-based): {analysis['beg_ref']+offset}  end (1-based): {analysis['end_ref']+offset}")
    print("Matches (count):", analysis['matches'])
    if 'matches_parasail' in analysis:
        print("Matches (parasail reported):", analysis['matches_parasail'])
    print("Mismatches (count):", len(analysis['mismatches']))
    if analysis['mismatches']:
        print(" Mismatches (ref_idx, query_idx, ref_base, query_base):")
        for mm in analysis['mismatches']:
            r,q,rb,qb = mm
            if one_based:
                print(f"  ref:{r+offset}  query:{q+offset}   {rb}->{qb}")
            else:
                print(f"  ref:{r}  query:{q}   {rb}->{qb}")
    print("Insertions relative to reference (count):", len(analysis['insertions']))
    for ins in analysis['insertions']:
        ref_pos, qpos, L, seq = ins
        if one_based:
            print(f"  insertion at ref:{ref_pos+offset} (query {qpos+offset}-{qpos+L-1+offset}) len={L} seq={seq}")
        else:
            print(f"  insertion at ref:{ref_pos} (query {qpos}-{qpos+L-1}) len={L} seq={seq}")
    print("Deletions relative to reference (count):", len(analysis['deletions']))
    for dele in analysis['deletions']:
        rpos, L, seq = dele
        if one_based:
            print(f"  deletion at ref:{rpos+offset} len={L} seq={seq}")
        else:
            print(f"  deletion at ref:{rpos} len={L} seq={seq}")
    print("Parasail raw score:", analysis.get('parasail_score'))

    
import parasail

def recursive_find_alignments(query, ref, ref_offset=0, results=None, min_len=None):
    """
    Recursively find best semi-global alignments of `query` in `ref`.

    Parameters
    ----------
    query : str
        Query sequence.
    ref : str
        Reference sequence segment to search.
    ref_offset : int
        Offset of `ref` relative to the original reference string.
    results : list
        Accumulates alignment results (modified in-place).
    min_len : int
        Minimum reference length to attempt alignment (defaults to len(query)).

    Returns
    -------
    results : list of dict
        Each dict is the output from analyze_alignment_from_cigar() with
        coordinates mapped to original reference.
    """

    
    if results is None:
        results = []
    if min_len is None:
        min_len = len(query) - 5

    if len(ref) < min_len:
        return results  # too short to align

    # Run semi-global alignment
    aln = parasail.sg_trace_scan_16(query, ref, 10, 1, parasail.blosum62)  # adjust matrix & penalties

    # Skip if score is 0 or no useful alignment
    if not aln.score or aln.score == 0:
        print("<<<<<<<<<<<<reutnring>>>>>>>>>>>>")
        return results

    cigar_str = try_get_cigar_string(aln)

    # Analyze alignment with your function
    aln_info = analyze_alignment_from_cigar(
        query, ref,
        cigar_str,
        beg_ref=None,
        beg_query=None,
        end_ref=aln.end_ref,
        end_query=aln.end_query,
        clip_ref_edges=True
    )



    # Adjust coordinates to original reference space
    aln_info['beg_ref'] += ref_offset
    aln_info['end_ref'] += ref_offset

    counter = 0
    #for i in aln_info['insertions']:
        ##counter += len(i[3])
    #for i in aln_info['deletions']:
        ##counter += len(i[2])
    #for i in aln_info['mismatches']:
        ##counter += len(i[2])
    
    #if len(aln_info['insertions']) + len(aln_info['deletions']) + len(aln_info['mismatches']) <= len(query) // 12:
    #if counter <= len(query)//12:
    if len(aln_info['insertions']) + len(aln_info['deletions']) + len(aln_info['mismatches']) <= len(query) // 8:
        results.append(aln_info)
    else:
        return

    # Split into left and right flanks
    left_start = 0
    #left_end = aln_info['beg_ref'] - ref_offset
    #right_start = aln_info['end_ref'] - ref_offset + 1
    left_end = aln_info['beg_ref'] - ref_offset
    right_start = aln_info['end_ref'] + 1 - ref_offset
    right_end = len(ref)

    left_ref = ref[left_start:left_end]
    right_ref = ref[right_start:right_end]

    # Recurse if flanks are large enough
    if len(left_ref) >= min_len:
        recursive_find_alignments(query, left_ref, ref_offset, results, min_len)
    if len(right_ref) >= min_len:
        recursive_find_alignments(query, right_ref, ref_offset + right_start, results, min_len)

    return results


# ---------------------------
# CLI entrypoint
# ---------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find best SEMI-GLOBAL match of query inside reference using parasail (trace variant).")
    parser.add_argument('query', help='query sequence (short)')
    parser.add_argument('reference', help='reference sequence (long)')
    parser.add_argument('--gap-open', type=int, default=10, help='gap open penalty (default 10)')
    parser.add_argument('--gap-extend', type=int, default=1, help='gap extend penalty (default 1)')
    parser.add_argument('--match', type=int, default=2, help='match score for nucleotide matrix')
    parser.add_argument('--mismatch', type=int, default=-2, help='mismatch penalty for nucleotide matrix')
    parser.add_argument('--one-based', action='store_true', help='print positions 1-based instead of 0-based')
    args = parser.parse_args()

    #report = find_best_sg_match_and_report(
        #query=args.query,
        #ref=args.reference,
        #gap_open=args.gap_open,
        #gap_extend=args.gap_extend,
        #match=args.match,
        #mismatch=args.mismatch
    #)

    all_hits = recursive_find_alignments(
        query=args.query,
        ref=args.reference,
        #gap_open=args.gap_open,
        #gap_extend=args.gap_extend,
        #match=args.match,
        #mismatch=args.mismatch
    )

    #for hit in all_hits:
        
        #print_report(hit, one_based=args.one_based)
        #print("\n")
        

    print(f"total alignments: {len(all_hits)}")
    for hit in all_hits:
        print(f"{hit['beg_ref']}, insertions: {hit['insertions']}, deletions: {hit['deletions']}, mismatches: {hit['mismatches']}")

