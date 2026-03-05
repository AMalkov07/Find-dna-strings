from utils.data_structures import TelomereSequence
from typing import Dict, List, Tuple, Optional
from concurrent.futures import ProcessPoolExecutor
from Bio import SeqIO
from Bio.Seq import Seq
import re


def _process_read_ends(args: Tuple) -> Tuple[str, Optional[str], Optional[str]]:
    """
    Module-level worker for ProcessPoolExecutor — extract telomeres from both ends of one read.
    Must be module-level (not a method) so multiprocessing can pickle it.
    Returns (header, start_telomer, end_telomer).
    """
    import re as _re
    from Bio.Seq import Seq as _Seq

    header, sequence, threshold = args

    def ac_tg_fraction(subseq):
        ac = sum(1 for b in subseq if b in ('A', 'C'))
        tg = sum(1 for b in subseq if b in ('T', 'G'))
        if ac >= tg:
            return ac / len(subseq), "AC"
        return tg / len(subseq), "TG"

    def count_opposite(subseq, dom):
        bad = ('T', 'G') if dom == "AC" else ('A', 'C')
        return sum(1 for b in subseq if b in bad)

    def first_opposite_idx(subseq, dom):
        bad = ('T', 'G') if dom == "AC" else ('A', 'C')
        for i, b in enumerate(subseq):
            if b in bad:
                return i
        return None

    def extract_telomer(seq, window_size=40, ignore_last=10, opp_stop_win=15, opp_stop_count=3):
        seq = seq.upper()
        n = len(seq)
        if n < window_size + ignore_last:
            return None
        detect_start = max(0, n - ignore_last - window_size)
        detect_end = n - ignore_last
        frac, dominant_type = ac_tg_fraction(seq[detect_start:detect_end])
        if frac < threshold:
            return None
        start = detect_start
        while start > 0:
            cs = start - 1
            window = seq[cs:cs + window_size]
            if len(window) < window_size:
                break
            frac, _ = ac_tg_fraction(window)
            if frac < threshold:
                break
            opp_win = seq[cs:cs + opp_stop_win]
            if count_opposite(opp_win, dominant_type) >= opp_stop_count:
                rev_idx = first_opposite_idx(opp_win[::-1], dominant_type)
                if rev_idx is not None:
                    start = cs + (len(opp_win) - rev_idx)
                break
            start = cs
        end = detect_end
        while end < n:
            ws = max(end - window_size, start)
            window = seq[ws:end]
            if len(window) < window_size:
                window = seq[max(0, n - window_size):n]
            frac, _ = ac_tg_fraction(window)
            if frac < threshold:
                break
            ows = max(end - opp_stop_win, start)
            opp_win = seq[ows:end]
            if count_opposite(opp_win, dominant_type) >= opp_stop_count:
                idx = first_opposite_idx(opp_win, dominant_type)
                if idx is not None:
                    end = ows + idx
                break
            end += 1
        out = seq[start:end]
        if dominant_type == "AC":
            out = str(_Seq(out).reverse_complement())[::-1]
        return out

    start_tel = extract_telomer(sequence)
    if start_tel and _re.search(r'(.)\1{19,}', start_tel):
        start_tel = None
    end_tel = extract_telomer(sequence[::-1])
    if end_tel and _re.search(r'(.)\1{19,}', end_tel):
        end_tel = None
    return header, start_tel, end_tel

class FastaReader:
    def __init__(self, fasta_file_path: str, max_ends: int):
        self.fasta_file_path = fasta_file_path
        self.max_ends = max_ends
    
    def _read_fasta(self) -> List[Tuple[str, str]]:
        fasta_file_path = self.fasta_file_path
    
        records = list(SeqIO.parse(fasta_file_path, "fasta"))
    
        # Check if we actually parsed any records
        if not records:
            raise ValueError("No valid FASTA records found in file.")
    
        # Extract headers and sequences
        results = [(record.id, str(record.seq)) for record in records]
        return results

    def _extract_header_info(self, header: str) -> Tuple[str, str, int]:
        #fix: add user input option for this regular expression and harded if statement below
        match = re.search(r'_(\d+)([LR])', header)
        if match:
            survivor_name = match.group(0)
            index = int(match.group(1))  # Extract the matched number
            # Extract 'L' or 'R', or default to empty
            modifier = match.group(2)
            if modifier not in "LR":
                raise ValueError(f"Invalid header format: {header}")
            index *= 2
            if modifier == "L":
                index -= 1

            if 1 <= index <= self.max_ends:
                return (survivor_name, match.group(1)+match.group(2), index - 1)  # Convert to 0-based indexing
        raise ValueError(f"Invalid header format: {header}")

    def _ac_tg_fraction(self, subseq) -> Tuple[float, str]:
        ac_count = sum(1 for b in subseq if b in ('A','C'))
        tg_count = sum(1 for b in subseq if b in ('T','G'))
        if ac_count >= tg_count:
            return ac_count / len(subseq), "AC"
        else:
            return tg_count / len(subseq), "TG"

    def _count_opposite_bases(self, subseq: str, dominant_type: str) -> int:
        #Count number of bases that are opposite to the current dominant type.
        if dominant_type == "AC":
            return sum(1 for b in subseq if b in ('T', 'G'))
        else:  # dominant_type == "TG"
            return sum(1 for b in subseq if b in ('A', 'C'))

    def _find_first_opposite_index(self, subseq: str, dominant_type: str) -> Optional[int]:
        """Return index of first opposite base in subseq, or None if none found."""
        opposite_set = ('T', 'G') if dominant_type == "AC" else ('A', 'C')
        for i, b in enumerate(subseq):
            if b in opposite_set:
                return i
        return None

    #fix add user variables for most of the inputs to below function
    def _extract_telomer(self, dna_sequence: str, window_size: int = 40, threshold: float = .9, ignore_last: int = 10, opposite_stop_window: int = 15, opposite_stop_count: int = 3) -> Optional[str]:

        seq_str = str(dna_sequence).upper()
        n = len(seq_str)
        if n < window_size + ignore_last:
            return None  # Too short to evaluate

        
        # --- Step 1: Check detection window ([-50:-10])
        detect_start = max(0, n - ignore_last - window_size)
        detect_end = n - ignore_last
        detect_window = seq_str[detect_start:detect_end]

        frac, dominant_type = self._ac_tg_fraction(detect_window)
        if frac < threshold:
            return None  # No telomer detected

        # --- Step 2: Extend backward
        start = detect_start
        while start > 0:
            candidate_start = start - 1
            window = seq_str[candidate_start:candidate_start + window_size]
            if len(window) < window_size:
                break
            frac, _ = self._ac_tg_fraction(window)
            if frac < threshold:
                break

            # Check for opposite bases stop condition
            opp_window = seq_str[candidate_start:candidate_start + opposite_stop_window]
            if self._count_opposite_bases(opp_window, dominant_type) >= opposite_stop_count:
                reversed_idx = self._find_first_opposite_index(opp_window[::-1], dominant_type)
                if reversed_idx is not None:
                    start = candidate_start + (len(opp_window) - reversed_idx)
                break
            
            start = candidate_start

        # --- Step 3: Extend forward (into ignored region if possible)
        end = detect_end
        while end < n:
            window_start = max(end - window_size, start)
            window = seq_str[window_start:end]
            if len(window) < window_size:
                # Pad with earlier bases if near end
                window = seq_str[max(0, n - window_size):n]
            frac, _ = self._ac_tg_fraction(window)
            if frac < threshold:
                break

            # Check for opposite bases stop condition
            opp_window_start = max(end - opposite_stop_window, start)
            opp_window = seq_str[opp_window_start:end]
            if self._count_opposite_bases(opp_window, dominant_type) >= opposite_stop_count:
                idx = self._find_first_opposite_index(opp_window, dominant_type)
                if idx is not None:
                    end = opp_window_start + idx  # stop just before first opposite base
                break


            end += 1

        output = dna_sequence[start:end]
        if dominant_type == "AC":
            output = str(Seq(output).reverse_complement())
            output = output[::-1] # <<<<<<<<<< important line about whether or not to reverse an ac string after getting the reverse compliment
        return output

    def _get_better_end(self, start_telomer: Optional[str], end_telomer: Optional[str]) -> Optional[str]:
        if not start_telomer and not end_telomer:
            return None
        if not start_telomer:
            return end_telomer
        if not end_telomer:
            return start_telomer
        if len(start_telomer) >= len(end_telomer):
            return start_telomer
        return end_telomer
    
        
    def parse_fasta_population(self, threshold: float = 0.75,
                               max_workers: Optional[int] = None) -> Tuple[Dict[str, str], int]:
        """
        Population mode: check both ends of each read for telomeric sequence.
        Returns (telomere_dict, total_reads) where telomere_dict maps
        read_id -> telomere_sequence (all in TG orientation).
        Reads with telomeres on both ends produce two entries (_start / _end).

        threshold: minimum TG or AC fraction to qualify as telomeric.
          Use a lower value than the default 0.9 for raw nanopore reads
          since sequencing errors reduce the apparent fraction.
        """
        fasta_extract = self._read_fasta()
        telomeres: Dict[str, str] = {}

        args = [(header, sequence, threshold) for header, sequence in fasta_extract]
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            results = list(executor.map(_process_read_ends, args))

        for header, start_telomer, end_telomer in results:
            if start_telomer and end_telomer:
                telomeres[f"{header}_start"] = start_telomer
                telomeres[f"{header}_end"] = end_telomer
            elif start_telomer:
                telomeres[header] = start_telomer
            elif end_telomer:
                telomeres[header] = end_telomer
        return telomeres, len(fasta_extract)

    def parse_fasta(self) -> List[Optional[TelomereSequence]]:
        """
        Entry point: fasta file path (from constructor)
        Exit point: List of TelomereSequence objects
        """
        fasta_extract = self._read_fasta()
        telomers: List[Optional[TelomereSequence]] = [None] * self.max_ends
        for header, sequence in fasta_extract:
            survivor_name, chr_match, chr_index_converted = self._extract_header_info(header)
            if telomers[chr_index_converted]:
                raise ValueError(f"The FASTA file contains the same {chr_match} end multiple times")
            start_telomer: Optional[str] = self._extract_telomer(sequence)
            if start_telomer and len(start_telomer) >= len(sequence) / 2:
                best_telomer: Optional[str] = start_telomer
            else:
                sequence_reverse = sequence[::-1]
                end_telomer: Optional[str] = self._extract_telomer(sequence_reverse)
                best_telomer: Optional[str] = self._get_better_end(start_telomer, end_telomer)
            telomers[chr_index_converted] = TelomereSequence(survivor_name=survivor_name,
                                                              sequence=best_telomer,
                                                              chromosome_end_id=chr_match,
                                                              analysis=None)

        return telomers
    