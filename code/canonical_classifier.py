#!/usr/bin/env python3
"""
code/canonical_classifier.py

Single authoritative sequence-processing, FASTQ-validation, and classification
primitives for the open-source computational replication of SADDLE (Xie et al.,
Nature Communications 2022, 13:1881).

Implements both:
  1. Open-source reimplementation of the reported Data_process.m and
     Dimer_Identify.m procedures (28-bp adapters, 60-bp length routing,
     length > 15 filter, Smith-Waterman local alignment with NUC44 +5/-4/-8
     scoring ratio > 0.6).
  2. Standalone Python/edlib classification primitives (17-bp core TruSeq
     prefix, trimmed-prefix primer matching, edlib semi-global gapped alignment)
     requiring no external command-line aligner.
"""

import gzip
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Tuple

try:
    import edlib
except ImportError:  # pragma: no cover
    edlib = None

# Exact 28-bp adapter sequences from SADDLE Data_process.m (lines 23-24)
ORIG_ADAPTER_R1 = "AGATCGGAAGAGCACACGTCTGAACTCC"
ORIG_ADAPTER_R2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAG"

# Core 17-bp Illumina TruSeq prefixes used in the standalone Python/edlib stream
CORE_ADAPTER_R1 = "AGATCGGAAGAGCACAC"
CORE_ADAPTER_R2 = "AGATCGGAAGAGCGTCG"

_RC_TABLE = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def write_utf8_lf(path: Path, text: str) -> None:
    """Write UTF-8 text strictly with Unix LF line endings across all platforms."""
    normalized = text.replace("\r\n", "\n").replace("\r", "\n")
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(normalized)


def revcomp(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return seq.translate(_RC_TABLE)[::-1]


def normalize_mate_id(header_line: str) -> str:
    """Extract canonical read identifier from a FASTQ header, stripping /1 or /2 suffix."""
    token = header_line.split()[0]
    if token.endswith("/1") or token.endswith("/2"):
        return token[:-2]
    return token


def _open_fastq(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="ascii", errors="strict")
    return open(path, "r", encoding="ascii", errors="strict")


def iter_validated_fastq_pairs(
    r1_path: Path, r2_path: Path, max_reads: Optional[int] = None
) -> Iterator[Tuple[str, str, str, str, str]]:
    """
    Stream paired-end FASTQ records with strict structural validation:
      - Headers must begin with '@' and separator lines with '+'.
      - R1 and R2 must contain the exact same number of complete 4-line records.
      - Normalized mate identifiers must match on every read pair.
      - Sequence and quality string lengths must match for both R1 and R2.

    Yields (read_id, seq1, qual1, seq2, qual2).
    """
    rec_num = 0
    with _open_fastq(r1_path) as f1, _open_fastq(r2_path) as f2:
        while True:
            h1 = f1.readline()
            h2 = f2.readline()
            if not h1 and not h2:
                break
            if bool(h1) != bool(h2):
                raise ValueError(
                    f"FASTQ EOF desynchronization after {rec_num} records between {r1_path.name} and {r2_path.name}"
                )
            s1 = f1.readline()
            p1 = f1.readline()
            q1 = f1.readline()
            s2 = f2.readline()
            p2 = f2.readline()
            q2 = f2.readline()
            if not (s1 and p1 and q1 and s2 and p2 and q2):
                raise ValueError(f"Truncated FASTQ record at pair {rec_num + 1}")
            if not (h1.startswith("@") and h2.startswith("@")):
                raise ValueError(f"Invalid FASTQ header at pair {rec_num + 1}: {h1[:30]!r} / {h2[:30]!r}")
            if not (p1.startswith("+") and p2.startswith("+")):
                raise ValueError(f"Invalid FASTQ '+' separator at pair {rec_num + 1}")

            id1 = normalize_mate_id(h1)
            id2 = normalize_mate_id(h2)
            if id1 != id2:
                raise ValueError(f"Mismatched mate IDs at pair {rec_num + 1}: {id1!r} != {id2!r}")

            seq1 = s1.rstrip("\r\n")
            qual1 = q1.rstrip("\r\n")
            seq2 = s2.rstrip("\r\n")
            qual2 = q2.rstrip("\r\n")
            if len(seq1) != len(qual1) or len(seq2) != len(qual2):
                raise ValueError(
                    f"Sequence/quality length mismatch at pair {rec_num + 1} ({id1}): "
                    f"R1=({len(seq1)},{len(qual1)}) R2=({len(seq2)},{len(qual2)})"
                )

            rec_num += 1
            yield id1, seq1, qual1, seq2, qual2
            if max_reads is not None and rec_num >= max_reads:
                break


def trim_adapter_pair(
    seq1: str,
    qual1: str,
    seq2: str,
    qual2: str,
    adp1: str = ORIG_ADAPTER_R1,
    adp2: str = ORIG_ADAPTER_R2,
) -> Tuple[Optional[str], Optional[str], Optional[str], Optional[str], str]:
    """
    Locate adapters in R1 and R2 and trim strictly before the first match.

    Returns (t1, q1, t2, q2, status) where status is one of:
      - 'ok', 'no_adapter_both', 'no_adapter_r1_only', 'no_adapter_r2_only'
    """
    p1 = seq1.find(adp1)
    p2 = seq2.find(adp2)
    if p1 == -1 and p2 == -1:
        return None, None, None, None, "no_adapter_both"
    if p1 == -1:
        return None, None, None, None, "no_adapter_r1_only"
    if p2 == -1:
        return None, None, None, None, "no_adapter_r2_only"
    return seq1[:p1], qual1[:p1], seq2[:p2], qual2[:p2], "ok"


def qc_check_pair(
    t1: str,
    q1: str,
    t2: str,
    q2: str,
    min_len_exclusive: int = 10,
    q30_frac_exclusive: float = 0.7,
) -> str:
    """
    Apply SADDLE Data_process.m (lines 109-111) quality filters:
      - trimmed length > min_len_exclusive (default 10 bp) on both R1 and R2
      - no ambiguous 'N' bases in either read
      - > 70% of bases with Phred > 30 (ASCII > 63) on both R1 and R2
    """
    l1 = len(t1)
    l2 = len(t2)
    if l1 <= min_len_exclusive or l2 <= min_len_exclusive:
        return "qc_short_len"
    if "N" in t1 or "N" in t2 or "n" in t1 or "n" in t2:
        return "qc_ambiguous_n"
    hq1 = sum(1 for c in q1 if ord(c) > 63)
    if hq1 / l1 <= q30_frac_exclusive:
        return "qc_low_quality"
    hq2 = sum(1 for c in q2 if ord(c) > 63)
    if hq2 / l2 <= q30_frac_exclusive:
        return "qc_low_quality"
    return "ok"


def build_prefix_index(primers: List[Tuple[str, str]], k: int = 8) -> Dict[str, List[int]]:
    """Index all primers by their lowercase k-bp 5' prefix."""
    index: Dict[str, List[int]] = {}
    for idx, (_, seq) in enumerate(primers):
        key = seq[:k].lower()
        index.setdefault(key, []).append(idx)
    return index


def sw_score_nuc44(s1: str, s2: str, match: int = 5, mismatch: int = -4, gap: int = -8) -> int:
    """
    Smith-Waterman local alignment score using the documented NUC44 nucleotide
    scoring defaults (match=+5, mismatch=-4, linear gap=-8).
    """
    if s1 == s2:
        return match * len(s1)
    n = len(s2)
    prev = [0] * (n + 1)
    best = 0
    for c1 in s1:
        curr = [0] * (n + 1)
        left = 0
        for j, c2 in enumerate(s2, 1):
            diag = prev[j - 1] + (match if c1 == c2 else mismatch)
            up = prev[j] + gap
            left_gap = left + gap
            val = diag if diag > up else up
            if left_gap > val:
                val = left_gap
            if val < 0:
                val = 0
            curr[j] = val
            left = val
            if val > best:
                best = val
        prev = curr
    return best


def match_primer_sw_nuc44(
    trimmed_read_lower: str,
    candidate_indices: List[int],
    primers_lower: List[str],
    threshold: float = 0.6,
) -> Optional[int]:
    """
    Open-source reimplementation of SADDLE Dimer_Identify.m (lines 89-97):
      - requires len(primer) <= len(trimmed_read)
      - computes Smith-Waterman NUC44 score(primer, read[0:len(primer)]) / (5 * len(primer)) > 0.6
    """
    rlen = len(trimmed_read_lower)
    for pid in candidate_indices:
        pseq = primers_lower[pid]
        plen = len(pseq)
        if plen <= rlen:
            sub = trimmed_read_lower[:plen]
            if sub == pseq:
                return pid
            score = sw_score_nuc44(pseq, sub)
            if score / (5.0 * plen) > threshold:
                return pid
    return None


def match_primer_prefix(
    trimmed_read_lower: str,
    candidate_indices: List[int],
    primers_lower: List[str],
    min_prefix: int = 8,
    threshold: float = 0.6,
) -> Optional[int]:
    """
    Standalone Python/edlib primer prefix matcher.
    Evaluates positional identity over check_len = min(len(primer), len(read)).
    Requires check_len >= min_prefix and (matches / check_len) > threshold.
    """
    rlen = len(trimmed_read_lower)
    if rlen < min_prefix:
        return None
    for pid in candidate_indices:
        pseq = primers_lower[pid]
        check_len = len(pseq) if len(pseq) < rlen else rlen
        if check_len >= min_prefix:
            matches = sum(1 for a, b in zip(trimmed_read_lower[:check_len], pseq[:check_len]) if a == b)
            if matches / check_len > threshold:
                return pid
    return None


def is_dimer_insert(insert_len: int, len_p1: int, len_p2: int) -> bool:
    """SADDLE primer-dimer length rule: insert < L1 + L2 - 2 bp."""
    return insert_len < (len_p1 + len_p2 - 2)


def edlib_shw_identity(trimmed_read_lower: str, ref_amplicon_lower: str) -> float:
    """
    Compute Myers' bit-vector semi-global (prefix/SHW) edit-distance similarity:
        similarity = 1 - d_E / len(trimmed_read)
    """
    if edlib is None:  # pragma: no cover
        raise RuntimeError("edlib package is required for gapped sequence alignment.")
    if not trimmed_read_lower:
        return 0.0
    res = edlib.align(trimmed_read_lower, ref_amplicon_lower, mode="SHW", task="distance")
    dist = res["editDistance"]
    if dist < 0:
        return 0.0
    return 1.0 - (dist / len(trimmed_read_lower))


def ungapped_positional_identity(trimmed_read_lower: str, ref_amplicon_lower: str) -> float:
    """Compute ungapped 5'-anchored positional character match fraction."""
    if not trimmed_read_lower or not ref_amplicon_lower:
        return 0.0
    comp_len = min(len(trimmed_read_lower), len(ref_amplicon_lower))
    matches = sum(
        1 for a, b in zip(trimmed_read_lower[:comp_len], ref_amplicon_lower[:comp_len]) if a == b
    )
    return matches / len(trimmed_read_lower)


def classify_best_target_edlib(
    orient_read_lower: str,
    amplicons: List[str],
    amplicons_rc: List[str],
    assigned_amp_idx: int,
) -> dict:
    """
    Evaluate dual-strand semi-global (SHW) edit distance across all N amplicons,
    explicitly collapsing forward (+) and reverse-complement (-) strands per amplicon
    before comparing the assigned cognate amplicon against all non-cognate amplicons.

    This prevents a reverse-strand hit on the SAME amplicon from being miscounted as a
    non-cognate match and cleanly separates:
      - strictly_better_noncognate: min_{j != assigned} d_amp(j) < d_amp(assigned)
      - tied_best_noncognate:       min_{j != assigned} d_amp(j) == d_amp(assigned)
    """
    if edlib is None:  # pragma: no cover
        raise RuntimeError("edlib package is required for gapped sequence alignment.")
    rlen = len(orient_read_lower)
    plex = len(amplicons)
    if rlen == 0 or plex == 0:
        return {
            "assigned_fwd_sim": -1.0,
            "best_overall_sim": 0.0,
            "best_amp_1based": 0,
            "best_strand": ".",
            "strictly_better_noncognate": False,
            "tied_best_noncognate": False,
        }

    assigned_fwd_sim = -1.0
    assigned_best_dist = rlen + 1
    best_overall_dist = rlen + 1
    best_amp_idx = -1
    best_strand = "."
    best_noncog_dist = rlen + 1

    if assigned_amp_idx != -1:
        res_f = edlib.align(orient_read_lower, amplicons[assigned_amp_idx], mode="SHW", task="distance")
        d_f = res_f["editDistance"] if res_f["editDistance"] >= 0 else rlen
        assigned_fwd_sim = 1.0 - (d_f / rlen)
        res_r = edlib.align(orient_read_lower, amplicons_rc[assigned_amp_idx], mode="SHW", task="distance")
        d_r = res_r["editDistance"] if res_r["editDistance"] >= 0 else rlen
        if d_f <= d_r:
            assigned_best_dist = d_f
            best_overall_dist = d_f
            best_amp_idx = assigned_amp_idx
            best_strand = "+"
        else:
            assigned_best_dist = d_r
            best_overall_dist = d_r
            best_amp_idx = assigned_amp_idx
            best_strand = "-"

    for i in range(plex):
        if i == assigned_amp_idx:
            continue
        k_lim = best_overall_dist if best_overall_dist <= rlen else rlen
        rf = edlib.align(orient_read_lower, amplicons[i], mode="SHW", task="distance", k=k_lim)
        df = rf["editDistance"]
        rr = edlib.align(orient_read_lower, amplicons_rc[i], mode="SHW", task="distance", k=k_lim)
        dr = rr["editDistance"]

        if df == -1 and dr == -1:
            continue
        if df != -1 and (dr == -1 or df <= dr):
            d_amp, s_amp = df, "+"
        else:
            d_amp, s_amp = dr, "-"

        if assigned_amp_idx != -1 and d_amp < best_noncog_dist:
            best_noncog_dist = d_amp
        if d_amp < best_overall_dist:
            best_overall_dist = d_amp
            best_amp_idx = i
            best_strand = s_amp

    best_overall_sim = 1.0 - (best_overall_dist / rlen) if best_amp_idx != -1 else 0.0
    strictly_better = (assigned_amp_idx != -1) and (best_noncog_dist < assigned_best_dist)
    tied_best = (assigned_amp_idx != -1) and (best_noncog_dist == assigned_best_dist)
    return {
        "assigned_fwd_sim": assigned_fwd_sim,
        "best_overall_sim": best_overall_sim,
        "best_amp_1based": best_amp_idx + 1 if best_amp_idx != -1 else 0,
        "best_strand": best_strand,
        "strictly_better_noncognate": strictly_better,
        "tied_best_noncognate": tied_best,
    }
