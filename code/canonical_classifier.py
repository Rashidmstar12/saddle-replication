#!/usr/bin/env python3
"""
code/canonical_classifier.py

Canonical Primer Matching and Classification Functions for:
"[Re] Computational Replication of SADDLE: Empirical Next-Generation Sequencing
      Read Partitions and Primer Interaction Distribution in a 384-Plex PCR Panel"

Author: Rashid Kadayil (ORCID: 0009-0009-6398-4557)

This module provides the single authoritative implementation of primer prefix
matching used across all pipeline scripts in this replication package.
"""

from typing import List, Tuple, Optional


def match_primer_prefix(
    trimmed_read: str,
    candidate_indices: List[int],
    all_primers: List[Tuple[str, str]],
    min_prefix: int = 8,
    threshold: float = 0.6,
) -> Optional[int]:
    """
    Canonical Primer Prefix Classifier (Primary Pipeline).

    Matches candidate primers identified by 8-bp 5' prefix index against
    the trimmed read sequence.

    Formula:
        check_len = min(len(primer_seq), len(trimmed_read))
        score = matching_bases / check_len

    Decision criteria:
        Requires check_len >= min_prefix (default: 8 bp).
        Returns the first primer index with score > threshold (default: 0.6).
        Returns None if no candidate satisfies the threshold.

    Rationale:
        Dividing by check_len correctly evaluates sequence identity over
        the observed extent of extreme-short reads (e.g. 12-14 bp primer-dimers)
        where the insert sequence is shorter than full primer length (~22 bp).
    """
    if len(trimmed_read) < min_prefix:
        return None

    for pid in candidate_indices:
        pseq = all_primers[pid][1]
        check_len = min(len(pseq), len(trimmed_read))
        if check_len >= min_prefix:
            matches = sum(1 for a, b in zip(trimmed_read[:check_len], pseq[:check_len]) if a == b)
            if matches / check_len > threshold:
                return pid
    return None


def match_primer_prefix_historical(
    trimmed_read: str,
    candidate_indices: List[int],
    all_primers: List[Tuple[str, str]],
    min_prefix: int = 8,
    threshold: float = 0.6,
) -> Optional[int]:
    """
    Historical Baseline Primer Prefix Classifier (Phase 0D Baseline).

    Matches candidate primers by 5' prefix, but calculates identity by
    dividing matching bases by FULL primer length:
        score = matching_bases / len(primer_seq)

    Decision criteria:
        Returns first primer index with score > threshold (default: 0.6).
        Returns None if no candidate satisfies the threshold.

    Historical Context:
        This function reproduces the earlier Phase 0D partition of 2,860,658
        classified reads and 24,872 primer-dimers. For trimmed reads of 12-14 bp,
        dividing by full primer length (~22 bp) yields score <= 0.60 even for
        perfect matches, excluding exactly 6 reads as primer-unmatched.
    """
    if len(trimmed_read) < min_prefix:
        return None

    for pid in candidate_indices:
        pseq = all_primers[pid][1]
        check_len = min(len(pseq), len(trimmed_read))
        if check_len >= min_prefix:
            matches = sum(1 for a, b in zip(trimmed_read[:check_len], pseq[:check_len]) if a == b)
            if matches / len(pseq) > threshold:
                return pid
    return None
