#!/usr/bin/env python3
"""
code/render_manuscript.py

Deterministic UTF-8/LF template renderer for the ReScience C manuscript
(both LaTeX rescience.cls and synchronized root paper.md) and README.md.
Every quantitative value is populated strictly from data/canonical_results.json.
"""

import gzip
import json
import re
import sys
from pathlib import Path
from typing import Dict

sys.path.insert(0, str(Path(__file__).resolve().parent))
import canonical_classifier as cc


def flatten_metrics(canon: dict) -> Dict[str, str]:
    m: Dict[str, str] = {}

    env = canon.get("environment", {})
    m["env_python"] = str(env.get("python", "unknown"))
    m["env_bowtie2"] = str(env.get("bowtie2", "not run"))
    m["env_numpy"] = str(env.get("numpy", "unknown"))
    m["env_openpyxl"] = str(env.get("openpyxl", "unknown"))

    total_raw_all = 0
    not_opt_rs = []
    conc_opt_rs = []
    conc_opt_ratios = []
    conc_opt_l1s = []
    conc_opt_jaccards = []
    conc_opt_exact_all = []
    conc_opt_exact_active = []
    conc_opt_low12_cell_pcts = []
    conc_opt_low12_l1_pcts = []

    for pid in ("PS1", "PS2", "PS3", "E19"):
        if pid not in canon["panels"]:
            continue
        p = canon["panels"][pid]
        k = pid.lower()

        pa = p.get("pipeline_a_bowtie2", {})
        pb = p["pipeline_b_edlib"]
        pub = p["published_reference"]

        raw_n = pb["counts"]["total_raw_reads"]
        total_raw_all += raw_n
        m[f"{k}_raw"] = f"{raw_n:,}"
        m[f"{k}_plex"] = str(p["plex"])
        m[f"{k}_primers"] = str(p["n_primers"])

        if pub["dimer_pct"] is not None:
            m[f"{k}_pub_dimer"] = f"{pub['dimer_pct']:.1f}"
            m[f"{k}_pub_dimer_int"] = f"{pub['dimer_pct']:.0f}"
        if pub["ontarget_pct"] is not None:
            m[f"{k}_pub_ontarget"] = f"{pub['ontarget_pct']:.0f}"
        if pub["nonspecific_pct"] is not None:
            m[f"{k}_pub_nonspec"] = f"{pub['nonspecific_pct']:.0f}"

        if pa:
            pa_c = pa["counts"]
            pa_p = pa["pct_of_classified_ipynb"]
            pa_pboth = pa["pct_of_classified_both_mapped"]
            pa_raw = pa["pct_of_raw_reads"]
            pa_dd = pa["dimer_distribution"]

            m[f"{k}_pa_classified"] = f"{pa['classified_denominator_ipynb']:,}"
            m[f"{k}_pa_ontarget_n"] = f"{pa['ontarget_ipynb_formula_pairs']:,}"
            m[f"{k}_pa_ontarget_both_n"] = f"{pa['ontarget_both_mapped_pairs']:,}"
            m[f"{k}_pa_floor_extra"] = f"{pa['ontarget_floor_div_extra_pairs']:,}"
            m[f"{k}_pa_nonspec_n"] = f"{pa_c['nonspecific']:,}"
            m[f"{k}_pa_dimer_n"] = f"{pa_c['dimers']:,}"
            m[f"{k}_pa_ontarget_pct"] = f"{pa_p['ontarget_pct']:.2f}"
            m[f"{k}_pa_nonspec_pct"] = f"{pa_p['nonspecific_pct']:.2f}"
            m[f"{k}_pa_dimer_pct"] = f"{pa_p['dimer_pct']:.2f}"
            m[f"{k}_pa_both_ontarget_pct"] = f"{pa_pboth['ontarget_pct']:.2f}"
            m[f"{k}_pa_raw_ontarget_pct"] = f"{pa_raw['ontarget_ipynb_pct']:.2f}"
            m[f"{k}_pa_raw_nonspec_pct"] = f"{pa_raw['nonspecific_pct']:.2f}"
            m[f"{k}_pa_raw_dimer_pct"] = f"{pa_raw['dimer_pct']:.2f}"
            m[f"{k}_pa_raw_excl_pct"] = f"{pa_raw['excluded_pct']:.2f}"
            m[f"{k}_pa_mixed_dropped"] = f"{pa_c['mixed_length_60bp_dropped_matlab_rule']:,}"
            m[f"{k}_pa_one_mate_mapped"] = f"{pa_c['bt2_one_mate_mapped']:,}"
            m[f"{k}_pa_dimer_active"] = f"{pa_dd['active_dimer_pairs_ge1']:,}"
            m[f"{k}_pa_dimer_top1_pct"] = f"{pa_dd['top1_pct_of_dimer_reads']:.2f}"
            if pa_dd["top5_pairs"]:
                m[f"{k}_pa_top1_reads"] = f"{pa_dd['top5_pairs'][0]['dimer_reads']:,}"

            if pub["dimer_pct"] is not None:
                m[f"{k}_pa_dimer_diff"] = f"{pa_p['dimer_pct'] - pub['dimer_pct']:+.2f}"

            if "moesm7_matrix_validation" in pa:
                mv = pa["moesm7_matrix_validation"]
                co = mv["conc_optimized"]
                not_opt_rs.append(mv["not_optimized"]["pearson_r_unordered"])
                conc_opt_rs.append(co["pearson_r_unordered"])
                conc_opt_ratios.append(co["total_read_ratio_rec_over_pub"])
                conc_opt_l1s.append(co["normalized_l1_error"])
                conc_opt_jaccards.append(co["active_pair_jaccard"])
                conc_opt_exact_all.append(co["exact_cell_match_pct"])
                conc_opt_exact_active.append(co["exact_cell_match_active_union_pct"])
                conc_opt_low12_cell_pcts.append(co["disagreement_stratification"]["max_count_1_to_2"]["pct_of_disagreeing_cells"])
                conc_opt_low12_l1_pcts.append(co["disagreement_stratification"]["max_count_1_to_2"]["pct_of_total_l1"])

                for tag in ("not_optimized", "conc_optimized"):
                    t = mv[tag]
                    m[f"{k}_pa_mv_{tag}_r"] = f"{t['pearson_r_unordered']:.4f}"
                    m[f"{k}_pa_mv_{tag}_r5"] = f"{t['pearson_r_unordered']:.5f}"
                    m[f"{k}_pa_mv_{tag}_r6"] = f"{t['pearson_r_unordered']:.6f}"
                    m[f"{k}_pa_mv_{tag}_rho"] = f"{t['spearman_rho_unordered']:.4f}"
                    m[f"{k}_pa_mv_{tag}_pub_sum"] = f"{t['published_total_dimer_reads']:,}"
                    m[f"{k}_pa_mv_{tag}_rec_sum"] = f"{t['reconstructed_total_dimer_reads']:,}"
                    m[f"{k}_pa_mv_{tag}_ratio"] = f"{t['total_read_ratio_rec_over_pub']:.4f}"
                    m[f"{k}_pa_mv_{tag}_top10"] = str(t["top10_overlap_count"])
                    m[f"{k}_pa_mv_{tag}_exact_pct"] = f"{t['exact_cell_match_pct']:.2f}"
                    m[f"{k}_pa_mv_{tag}_exact_act_pct"] = f"{t['exact_cell_match_active_union_pct']:.1f}"
                    m[f"{k}_pa_mv_{tag}_jaccard"] = f"{t['active_pair_jaccard']:.4f}"
                    m[f"{k}_pa_mv_{tag}_mae_all"] = f"{t['mae_all_pairs']:.2f}"
                    m[f"{k}_pa_mv_{tag}_mae_union"] = f"{t['mae_active_union']:.2f}"
                    m[f"{k}_pa_mv_{tag}_norm_l1"] = f"{t['normalized_l1_error'] * 100.0:.2f}"
                    ds = t["disagreement_stratification"]
                    m[f"{k}_pa_mv_{tag}_dis_n"] = f"{t['disagreeing_cells_total']:,}"
                    m[f"{k}_pa_mv_{tag}_low12_pct"] = f"{ds['max_count_1_to_2']['pct_of_disagreeing_cells']:.1f}"
                    m[f"{k}_pa_mv_{tag}_low12_l1_pct"] = f"{ds['max_count_1_to_2']['pct_of_total_l1']:.1f}"
                    m[f"{k}_pa_mv_{tag}_ge10_l1_pct"] = f"{ds['max_count_ge_10']['pct_of_total_l1']:.1f}"

        pb_c = pb["counts"]
        pb_ex = pb["exclusions_summary"]
        pb_dd = pb["dimer_distribution"]

        m[f"{k}_pb_classified"] = f"{pb_c['classified_total']:,}"
        m[f"{k}_pb_excluded_pct_raw"] = f"{pb_ex['total_excluded_pct_raw']:.2f}"
        m[f"{k}_pb_unmatched_pct_raw"] = f"{pb_ex['primer_unmatched_pct_raw']:.2f}"
        m[f"{k}_pb_dimer_n"] = f"{pb_c['dimers']:,}"
        m[f"{k}_pb_dimer_pct"] = f"{pb['method1_primer_pairing']['dimer_pct_classified']:.2f}"
        m[f"{k}_pb_dimer_raw_pct"] = f"{pb['method1_primer_pairing']['dimer_pct_raw']:.2f}"

        for mkey, short in [
            ("method1_primer_pairing", "m1"),
            ("method2_edlib_75", "m2_75"),
            ("method2_edlib_80", "m2_80"),
            ("method2_edlib_85", "m2_85"),
            ("ungapped_positional_75", "ungap75"),
        ]:
            blk = pb[mkey]
            m[f"{k}_pb_{short}_on_n"] = f"{blk['ontarget_count']:,}"
            m[f"{k}_pb_{short}_on_pct"] = f"{blk['ontarget_pct_classified']:.2f}"
            m[f"{k}_pb_{short}_on_raw"] = f"{blk['ontarget_pct_raw']:.2f}"
            m[f"{k}_pb_{short}_non_pct"] = f"{blk['nonspecific_pct_classified']:.2f}"

        m[f"{k}_pb_total_pairs"] = f"{pb_dd['total_candidate_pairs']:,}"
        m[f"{k}_pb_zero_pairs"] = f"{pb_dd['zero_dimer_pairs']:,}"
        m[f"{k}_pb_zero_pct"] = f"{pb_dd['zero_dimer_pct']:.2f}"
        m[f"{k}_pb_active_pairs"] = f"{pb_dd['active_dimer_pairs_ge1']:,}"
        m[f"{k}_pb_active_pct"] = f"{pb_dd['active_dimer_pct_ge1']:.3f}"
        m[f"{k}_pb_severe_pairs"] = f"{pb_dd['severe_dimer_pairs_ge10']:,}"
        m[f"{k}_pb_severe_pct"] = f"{pb_dd['severe_dimer_pct_ge10']:.4f}"
        m[f"{k}_pb_top1_pct"] = f"{pb_dd['top1_pct_of_dimer_reads']:.2f}"
        m[f"{k}_pb_top10_pct"] = f"{pb_dd['top10_pct_of_dimer_reads']:.2f}"
        if pb_dd["top5_pairs"]:
            t1 = pb_dd["top5_pairs"][0]
            m[f"{k}_pb_top1_p1"] = t1["primer_1"]
            m[f"{k}_pb_top1_p2"] = t1["primer_2"]
            m[f"{k}_pb_top1_reads"] = f"{t1['dimer_reads']:,}"

        if "moesm7_matrix_validation" in pb:
            mv = pb["moesm7_matrix_validation"]
            for tag in ("not_optimized", "conc_optimized"):
                t = mv[tag]
                m[f"{k}_pb_mv_{tag}_r"] = f"{t['pearson_r_unordered']:.4f}"
                m[f"{k}_pb_mv_{tag}_r5"] = f"{t['pearson_r_unordered']:.5f}"

    m["all_panels_total_raw"] = f"{total_raw_all:,}"
    if not_opt_rs:
        m["pa_not_opt_r_min"] = f"{min(not_opt_rs):.4f}"
        m["pa_not_opt_r_max"] = f"{max(not_opt_rs):.4f}"
    if conc_opt_rs:
        m["pa_conc_opt_r6_min"] = f"{min(conc_opt_rs):.6f}"
        m["pa_conc_opt_r6_max"] = f"{max(conc_opt_rs):.6f}"
        m["pa_conc_opt_ratio_min"] = f"{min(conc_opt_ratios):.4f}"
        m["pa_conc_opt_ratio_max"] = f"{max(conc_opt_ratios):.4f}"
        m["pa_conc_opt_l1_min"] = f"{min(conc_opt_l1s) * 100.0:.2f}"
        m["pa_conc_opt_l1_max"] = f"{max(conc_opt_l1s) * 100.0:.2f}"
        m["pa_conc_opt_jacc_min"] = f"{min(conc_opt_jaccards):.4f}"
        m["pa_conc_opt_jacc_max"] = f"{max(conc_opt_jaccards):.4f}"
        m["pa_conc_opt_exact_all_min"] = f"{min(conc_opt_exact_all):.2f}"
        m["pa_conc_opt_exact_all_max"] = f"{max(conc_opt_exact_all):.2f}"
        m["pa_conc_opt_exact_act_min"] = f"{min(conc_opt_exact_active):.1f}"
        m["pa_conc_opt_exact_act_max"] = f"{max(conc_opt_exact_active):.1f}"
        m["pa_conc_opt_low12_cell_min"] = f"{min(conc_opt_low12_cell_pcts):.1f}"
        m["pa_conc_opt_low12_cell_max"] = f"{max(conc_opt_low12_cell_pcts):.1f}"
        m["pa_conc_opt_low12_l1_min"] = f"{min(conc_opt_low12_l1_pcts):.1f}"
        m["pa_conc_opt_low12_l1_max"] = f"{max(conc_opt_low12_l1_pcts):.1f}"

    ra = canon.get("random_audits_e19", {})
    if ra:
        m["audit_seed"] = str(ra["seed"])
        m["audit_pop_nondimer"] = f"{ra['population_nondimer_reads']:,}"
        m["audit_pop_cognate"] = f"{ra['population_cognate_nondimer_reads']:,}"
        for thr_key, tag in [("all_384_target_concordance_75", "c75"), ("all_384_target_concordance_80", "c80")]:
            c = ra[thr_key]
            m[f"{tag}_n"] = f"{c['sample_size']:,}"
            m[f"{tag}_cog_tot"] = f"{c['cognate_total']:,}"
            m[f"{tag}_both_on"] = f"{c['both_ontarget']:,}"
            m[f"{tag}_both_off"] = f"{c['both_offtarget']:,}"
            m[f"{tag}_cog_off_all_on"] = f"{c['a_off_b_on']:,}"
            m[f"{tag}_cog_on_all_off"] = f"{c['a_on_b_off']:,}"
            m[f"{tag}_cog_cross_best"] = f"{c['cognate_cross_target_best']:,}"
            m[f"{tag}_cog_tied_best"] = f"{c['cognate_tied_noncognate_best']:,}"
            m[f"{tag}_cog_on_pct"] = f"{c['method_a_ontarget_pct']:.2f}"
            m[f"{tag}_cog_on_ci"] = f"{c['method_a_ontarget_ci95'][0]:.2f}--{c['method_a_ontarget_ci95'][1]:.2f}"
            m[f"{tag}_all_on_pct"] = f"{c['method_b_ontarget_pct']:.2f}"
            m[f"{tag}_all_on_ci"] = f"{c['method_b_ontarget_ci95'][0]:.2f}--{c['method_b_ontarget_ci95'][1]:.2f}"
            m[f"{tag}_conc_pct"] = f"{c['concordance_pct']:.2f}"
            m[f"{tag}_conc_ci"] = f"{c['concordance_ci95'][0]:.2f}--{c['concordance_ci95'][1]:.2f}"
            m[f"{tag}_disc_pct"] = f"{c['discordance_pct']:.2f}"

        bm = ra.get("matched_pair_bowtie2_audit")
        if bm:
            m["bm_ge60_n"] = f"{bm['pairs_ge60_both_mates']:,}"
            m["bm_lt60_n"] = f"{bm['pairs_lt60_either_mate']:,}"
            m["bm_bt2_n"] = f"{bm['bt2_both_mapped_any_count']:,}"
            m["bm_bt2_conc_n"] = f"{bm['bt2_concordant_count']:,}"
            m["bm_bt2_pct"] = f"{bm['bt2_both_mapped_any_pct']:.2f}"
            m["bm_bt2_ci"] = f"{bm['bt2_both_mapped_any_ci95'][0]:.2f}--{bm['bt2_both_mapped_any_ci95'][1]:.2f}"
            m["bm_bt2_ge60_pct"] = f"{bm['bt2_ge60_both_mapped_pct']:.2f}"
            m["bm_bt2_ge60_ci"] = f"{bm['bt2_ge60_both_mapped_ci95'][0]:.2f}--{bm['bt2_ge60_both_mapped_ci95'][1]:.2f}"
            for key, prefix in [
                ("vs_cognate_edlib75", "bm_c75"),
                ("vs_cognate_edlib80", "bm_c80"),
                ("vs_all_target_edlib75", "bm_a75"),
                ("vs_all_target_edlib80", "bm_a80"),
            ]:
                cm = bm[key]
                m[f"{prefix}_conc"] = f"{cm['concordance_pct']:.2f}"
                m[f"{prefix}_conc_ci"] = f"{cm['concordance_ci95'][0]:.2f}--{cm['concordance_ci95'][1]:.2f}"
                m[f"{prefix}_diff"] = f"{cm['diff_b_minus_a_pp']:+.2f}"

        ind = ra["indel_recovery_audit"]
        m["ind_n"] = f"{ind['sample_size']:,}"
        m["ind_rec_n"] = f"{ind['recovered_by_gapped_count']:,}"
        m["ind_rec_pct"] = f"{ind['recovered_pct_of_cognate_sample']:.2f}"
        m["ind_rec_ci"] = f"{ind['recovered_pct_ci95'][0]:.2f}--{ind['recovered_pct_ci95'][1]:.2f}"
        m["ind_gaps_n"] = f"{ind['reads_with_gaps_count']:,}"
        m["ind_gaps_pct"] = f"{ind['reads_with_gaps_pct']:.1f}"
        m["ind_subonly_n"] = f"{ind['substitution_only_count']:,}"
        m["ind_le2_n"] = f"{ind['small_indels_le2bp_count']:,}"
        m["ind_le2_pct"] = f"{ind['small_indels_le2bp_pct']:.2f}"
        m["ind_gt2_n"] = f"{ind['larger_or_multi_indels_gt2bp_count']:,}"
        m["ind_gt2_pct"] = f"{ind['larger_or_multi_indels_gt2bp_pct']:.2f}"
        m["ind_med_bp"] = f"{ind['median_indel_bp']:.1f}"
        m["ind_iqr_lo"] = f"{ind['iqr_indel_bp'][0]:.1f}"
        m["ind_iqr_hi"] = f"{ind['iqr_indel_bp'][1]:.1f}"

    # Programmatically derive below-threshold non-cognate match summary from archived per-read audit TSV
    audit_tsv = Path(__file__).resolve().parent.parent / "data" / "matched_pair_audit_5000.tsv.gz"
    if audit_tsv.exists():
        below_cog = []
        below_alt = []
        with gzip.open(audit_tsv, "rt", encoding="utf-8") as f:
            _ = f.readline()
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 12 and parts[2] == "1" and parts[10] == "1" and float(parts[7]) < 0.75:
                    below_cog.append(float(parts[6]) * 100.0)
                    below_alt.append(float(parts[7]) * 100.0)
        m["audit_strict_below75_n"] = str(len(below_cog))
        if below_cog:
            m["audit_strict_below75_cog_range"] = f"{min(below_cog):.2f}--{max(below_cog):.2f}"
            m["audit_strict_below75_alt_range"] = f"{min(below_alt):.2f}--{max(below_alt):.2f}"
        else:
            m["audit_strict_below75_cog_range"] = "n/a"
            m["audit_strict_below75_alt_range"] = "n/a"

    return m


def render_template(src: Path, dst: Path, metrics: Dict[str, str]) -> None:
    text = src.read_text(encoding="utf-8")
    missing = set()

    def repl(match: re.Match) -> str:
        key = match.group(1)
        if key not in metrics:
            missing.add(key)
            return match.group(0)
        return metrics[key]

    rendered = re.sub(r"\{\{([a-zA-Z0-9_]+)\}\}", repl, text)
    if missing:
        raise KeyError(f"Unresolved template tokens in {src.name}: {sorted(missing)}")
    cc.write_utf8_lf(dst, rendered)


def main() -> None:
    base_dir = Path(__file__).resolve().parent.parent
    canon = json.loads((base_dir / "data" / "canonical_results.json").read_text(encoding="utf-8"))
    metrics = flatten_metrics(canon)

    pairs = [
        (base_dir / "article" / "content.tex.in", base_dir / "article" / "content.tex"),
        (base_dir / "article" / "metadata.yaml.in", base_dir / "article" / "metadata.yaml"),
        (base_dir / "paper.md.in", base_dir / "paper.md"),
        (base_dir / "README.md.in", base_dir / "README.md"),
    ]
    for src, dst in pairs:
        if src.exists():
            render_template(src, dst, metrics)
            print(f"Rendered {dst.relative_to(base_dir).as_posix()} from {src.name} (UTF-8/LF)")


if __name__ == "__main__":
    main()
