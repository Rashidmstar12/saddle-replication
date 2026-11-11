#!/usr/bin/env python3
"""
code/run_baseline_benchmark.py

PHASE 0D: Controlled Baseline Benchmark
Trains baseline models on SADDLE 96-plex (PS1-PS3: 55,584 pairs) and evaluates
on the independent held-out SADDLE 384-plex (295,296 pairs).

Baselines:
1. SADDLE Badness heuristic
2. Conventional nearest-neighbor thermodynamics (SantaLucia delta G)
3. Logistic Regression on physical sequence features
4. XGBoost on the exact same physical sequence features

Evaluations:
- PR-AUC, ROC-AUC, Precision@Top-10, Top-20, Top-50, Top-100, Brier score
- Stratified sensitivity without dominant pair rp213 x rp371
- High-confidence dimer threshold (>= 10 reads)
"""

import sys
import json
import time
from pathlib import Path
import numpy as np
import openpyxl
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    average_precision_score,
    brier_score_loss,
    roc_auc_score,
)
from sklearn.preprocessing import StandardScaler
from xgboost import XGBClassifier

ROOT = Path(__file__).resolve().parent.parent
RAW_DIR = ROOT / "data" / "raw" if (ROOT / "data" / "raw").exists() else ROOT / "raw"
DERIVED_DIR = ROOT / "data" / "derived" if (ROOT / "data" / "derived").exists() else ROOT / "data"

MOESM5_PATH = RAW_DIR / "41467_2022_29500_MOESM5_ESM.xlsx"
MOESM6_PATH = RAW_DIR / "41467_2022_29500_MOESM6_ESM.xlsx"
MOESM7_PATH = RAW_DIR / "41467_2022_29500_MOESM7_ESM.xlsx"
UNORDERED_384_PATH = DERIVED_DIR / "saddle_384_matrix_unordered.npy"

# SantaLucia unified parameters (delta G at 37C in kcal/mol)
SANTALUCIA_DG37 = {
    'AA': -1.00, 'TT': -1.00,
    'AT': -0.88,
    'TA': -0.58,
    'CA': -1.45, 'TG': -1.45,
    'GT': -1.44, 'AC': -1.44,
    'CT': -1.28, 'AG': -1.28,
    'GA': -1.30, 'TC': -1.30,
    'CG': -2.17,
    'GC': -2.24,
    'GG': -1.84, 'CC': -1.84,
}

def revcomp(s: str) -> str:
    rc = str.maketrans("ACGTacgt", "TGCAtgca")
    return s.translate(rc)[::-1]

def compute_pairwise_physics(s1: str, s2: str):
    """
    Computes thermodynamic and physical sequence features between two primers:
    1. min_dg: Minimum hybridization free energy (kcal/mol, negative = stable)
    2. dg_3prime: Minimum free energy at 3' terminal overlap (last 5 bp)
    3. max_match: Longest continuous Watson-Crick match length
    4. max_3prime_match: Longest continuous match ending at 3' terminus
    5. gc_mean: Average GC content
    6. len_diff: Difference in length
    7. min_len: Shorter primer length
    """
    s1 = s1.upper()
    s2 = s2.upper()
    l1, l2 = len(s1), len(s2)
    rc2 = revcomp(s2)

    max_match = 0
    max_3prime_match = 0
    min_dg = 0.0
    min_3prime_dg = 0.0

    gc1 = sum(1 for c in s1 if c in "GC") / l1
    gc2 = sum(1 for c in s2 if c in "GC") / l2
    gc_mean = (gc1 + gc2) / 2.0
    len_diff = abs(l1 - l2)
    min_len = min(l1, l2)

    # Shift rc2 across s1: offset from -(l2 - 1) to (l1 - 1)
    for offset in range(-l2 + 1, l1):
        start_1 = max(0, offset)
        end_1 = min(l1, offset + l2)
        start_2 = start_1 - offset

        sub1 = s1[start_1:end_1]
        sub2 = rc2[start_2:start_2 + len(sub1)]

        cur_run = 0
        cur_dg = 0.0
        in_run = False

        for idx in range(len(sub1)):
            if sub1[idx] == sub2[idx]:
                cur_run += 1
                if cur_run > max_match:
                    max_match = cur_run
                if in_run:
                    dinuc = sub1[idx-1:idx+1]
                    cur_dg += SANTALUCIA_DG37.get(dinuc, -1.0)
                else:
                    in_run = True
                    cur_dg = 1.96  # Initiation penalty
            else:
                if in_run and cur_dg < min_dg:
                    min_dg = cur_dg
                cur_run = 0
                cur_dg = 0.0
                in_run = False

        if in_run and cur_dg < min_dg:
            min_dg = cur_dg

        # 3' end check for s1 (end_1 == l1)
        if end_1 == l1 and len(sub1) > 0 and sub1[-1] == sub2[-1]:
            run_3p = 0
            dg_3p = 1.96
            for idx in range(len(sub1)-1, -1, -1):
                if sub1[idx] == sub2[idx]:
                    run_3p += 1
                    if run_3p > 1:
                        dinuc = sub1[idx:idx+2]
                        dg_3p += SANTALUCIA_DG37.get(dinuc, -1.0)
                else:
                    break
            if run_3p > max_3prime_match:
                max_3prime_match = run_3p
            if dg_3p < min_3prime_dg:
                min_3prime_dg = dg_3p

        # 3' end check for s2 (start_2 == 0 since rc2 5' corresponds to s2 3')
        if start_2 == 0 and len(sub1) > 0 and sub1[0] == sub2[0]:
            run_3p = 0
            dg_3p = 1.96
            for idx in range(len(sub1)):
                if sub1[idx] == sub2[idx]:
                    run_3p += 1
                    if run_3p > 1:
                        dinuc = sub1[idx-1:idx+1]
                        dg_3p += SANTALUCIA_DG37.get(dinuc, -1.0)
                else:
                    break
            if run_3p > max_3prime_match:
                max_3prime_match = run_3p
            if dg_3p < min_3prime_dg:
                min_3prime_dg = dg_3p

    return min_dg, min_3prime_dg, max_match, max_3prime_match, gc_mean, len_diff, min_len

def compute_pairwise_badness(s1: str, s2: str) -> float:
    """
    Computes SADDLE pairwise interaction heuristic score.
    
    METHODOLOGICAL QUALIFICATION:
    This is an approximate 4-to-7-mer heuristic variant of the SADDLE Badness formula.
    The published SADDLE formula (Xie et al., Nat Commun 2022) describes complementary match
    lengths up to 8 bp. In this benchmark code, MID8 (256.0) was declared but internal 8-mer
    matches were not iterated over. Furthermore, in the baseline models (Logistic Regression
    and XGBoost), this approximate Badness score was included as an input feature alongside
    physical parameters. The benchmark therefore evaluates ML models augmenting this approximate
    heuristic, rather than a de novo replacement of the exact published SADDLE implementation.
    """
    p1 = s1.lower()
    p2 = s2.lower()
    score = 0.0
    END4, END5, END6 = 16.0, 32.0, 64.0
    MID7, MID8 = 128.0, 256.0

    for a, b in [(p1, p2), (p2, p1)]:
        rc_b = revcomp(b)
        l_rc = len(rc_b)
        l_a = len(a)

        k4 = a[-4:]
        k5 = a[-5:]
        k6 = a[-6:]

        for j in range(l_rc - 3):
            if rc_b[j:j+4] == k4:
                ngc = sum(1 for c in k4 if c in 'gc')
                score += (END4 / (j + 1.0)) * (2.0 ** ngc)
        for j in range(l_rc - 4):
            if rc_b[j:j+5] == k5:
                ngc = sum(1 for c in k5 if c in 'gc')
                score += (END5 / (j + 1.0)) * (2.0 ** ngc)
        for j in range(l_rc - 5):
            if rc_b[j:j+6] == k6:
                ngc = sum(1 for c in k6 if c in 'gc')
                score += (END6 / (j + 1.0)) * (2.0 ** ngc)

        for i_a in range(l_a - 6):
            m7 = a[i_a:i_a+7]
            dist_a = l_a - (i_a + 7)
            for j in range(l_rc - 6):
                if rc_b[j:j+7] == m7:
                    ngc = sum(1 for c in m7 if c in 'gc')
                    score += (MID7 / ((dist_a + 1.0) * (j + 1.0))) * (2.0 ** ngc)

    return score

def load_development_data():
    """Loads and computes features for SADDLE 96-plex (PS1-PS3: 55,584 pairs)."""
    print("Loading 96-plex development data (PS1, PS2, PS3)...")
    wb5 = openpyxl.load_workbook(MOESM5_PATH, data_only=True)
    wb6 = openpyxl.load_workbook(MOESM6_PATH, data_only=True)
    wb7 = openpyxl.load_workbook(MOESM7_PATH, data_only=True)

    panels = [
        ("PS1", "96plex PS1", "Badness_PS1", "PS1_Dimer NGS(not optimized)"),
        ("PS2", "96plex PS2", "Badness_PS2", "PS2_Dimer NGS(not optimized)"),
        ("PS3", "96plex PS3", "Badness_PS3", "PS3_Dimer NGS(not optimized)")
    ]

    X_list = []
    y_list = []
    triu_i, triu_j = np.triu_indices(192)

    for name, s5, s6, s7 in panels:
        ws5 = wb5[s5]
        fps = [r[1].upper() for r in ws5.iter_rows(values_only=True)]
        rps = [r[3].upper() for r in ws5.iter_rows(values_only=True)]
        primers = fps + rps
        assert len(primers) == 192

        ws6 = wb6[s6]
        bad_mat = np.array([[float(c or 0) for c in r] for r in ws6.iter_rows(values_only=True)])

        ws7 = wb7[s7]
        dim_mat = np.array([[float(c or 0) for c in r] for r in ws7.iter_rows(values_only=True)])

        dim_unordered = dim_mat[triu_i, triu_j] + np.where(triu_i != triu_j, dim_mat[triu_j, triu_i], 0)
        bad_unordered = bad_mat[triu_i, triu_j] + np.where(triu_i != triu_j, bad_mat[triu_j, triu_i], 0)

        y = (dim_unordered > 0).astype(np.int32)

        # Compute physics features for each pair
        print(f"  Extracting physics features for {name} ({len(dim_unordered)} pairs)...")
        panel_feats = []
        for idx in range(len(dim_unordered)):
            i, j = triu_i[idx], triu_j[idx]
            p1, p2 = primers[i], primers[j]
            sb = bad_unordered[idx]
            min_dg, dg_3p, max_m, max_3p, gc, ldiff, minl = compute_pairwise_physics(p1, p2)
            panel_feats.append([sb, min_dg, dg_3p, max_m, max_3p, gc, ldiff, minl])

        X_list.append(np.array(panel_feats, dtype=np.float32))
        y_list.append(y)

    X_train = np.vstack(X_list)
    y_train = np.concatenate(y_list)
    print(f"Development Set: {X_train.shape[0]:,} pairs, {np.sum(y_train):,} positives ({np.mean(y_train)*100:.2f}%).")
    return X_train, y_train

def load_test_data():
    """Loads and computes features for held-out SADDLE 384-plex (295,296 pairs)."""
    print("Loading 384-plex held-out test data (E19)...")
    wb5 = openpyxl.load_workbook(MOESM5_PATH, data_only=True)
    ws = wb5["384 Plex Primer"]
    fps, rps = [], []
    for r in ws.iter_rows(values_only=True):
        fps.append((r[0], str(r[1]).upper()))
        rps.append((r[2], str(r[3]).upper()))
    all_primers = fps + rps
    assert len(all_primers) == 768

    dim_mat_unordered = np.load(UNORDERED_384_PATH)
    triu_i, triu_j = np.triu_indices(768)
    reads_unordered = dim_mat_unordered[triu_i, triu_j]
    y_test = (reads_unordered > 0).astype(np.int32)
    y_test_ge10 = (reads_unordered >= 10).astype(np.int32)

    # Dominant pair check: rp213 (index 212+384=596) x rp371 (index 370+384=754)
    dominant_mask = np.ones(len(reads_unordered), dtype=bool)
    dominant_idx = None
    for idx in range(len(reads_unordered)):
        i, j = triu_i[idx], triu_j[idx]
        if (i == 596 and j == 754) or (i == 754 and j == 596):
            dominant_idx = idx
            dominant_mask[idx] = False
            break

    print(f"  Identified dominant pair at index {dominant_idx} with {reads_unordered[dominant_idx]} reads.")
    print(f"  Extracting features for 384-plex ({len(reads_unordered):,} pairs)...")
    t0 = time.time()
    feats = []
    for idx in range(len(reads_unordered)):
        i, j = triu_i[idx], triu_j[idx]
        p1 = all_primers[i][1]
        p2 = all_primers[j][1]
        sb = compute_pairwise_badness(p1, p2)
        min_dg, dg_3p, max_m, max_3p, gc, ldiff, minl = compute_pairwise_physics(p1, p2)
        feats.append([sb, min_dg, dg_3p, max_m, max_3p, gc, ldiff, minl])
        if (idx + 1) % 50000 == 0:
            print(f"    Computed {idx+1:,} / {len(reads_unordered):,} pairs ({time.time()-t0:.1f}s)...")

    X_test = np.array(feats, dtype=np.float32)
    print(f"Test Set: {X_test.shape[0]:,} pairs, {np.sum(y_test):,} positives ({np.mean(y_test)*100:.3f}%).")
    return X_test, y_test, y_test_ge10, dominant_mask, reads_unordered

def evaluate_predictions(scores, y_true, name="Model"):
    """Evaluates ranking, precision-recall, and top-k precision using deterministic tie-breaking."""
    roc = roc_auc_score(y_true, scores)
    pr = average_precision_score(y_true, scores)

    # Deterministic stable sort: -scores descending, index ascending
    sorted_order = np.lexsort((np.arange(len(scores)), -scores))
    sorted_y = y_true[sorted_order]

    p10 = np.mean(sorted_y[:10])
    p20 = np.mean(sorted_y[:20])
    p50 = np.mean(sorted_y[:50])
    p100 = np.mean(sorted_y[:100])

    return {
        "model": name,
        "roc_auc": round(float(roc), 4),
        "pr_auc": round(float(pr), 4),
        "precision_top10": round(float(p10), 4),
        "precision_top20": round(float(p20), 4),
        "precision_top50": round(float(p50), 4),
        "precision_top100": round(float(p100), 4),
    }

def run_benchmark():
    print("=" * 70)
    print("PHASE 0D: CONTROLLED BASELINE BENCHMARK PILOT")
    print("=" * 70)

    DERIVED_DIR.mkdir(parents=True, exist_ok=True)
    required = [MOESM5_PATH, MOESM6_PATH, MOESM7_PATH, UNORDERED_384_PATH]
    missing = [p.name for p in required if not p.exists()]
    if missing:
        print(f"SKIPPED_MISSING_INPUT: External raw or generated inputs missing: {missing}")
        print("This baseline modeling script requires external supplementary spreadsheets and raw matrix.")
        print("To verify benchmark metrics from bundled artifacts without retraining, run:")
        print("    python code/verify_replication.py")
        print("For raw data download instructions, see README.md.")
        sys.exit(2)

    X_train, y_train = load_development_data()
    X_test, y_test, y_test_ge10, dom_mask, reads = load_test_data()

    # Feature 0: SADDLE Badness
    # Feature 1: delta_g_min (lower = stronger, so negative delta_g for ranking)
    scores_badness = X_test[:, 0]
    scores_thermo = -X_test[:, 1]  # Invert so higher score = stronger binding

    # Model 3: Logistic Regression
    print("\nTraining Model 3: Logistic Regression...")
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)

    lr = LogisticRegression(class_weight="balanced", C=1.0, max_iter=1000, random_state=42)
    lr.fit(X_train_scaled, y_train)
    scores_lr = lr.predict_proba(X_test_scaled)[:, 1]
    brier_lr = brier_score_loss(y_test, scores_lr)

    # Model 4: XGBoost
    print("Training Model 4: XGBoost Classifier...")
    # Scale positive weight to handle 98% negative class
    scale_pos = (len(y_train) - np.sum(y_train)) / np.sum(y_train)
    xgb = XGBClassifier(
        n_estimators=100,
        max_depth=3,
        learning_rate=0.05,
        scale_pos_weight=scale_pos,
        random_state=42,
        eval_metric="logloss",
        n_jobs=-1
    )
    xgb.fit(X_train, y_train)
    scores_xgb = xgb.predict_proba(X_test)[:, 1]
    brier_xgb = brier_score_loss(y_test, scores_xgb)

    models = [
        ("SADDLE Badness (Heuristic)", scores_badness),
        ("SantaLucia Thermodynamics (-dG)", scores_thermo),
        ("Logistic Regression (Physics)", scores_lr),
        ("XGBoost (Physics Features)", scores_xgb),
    ]

    print("\n" + "=" * 70)
    print("EVALUATION 1: FULL HELD-OUT 384-PLEX (N = 295,296 pairs)")
    print(f"Prevalence: {np.mean(y_test)*100:.3f}% (1,172 / 295,296)")
    print("-" * 70)
    full_results = []
    for name, s in models:
        res = evaluate_predictions(s, y_test, name=name)
        if "Logistic" in name:
            res["brier_score"] = round(float(brier_lr), 4)
        elif "XGBoost" in name:
            res["brier_score"] = round(float(brier_xgb), 4)
        full_results.append(res)
        print(f"{name:<35}: PR-AUC={res['pr_auc']:.4f} | ROC-AUC={res['roc_auc']:.4f} | P@10={res['precision_top10']:.2f} | P@100={res['precision_top100']:.2f}")

    print("\n" + "=" * 70)
    print("EVALUATION 2: SENSITIVITY WITHOUT DOMINANT PAIR (N = 295,295 pairs)")
    print(f"Prevalence: {np.mean(y_test[dom_mask])*100:.3f}% (1,171 / 295,295)")
    print("-" * 70)
    sens_results = []
    for name, s in models:
        res = evaluate_predictions(s[dom_mask], y_test[dom_mask], name=name)
        sens_results.append(res)
        print(f"{name:<35}: PR-AUC={res['pr_auc']:.4f} | ROC-AUC={res['roc_auc']:.4f} | P@10={res['precision_top10']:.2f} | P@100={res['precision_top100']:.2f}")

    print("\n" + "=" * 70)
    print("EVALUATION 3: HIGH-CONFIDENCE DIMERS (Reads >= 10, N_pos = 98)")
    print(f"Prevalence: {np.mean(y_test_ge10)*100:.4f}% (98 / 295,296)")
    print("-" * 70)
    high_results = []
    for name, s in models:
        res = evaluate_predictions(s, y_test_ge10, name=name)
        high_results.append(res)
        print(f"{name:<35}: PR-AUC={res['pr_auc']:.4f} | ROC-AUC={res['roc_auc']:.4f} | P@10={res['precision_top10']:.2f} | P@100={res['precision_top100']:.2f}")

    benchmark_summary = {
        "dataset": "SADDLE 384-plex Held-Out Independent Benchmark",
        "development_dataset": "SADDLE 96-plex PS1-PS3 (55,584 pairs)",
        "test_dataset": "SADDLE 384-plex E19 (295,296 pairs)",
        "full_evaluation": full_results,
        "sensitivity_evaluation_no_dominant": sens_results,
        "high_confidence_evaluation_ge10": high_results,
        "logistic_regression_weights": {
            "features": ["saddle_badness", "min_dg", "dg_3prime", "max_match", "max_3prime_match", "gc_mean", "len_diff", "min_len"],
            "coefficients": [round(float(c), 4) for c in lr.coef_[0]],
            "intercept": round(float(lr.intercept_[0]), 4),
        },
        "xgboost_feature_importances": {
            "features": ["saddle_badness", "min_dg", "dg_3prime", "max_match", "max_3prime_match", "gc_mean", "len_diff", "min_len"],
            "importance": [round(float(imp), 4) for imp in xgb.feature_importances_]
        }
    }

    out_path = DERIVED_DIR / "baseline_benchmark_results.json"
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(benchmark_summary, f, indent=2)
    print(f"\nBenchmark results written to {out_path}")

if __name__ == "__main__":
    run_benchmark()
