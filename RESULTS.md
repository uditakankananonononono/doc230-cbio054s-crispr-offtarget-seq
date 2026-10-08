# CBIO054-S results - partial: 3D core blocked (no coordinates in public encoded data). Sequence-only. No 3D/chromatin claim.
All numbers from results/results.json (src/run.py, results/run.log). CIRCLE-seq: 584,949 pairs, 10 guides, 7,371 positives (base rate 1.3%). Mean over 5 folds.
- S1 leave-guide-out: AUPRC B1 mismatch count 0.054, B2 position LR 0.097, CNN 0.282 (AUROC 0.74 / 0.88 / 0.93). Gate (CNN >= best baseline + 0.02) MET.
- S2 random row split: CNN AUPRC 0.606 vs 0.282 cold-guide, gap 0.32. Gate (<= 0.05) NOT met: random splits strongly overstate performance (near-duplicate off-target sites of the same guide leak across folds).
- S3 position ablation (B2 only): seed (PAM-proximal) mismatches AUPRC 0.029 AUROC 0.713; distal 0.037 / 0.661. Mixed result, no clean seed dominance claim.
- S4 transfer CIRCLE-seq -> GUIDE-seq: Kleinstiver (54 positives) AUPRC B1/B2/CNN 0.116/0.113/0.118; Listgarten (56 positives) 0.037/0.057/0.081 (AUROC CNN 0.992). Positives are few, differences within noise except possibly Listgarten; no robust CNN advantage claimed.
- S5 shuffle control (CNN, fold 0): AUPRC 0.011 vs base 0.007, AUROC 0.589 (slightly above chance, one fold only, noted). B2 negative ratio 1:1, 10:1, 50:1 AUPRC 0.095, 0.090, 0.084 (small effect).
Limits: negatives are in-silico candidate sites; CIRCLE-seq is in vitro; guides few (10); 3 CNN epochs; no 3D or epigenomic features, which the parent project's title implies - not tested here at all.

## Addendum: 5-fold shuffle control (requested after review; scripts src/shuffle5.py, src/shuffle5_seeds.py)
Labels of the training set permuted, CNN trained as in S1, evaluated on the same held-out guides.
- Seed set 0: per-fold AUROC 0.589, 0.408, 0.442, 0.533, 0.694 (mean 0.533); AUPRC 0.011, 0.010, 0.012, 0.010, 0.043 (base 0.007-0.023).
- Seed set 1: AUROC 0.547, 0.422, 0.638, 0.438, 0.689 (mean 0.547); AUPRC 0.008-0.039.
- A third seed set was killed by memory pressure (two jobs in parallel) and is not reported.
Reading: the control does NOT return cleanly to 0.5. Fold-to-fold spread is large (0.41-0.69) and fold 4 gives about 0.69 in both seed sets, so part of the deviation is tied to that fold's test guides, not training noise. The cause is not established. For scale: real S1 AUPRC (0.28 mean) is far above every shuffled AUPRC (max 0.043), but the AUROC control means AUROC differences of ~0.1 on single folds should not be read as signal. No refit was done after seeing this.

## Audit (second pass): fold-4 shuffle anomaly - NO LEAKAGE FOUND, residual unexplained
Prereg: AUDIT_PREREG.md (committed before the audit ran; A5 addendum committed after A1-A4 and before A5). Scripts src/audit.py, src/audit_a5.py; outputs results/audit_fold4.json, results/audit_a5.json. Re-run on a re-fetched copy of the data (SHA256 identical) with a newer torch than the original run (version differs; seeds are not bit-comparable).
- A1 guide overlap: min Hamming distance train-vs-test guides 9-12 in all folds (fold 4: 12); zero identical guides.
- A2 site overlap: 0% of test-positive (guide, site) pairs occur in train in any fold; off-site sequence alone appears in train positives for 0%.
- A3 untrained CNN, fold 4: AUROC 0.499 (5 seeds). So the elevated fold-4 shuffle AUROC is not an architecture prior. Mismatch-count AUROC is lowest in fold 4 (0.596), so it is not driven by an easy fold either.
- A4 within-guide label shuffle, fold 4: AUROC 0.618 / 0.503 / 0.564 (3 seeds); fold 0: 0.540 / 0.472 / 0.480.
- A5 global label shuffle (original style), 8 seeds: fold 4 mean AUROC 0.608 (range 0.46-0.70, sd ~0.08); fold 0 mean 0.540 (0.47-0.62). Shuffled AUPRC fold 4 0.019-0.046 (base 0.023).
Decision rules (pre-registered): leakage NOT shown; architecture explanation NOT supported; A5 rule gives "persistent fold-specific effect, unexplained" (mean 0.608 > 0.58). Honest reading: no train/test overlap found by any test I could define, shuffle AUROC on fold 4 is mildly and noisily above 0.5 (2 test guides, large seed variance), cause unexplained. The headline S1 gate is on AUPRC: real cold-guide AUPRC per fold 0.20/0.35/0.26/0.18/0.43 vs the highest shuffled AUPRC 0.046 (fold 4), a >4x gap in every fold, so the verdict is unaffected. AUROC differences of ~0.1 on single folds still should not be read as signal.
