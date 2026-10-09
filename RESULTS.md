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

## Extension batch E6-E10 (2026-10-09; PREREG_EXT.md committed first, 54d7d3f; src/ext_run.py; results/ext_results.json, ext_run.log)
Same three public CSVs (SHA256 re-verified identical), same CNN, B2 position-LR and the same 5 cold-guide folds (seed 1) as S1. All metrics here are AUPRC pooled over all cold-guide out-of-fold rows (S1 reported the mean of per-fold AUPRC), so absolute numbers are NOT directly comparable to S1; all comparisons below are within the same pool. One model seed per direction, no repeats: differences under about 0.02 are noise.
| Dir | Result | Gate |
|---|---|---|
| E6 mismatch strata | <=3 mismatches (n 4,337, 757 pos): CNN 0.561 vs B2 0.189; >3 mismatches (n 580,612, 6,614 pos): CNN 0.220 vs B2 0.068 | MET (CNN >= B2 + 0.02 in <=3 stratum) |
| E7 per-guide (10 CIRCLE-seq guides) | CNN beats B2 on 10 of 10 guides | MET (>= 7) |
| E8 training-size curve | 25% of guides 0.026, 50% 0.109, 100% 0.237 | MET (100% >= 25% + 0.02) |
| E9 calibration | ECE 0.027; mean predicted probability 0.039 vs true positive rate 0.0126 | informational; mean prediction is about 3x the true rate (observed fact). The cause is not established; 10:1 downsampling would imply a prior near 0.09, and 0.039 is below that, so that explanation does not fit. ECE 0.027 on a 1.3% base rate is weak calibration evidence; the relative overshoot is the better number |
| E10 guide-shuffled training | 0.082 vs original 0.237, drop 0.155 | MET (drop >= 0.05) |
Reading: within the sequence-only data, the CNN beats the mismatch-position LR in every stratum tested (this shows base identity or interaction matters beyond mismatch position, not where; the <=3-mismatch stratum has 17% positives vs 1.1% elsewhere, so AUPRC is not comparable across strata), including the low-mismatch stratum where the baseline has the most to work with, and the gain is on every guide. Alignment between guide and candidate matters (E10: the control permutes the guide array among training rows while the candidate stays, which destroys the guide-candidate alignment, i.e. the mismatch pattern; it cannot separate guide context from mismatch information) and from more training guides (E8: 4x guides roughly 9x AUPRC, which also says the 10-guide training set is small and the model is data-limited). This is a model-vs-baseline comparison on one fixed cold-guide split, with no cross-split claim. Limits: single seed, only 10 guides, CIRCLE-seq positives are in-vitro cleavage sites, no 3D or chromatin information exists in these data, and the 3D-aware core of the parent remains unbuilt. With this batch 054-S has 10 directions run (S1-S5 + E6-E10); the 3D core is still blocked.

### Gate review corrections (2026-10-09, wording only; no numbers changed)
- Pooled AUPRC pools out-of-fold scores from 5 separately trained fold models on different guides; this is another reason it is not comparable to S1 (mean of per-fold AUPRC).
- E8 subsamples whole GUIDES (a guide-count curve: 25/50/100% of training guides), not "rows" as PREREG_EXT.md says. Read E8 as a guide-count curve.
- E7 as prereg'd skips guides without positives; all 10 of 10 guides had positives here, so none were skipped.
- PREREG_EXT.md's "~600-700 CIRCLE-seq positives across 10 guides" does not match the data: 7,371 positives total (about 737 per guide). The prereg figure was a misstatement; it may have meant per guide, which I did not verify.
- Implementation (src/ext_run.py) was committed together with the results, so the prereg froze gates but not implementation choices.
- Accepted reading: E6-E10 ran as preregistered gates and all were met on one seed; E10 is an alignment-destroying control, not a guide-context test; E8 is a guide-count curve; pooled AUPRC is across five fold models.
