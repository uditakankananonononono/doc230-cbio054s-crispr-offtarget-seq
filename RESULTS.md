# CBIO054-S results - partial: 3D core blocked (no coordinates in public encoded data). Sequence-only. No 3D/chromatin claim.
All numbers from results/results.json (src/run.py, results/run.log). CIRCLE-seq: 584,949 pairs, 10 guides, 7,371 positives (base rate 1.3%). Mean over 5 folds.
- S1 leave-guide-out: AUPRC B1 mismatch count 0.054, B2 position LR 0.097, CNN 0.282 (AUROC 0.74 / 0.88 / 0.93). Gate (CNN >= best baseline + 0.02) MET.
- S2 random row split: CNN AUPRC 0.606 vs 0.282 cold-guide, gap 0.32. Gate (<= 0.05) NOT met: random splits strongly overstate performance (near-duplicate off-target sites of the same guide leak across folds).
- S3 position ablation (B2 only): seed (PAM-proximal) mismatches AUPRC 0.029 AUROC 0.713; distal 0.037 / 0.661. Mixed result, no clean seed dominance claim.
- S4 transfer CIRCLE-seq -> GUIDE-seq: Kleinstiver (54 positives) AUPRC B1/B2/CNN 0.116/0.113/0.118; Listgarten (56 positives) 0.037/0.057/0.081 (AUROC CNN 0.992). Positives are few, differences within noise except possibly Listgarten; no robust CNN advantage claimed.
- S5 shuffle control (CNN, fold 0): AUPRC 0.011 vs base 0.007, AUROC 0.589 (slightly above chance, one fold only, noted). B2 negative ratio 1:1, 10:1, 50:1 AUPRC 0.095, 0.090, 0.084 (small effect).
Limits: negatives are in-silico candidate sites; CIRCLE-seq is in vitro; guides few (10); 3 CNN epochs; no 3D or epigenomic features, which the parent project's title implies - not tested here at all.
