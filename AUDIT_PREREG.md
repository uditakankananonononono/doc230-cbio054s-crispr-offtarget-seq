# AUDIT PREREG - 054-S fold-4 shuffle-control anomaly (written before the audit ran)

Question: in the 5-fold shuffle control (CNN trained on permuted labels, evaluated on held-out guides), fold 4 gave AUROC ~0.69 in both seed sets (others 0.41-0.64). Is this leakage between train and test, or something else? Status was PROVISIONAL; it flips to NOT COUNTED if leakage is shown.

Same data (hashes in data_manifest/sha256.txt re-verified), same fold definition (guides sorted, rng(1) shuffled, folds[i::5]; "fold 4" = index 4).

## Tests
A1 Guide overlap: for each fold, the minimum Hamming distance (sgRNA_seq) between any test guide and any train guide; any identical guide sequences across groups.
A2 Site overlap: for each fold, the fraction of test POSITIVE off-target sites (sgRNA_seq+off_seq pair) and of test positive off_seq alone that also occur in the train set (as any row, and as positive).
A3 Untrained-network control: random-init CNN (same architecture, no training), 5 seeds, AUROC/AUPRC per fold; and B1 (negative mismatch count) AUROC per fold.
A4 Within-guide shuffle: training labels permuted within each training guide (preserves per-guide base rates), 3 seeds, fold 4 only plus fold 0 for comparison.

## Decision rules (fixed now)
LEAKAGE if either: (i) fold 4 has an identical guide sequence across train/test, or min guide Hamming distance < 3 where other folds have none; or (ii) fold 4 fraction of test positive pairs present in train is > 5% and at least 2x every other fold. Then 054-S flips to NOT COUNTED.
ARCHITECTURE/PRIOR EXPLANATION if leakage rules are not met AND the untrained network gives fold-4 AUROC >= 0.60 and its per-fold AUROC ranks fold 4 at or near the top (correlation across folds with B1 AUROC > 0.5). Then the shuffle control is mis-specified (a feature-driven prior, not a label signal), 054-S remains counted with a corrected control description.
UNRESOLVED otherwise: stays provisional, say so.
No re-fit of the original S1-S5 gates.

## Addendum A5 (written after A1-A4 results were seen, before A5 ran)
A1-A4 results: no guide overlap (min Hamming 9-12, none identical), no site overlap (0 test positive pairs in train in any fold); untrained network AUROC on fold 4 = 0.499 (so not an architecture prior); within-guide shuffle on fold 4 gave AUROC 0.618/0.503/0.564 across 3 seeds, not a stable 0.69. The decision rules above therefore give: NOT leakage, architecture explanation NOT supported, remainder UNRESOLVED. A5 tests seed noise: the original-style global label shuffle, fold 4 and fold 0, 8 seeds each (seeds 10-17). Rule: if fold-4 mean AUROC <= 0.58 and the range includes a value <= 0.52, call the original 0.69 seed noise from only 2 test guides; otherwise a persistent fold-specific effect, still unexplained.
