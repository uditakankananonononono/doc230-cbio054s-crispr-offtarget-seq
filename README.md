partial: 3D core blocked (no coordinates in public encoded data) Audit: no leakage found (guide + site overlap 0, RESULTS.md); fold-4 shuffle AUROC residual (mean 0.61) unexplained. Extended 2026-10-09 with E6-E10 (mismatch strata, per-guide, training-size curve, calibration, guide-shuffled control; RESULTS.md): sequence-only, single seed.

# doc230-cbio054s-crispr-offtarget-seq
Sequence-only CRISPR off-target benchmark (5 directions) for parent CBIO054 "3D-Aware CRISPR Off-Target Prediction". No 3D/chromatin features used or claimed. See PREREG.md, RESULTS.md, src/run.py, data_manifest.

## Reproduce / readability audit (2026-10-09)
Data: sources and SHA256 are in data_manifest/. Scripts read fixed paths under /tmp (download the files there first, check the SHA256, then run python3 src/run.py and python3 src/audit.py). Not independently re-run from a clean machine as part of this readability check; the numbers in RESULTS.md come from the original runs on 2 CPU / 2 GB. Three CSVs from bio.info.uqam.ca go in /tmp/crispr/.
