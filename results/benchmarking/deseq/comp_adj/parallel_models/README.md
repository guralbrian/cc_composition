# DESeq2 Parallel Models Benchmarking Results

This directory contains the results from running DESeq2 differential expression analyses with various model specifications in parallel using SLURM.

## Subdirectories and Files

- `clr_4`, `clr_cm`, `clr_cm_fb`, `clr_cm_fb_im`: Results directories for CLR-adjusted models.
- `raw_4`, `raw_cm`, `raw_cm_fb`, `raw_cm_fb_im`, `unadjusted`: Results directories for unadjusted and raw models.
- `pc1`, `pc2`: Results directories for principal component adjusted models.
- `counts_10000g.csv`: Counts matrix used for the analysis.
- `counts_dict_10000g.csv`: Dictionary file for counts matrix.
- `ratios_10000g.csv`: Ratios file.

## Metrics

Metrics such as false positives, true positives, etc., are saved within the respective subdirectories.

To understand the different metrics and models used, refer to the README in the `scripts/benchmarking/deseq/parallel_models/` directory.