# DESeq2 Benchmarking Scripts

This directory contains scripts for running benchmarking analyses using DESeq2.

## Subdirectories

### `parallel_models`

This subdirectory contains scripts for running DESeq2 differential expression analyses with various model specifications in parallel using SLURM. 

- `single_model.R`: R script for running a single DESeq2 analysis with a specific model. Takes a model index as an argument.
- `single_model.sh`: SLURM script for submitting a single DESeq2 analysis job. Uses `single_model.R`.
- `many_models.sh`: Script for submitting multiple DESeq2 analysis jobs, each with a different model, using SLURM.

## How to Run

1. Navigate to `scripts/benchmarking/deseq/parallel_models/`
2. Run `bash many_models.sh <start_model> <end_model>` to submit multiple jobs. Replace `<start_model>` and `<end_model>` with the range of model indices you want to run.

For example, to run all models from index 1 to 11, run:
```bash
bash many_models.sh 1 11
