# Data Loading Scripts

This directory contains scripts for loading datasets used in this project. All of the datasets described here are publically available scRNAseq from mouse hearts. Two versions of each script are present, one for workflows within RStudio processes, and another for Snakemake pipelines (denoted with the _snakemake suffix). 

## Datasets

### Dataset 1

- **Script**: `tabula_muris.R`
- **Paper**: [Single-cell transcriptomics of 20 mouse organs creates a Tabula Muris](https://www.nature.com/articles/s41586-018-0590-4)
- **Citation**: The Tabula Muris Consortium., Overall coordination., Logistical coordination. et al. Single-cell transcriptomics of 20 mouse organs creates a Tabula Muris. Nature 562, 367–372 (2018). https://doi.org/10.1038/s41586-018-0590-4
- **Summary**: scRNAseq of 53,760 cells from 20 organs and 8 mices, male and female. We're interested in the heart tissue cells. The goal is to increase the diversity in data sources within each cell type with the idea that it will result in more robust deconvolution marker selection. 
- **Data Description**: Currently, Tabula Muris droplet data accessed by `TabulaMurisSenisData::TabulaMurisSenisDroplet`(), but this project will eventually switch to accessing and processing the raw read files.

### Dataset 2

- **Script**: `wu.R`
- **Paper**: [Single-cell RNA sequencing of mouse left ventricle reveals cellular diversity and intercommunication](doi.org/10.1152/physiolgenomics.00016.2021)
- **Citation**: Wu, L., Li, Y., Shen, J., Zhu, Q., Jiang, J., Ma, S., … & Li, X. (2022). Single-cell rna sequencing of mouse left ventricle reveals cellular diversity and intercommunication. Physiological Genomics, 54(1), 11-21. https://doi.org/10.1152/physiolgenomics.00016.2021
- **Summary**: Wu et al isolated left ventricles from 12wk old C57BL/6J mice and sequences with 10x Chromium and Hiseq. 
- **Data Description**: They maintained 8,527 cells, filtered >10% mitocondrial read-mapped cells, and > +/- 2-fold standard deviations. This dataset is incredibly enriched for macrophages/immune cells. 

### Dataset 3

- **Script**: `martini.R`
- **Paper**: [Single-Cell Sequencing of Mouse Heart Immune Infiltrate in Pressure Overload–Driven Heart Failure Reveals Extent of Immune Activation](doi.org/10.1161/CIRCULATIONAHA.119.041694)
- **Citation**: Martini, E., Kunderfranco, P., Peano, C., Carullo, P., Cremonesi, M., Schorn, T., Carriero, R., Termanini, A., Colombo, F. S., Jachetti, E., Panico, C., Faggian, G., Fumero, A., Torracca, L., Molgora, M., Cibella, J., Pagiatakis, C., Brummelman, J., Alvisi, G., Mazza, E. M. C., … Kallikourdis, M. (2019). Single-Cell Sequencing of Mouse Heart Immune Infiltrate in Pressure Overload-Driven Heart Failure Reveals Extent of Immune Activation. Circulation, 140(25), 2089–2107. https://doi.org/10.1161/CIRCULATIONAHA.119.041694
- **Summary**: Martini et al mapped the cardiac immune composition in the standard murine nonischemic, pressure-overload heart failure mode. They focused  on CD45+ cells, looking for immune cell subsets in the heart, at early and late stages of disease and in controls.
- **Data Description**: Initial quality control in each biological replicate was assessed, by filtering out cells meeting any of the following criteria: less than 200 or more than 5,000 unique genes expressed, more than 20,000 UMIs, or more than 12.5% of reads mapping to mitochondria. Final dataset included 17,853 cardiac CD45+ (immune) cells


Certainly! Below is a draft README file for your project. Feel free to modify it as needed.

---

# Cardiac Cell Composition Analysis

Written by [Brian Gural]{https://www.linkedin.com/in/brian-gural-09bb60128/}
Last updated to v0.2.0 on September 5th 2023

## Project Goals

The primary aim of this project is to analyze the cellular composition of cardiac tissue using single-cell RNA sequencing data as a reference. Data is being generated and analyzed by the [Rau Lab at the University of North Carolina at Chapel Hill]{https://raulab.web.unc.edu/}. We are currently exploring pilot data, consisting of bulk RNAseq from the left ventricles of isoproterenol-treated mice from 16 strains of the Collaborative Cross. For reference-based deconvolution, we are using a combination of internal snRNAseq (LV, untreated, C57B/L6) and publically available scRNAseq (described in the relevant script sub-directory){https://github.com/guralbrian/cc_composition/tree/main/scripts/load_data/single_cell}. 

**This analysis pipeline currently is able to:**

- [Download]{https://github.com/guralbrian/cc_composition/tree/main/scripts/load_data/} and [preprocess]{https://github.com/guralbrian/cc_composition/tree/main/scripts/quality_control} the raw data
- [Merge sn/scRNAseq references and identify cell-type specific marker genes]{https://github.com/guralbrian/cc_composition/blob/main/scripts/analysis/merge_sn_old.Rmd}
- [Perform reference-based deconvolution]{https://github.com/guralbrian/cc_composition/blob/main/scripts/analysis/deconvolution_old.Rmd} 
- [Visualize compositions and perform/visualize dimensional reductions]{Phttps://github.com/guralbrian/cc_composition/blob/main/scripts/visualization/06222023/composition/comp_and_pca.Rmd}

**Immediate next steps for actual analysis include:**

- Unblind samples 
- Validate compositional estimates with ground truths
    - Ground truth one: Bulk RNAseq of experimentally purified major cell types
    - Ground truth two: Pseudobulk generated by summed expression from sn/scRNAseq
- Estimate phenotype and genotype effects on composition by Dirichlet regression
- Differential expression analysis adjusted w/ compositional covariates

**Further analysis include:**
- QTL mapping for compositional changes during isoproterenol treatment
- Composition-adjusted QTL mapping for heart failure phenotypes
- Further ground-truths after benchwork is completed. Involves *in-vitro* mixtures of RNA isolated from purified cardiac cell types, sequenced by bulk RNAseq *a la* [Tumor Deconvolution DREAM Challenge]{https://www.synapse.org/#!Synapse:syn15589870/wiki/592683}

## Expected changes to analysis infrastructure

This project is also meant to be a learning process for reproducible data analysis. Currently, all scripts are .R or .Rmd and meant to be manually opened and run in RStudio. Libraries and languages are not version controlled. Code is being pushed to a working branch and merged to main upon new, personally defined, functional versions. All code, data, and results are stored on UNC's HPC, Longleaf. 

I intend to implement the following:
    1. Snakemake automation for running analysis scripts
        - Will need to adjust .R scripts to accept Snakemake options 
        - .RMD files will need to be converted to .R files 
    2. Choose and implement a package version control system. 
        - Lightweight options include Conda and Renv
        - More robust and intensive options include Docker and Singularity 
            - Docker commands aren't natively supported by Longleaf
            - Singlularity may be able to provide a workaround.
    3. Interactive data visualization tools, like Shiny.
    4. Automated storage of slurm.out files
        - Also would like to merge .out files with R .log files for relevant scripts


## Directory Structure

```
.
├── data
│   ├── raw
│       ├── cc_counts
│       └── rau_sn
│   └── processed
│       ├── single_cell
│       └── bulk
├── scripts
│   └── load_data
│       └── single_cell
│           ├── wu.R
│           ├── martini.R
│           └── tabula_muris.R
├── snakemake
│   └── target_files
├── results
│   └── UNFINISHED
├── tools
│   └── UNFINISHED
├── logs
│   └── slurm
├── cc_composition.Rproj
├── environment.yml
├── Snakefile
└── README.md
```

- `data/raw`: Contains the raw 10x outputs and unprocessed counts matrix for bulk RNAseq (processed by salmon).
- `data/processed/single_cell`: Contains the processed data in h5 format.
- `scripts/load_data/single_cell`: Contains R scripts for data downloading and preprocessing.
- `snakemake/target_files`: Contains Snakemake target files. 
- `logs/slurm`: Contains SLURM log files.

