# Data Loading Scripts

This directory contains scripts for loading datasets used in this project. All of the datasets described here are publically available scRNAseq from mouse hearts. 

## Datasets

### Dataset 1

- **Script**: `tabula_muris.R`
- **Paper**: [Single-cell transcriptomics of 20 mouse organs creates a Tabula Muris](https://www.nature.com/articles/s41586-018-0590-4)
- **Citation**: The Tabula Muris Consortium., Overall coordination., Logistical coordination. et al. Single-cell transcriptomics of 20 mouse organs creates a Tabula Muris. Nature 562, 367–372 (2018). https://doi.org/10.1038/s41586-018-0590-4
- **Summary**: scRNAseq of 53,760 cells from 20 organs and 8 mices, male and female. We're interested in the heart tissue cells. The goal is to increase the diversity in data sources within each cell type with the idea that it will result in more robust deconvolution marker selection. 
- **Data Description**: Currently, Tabula Muris droplet data accessed by TabulaMurisSenisData::TabulaMurisSenisDroplet(), but this project will eventually switch to accessing and processing the raw read files.

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