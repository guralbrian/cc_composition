# Snakefile skeleton

rule load_single_cell_data:
    output:
        h5_file = "data/processed/single_cell/{dataset}_data.h5"
    params:
        script = "scripts/load_data/single_cell/{dataset}.R",
        outpath = "data/processed/single_cell/{dataset}_data.h5"
    shell:
        """
        Rscript {params.script} --outpath {params.outpath}
        """
