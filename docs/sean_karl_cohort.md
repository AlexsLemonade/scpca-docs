# The Sean Karl cohort 

This page contains information about additional processing specific to the Sean Karl cohort, which corresponds to samples available as part of `SCPCP000026`. 
In collaboration with the submitters, the Data Lab developed an open-source Nextflow pipeline, [`ews-nf`](https://github.com/AlexsLemonade/ews-nf), which was used to annotate all cells and identify recurrent gene expression programs, or metaprograms, present across all samples in the cohort. 
The results from this workflow have been incorporated into all samples available for download. 

## Custom cell type annotations 

Custom cell type annotations can be found in the `submitter_celltype_annotation` column of the cell metadata (`colData` for `SingleCellExperiment` objects or `obs` for `AnnData` objects). 
To generate these custom cell type annotations, tumor cells were identified using `ews-nf` with the default annotation workflow (`--workflow annotation`). 
Please see the [`ews-nf` documentation](https://github.com/AlexsLemonade/ews-nf/blob/main/README.md) for more information on how tumor cells were classified. 

All immune cells (identified by their {ref}`consensus cell type annotations <processing_information:Cell type annotation>`) were then further refined by the submitters.

## Metaprograms 

Recurrent gene expression programs, or metaprograms, were identified across all tumor cells and all samples within the Sean Karl cohort using the metaprogram workflow within `ews-nf` (`--workflow nmf_metaprograms`). 
Briefly, `cNMF` ([Kotlier _et al._ (2019)](https://doi.org/10.7554/eLife.43803)) was used to run non-negative matrix factorization (NMF) on each sample across a range of ranks. 
The resulting NMF programs across all samples were then clustered into metaprograms based on their Pearson correlation coefficient. 
This process was repeated across a wide range of k values (number of clusters) and a set of unbiased metrics were calculated and used to determine the optimal value of k, corresponding to the final number of metaprograms. 

Each resulting metaprogram contains a list of genes and gene weights that were then used to score individual cells. 
Scores were calculated by quantile-normalizing the raw counts matrix and then taking the dot product between the quantile-normalized counts and each metaprogram. 

The objects available on the Portal contain the scores for all cells across all metaprograms. 
The scores for each metaprogram are present as an individual column in the cell metadata (`colData` for `SingleCellExperiment` objects or `obs` for `AnnData` objects). 
<!--TODO: Update the naming scheme here with anything decided in https://github.com/AlexsLemonade/scpca-nf/pull/1307-->
Column names are formatted with: `submitter_data_<metaprogram name>_score`, where the `metaprogram name` is the designated name provided by the submitters. 

For more detailed information on how metaprograms were generated, see the [`ews-nf` documentation](https://github.com/AlexsLemonade/ews-nf/blob/main/README.md)
