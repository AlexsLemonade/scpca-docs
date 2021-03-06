# Frequently Asked Questions

## Why did we use Alevin-fry for processing?

We aimed to process all of the data in the portal such that it is comparable to widely used pipelines, namely Cell Ranger from 10x Genomics.
In our own benchmarking, we found that [Alevin-fry](https://github.com/COMBINE-lab/alevin-fry) produces very similar results to [Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count), while allowing faster, more memory efficient processing of single-cell and single-nuclei RNA-sequencing data.
In the configuration that we are using ("selective alignment" mapping to a human transcriptome that includes introns), Alevin-fry uses approximately 12-16 GB of memory per sample and completes mapping and quantification in less than an hour.
By contrast, Cell Ranger uses up to 25-30 GB of memory per sample and takes anywhere from 2-8 hours to align and quantify one sample.
Quantification of samples processed with both Alevin-fry and Cell Ranger resulted in similar distributions of mapped UMI count per cell and genes detected per cell for both tools.

![](https://github.com/AlexsLemonade/alsf-scpca/blob/c0c2442d7242f6e06a5ac6d1e45bd1951780da14/analysis/docs-figures/plots/total_umi_per_cell.png?raw=true)

![](https://github.com/AlexsLemonade/alsf-scpca/blob/c0c2442d7242f6e06a5ac6d1e45bd1951780da14/analysis/docs-figures/plots/total_genes_per_cell.png?raw=true)

We also compared the mean gene expression reported for each gene by both methods and observed a high correlation with a Pearson R correlation coefficient of 0.98.

![](https://github.com/AlexsLemonade/alsf-scpca/blob/c0c2442d7242f6e06a5ac6d1e45bd1951780da14/analysis/docs-figures/plots/gene_exp_correlation.png?raw=true)

Recent reports from others support our findings.
[He _et al._ (2021)](https://doi.org/10.1101/2021.06.29.450377)) demonstrated that Alevin-fry can process single-cell and single-nuclei data more quickly and efficiently then other available methods, while also decreasing the false positive rate of gene detection that is commonly seen in methods that utilize transcriptome alignment.
[You _et al._ (2021)](https://doi.org/10.1101/2021.06.17.448895) and [Tian _et al._ (2019)](https://doi.org/10.1038/s41592-019-0425-8) have also noted that results from different pre-processing workflows for single-cell RNA-sequencing analysis tend to result in compatible results downstream.

## How do I use the provided RDS files in R?

We are providing the gene expression data to you as a [`SingleCellExperiment` object](http://bioconductor.org/books/3.13/OSCA.intro/the-singlecellexperiment-class.html) in an RDS file.

_Note: You will need to install and load the [`SingleCellExperiment` package](https://bioconductor.org/packages/3.13/bioc/html/SingleCellExperiment.html) from Bioconductor to work with the provided files._

To read in the RDS files you can use the `readRDS` command in base R.

```r
library(SingleCellExperiment)
scpca_sample <- readRDS("SCPCL000000_filtered.rds")
```

## What is the difference between samples and libraries?

A sample ID, labeled as `scpca_sample_id` and indicated by the prefix `SCPCS`, represents a unique tissue that was collected from a participant.

The library ID, labeled as `scpca_library_id` and indicated by the prefix `SCPCL`, represents a single set of cells from a tissue sample, or a particular combination of samples in the case of multiplexed libraries.
For single-cell or single-nuclei experiments, this will be the result of emulsion and droplet generation using the 10X Genomics workflow, potentially including both RNA-seq, CITE-seq and cell hashing sequencing libraries.
Multiplexed libraries will have more than one sample ID corresponding to each library ID.

In most cases, each sample will only have one corresponding single-cell or single-nuclei library, and may also have an associated bulk RNA-seq library.
However, in some cases multiple libraries were created by separate droplet generation and sequencing from the same sample, resulting in more than one single-cell or single-nuclei library ID being associated with the same sample ID.

## Why do some samples have missing participant IDs?

The `participant_id`, when present, indicates the participant from which a collection of samples was obtained.
For example, one participant may have a sample collected both at initial diagnosis and at relapse.
This would result in two different sample ID's, but the same participant ID.
However, for most participants, only a single sample was collected and submitted for sequencing.
Because of this, many of the samples do not have a separate participant ID.
Participant IDs are only present for samples that were derived from the same participant as at least one other sample.

## What is a multiplexed sample?

Multiplexed samples refer to samples that have been combined together into a single library using cell hashing ([Stoeckius _et al._ 2018](https://doi.org/10.1186/s13059-018-1603-1)) or a related technology and then sequenced together.
This means that a single library contains cells or nuclei that correspond to multiple samples.
Each sample has been tagged with a hashtag oligo (HTO) prior to mixing, and that HTO can be used to identify which cells or nuclei belong to which sample within a multiplexed library.
The libraries available for download on the portal have not been separated by sample (i.e. demultiplexed), and therefore contain data from multiple samples.
For more information on working with multiplexed samples, see the {ref}`special considerations for multiplexed samples section in getting started with an ScPCA dataset <getting_started:Special considerations for multiplexed samples>`.

## Why are demultiplexed samples not available?

Downloading a multiplexed sample on the portal will result in obtaining the gene expression files corresponding to the library containing the chosen multiplexed sample and any other samples that were multiplexed with that chosen sample.
This means that users will receive the gene expression data for all samples that were combined into a given library and will have to separate any cells corresponding to the sample of interest before proceeding with downstream analysis.

We have applied multiple {ref}`demultiplexing methods <processing_information:hto demultiplexing>` to multiplexed libraries and noticed that these demultiplexing methods can vary both in calls and confidence levels assigned.
[Here we have performed some exploratory analysis comparing demultiplexing methods within a single multiplexed library](https://htmlpreview.github.io/?https://github.com/AlexsLemonade/alsf-scpca/blob/main/analysis/quantifier-comparisons/15-demux-comparisons.nb.html).
Because of the inconsistency across demultiplexing methods used, the choice of demultiplexing method to use is up to the discretion of the user.
Rather than separating out each sample, the sample calls and any associated statistics regarding sample calls for multiple demultiplexing methods can be found in the `_filtered.rds` file for each multiplexed library.
See the {ref}`demultiplexing results section <sce_file_contents:demultiplexing results>` for instructions on how to access the demultiplexing results in the `SingleCellExperiment` objects for multiplexed libraries.
We also include the Hash Tag Oligo counts matrix to allow demultiplexing using other available methods.

## What are estimated demux cell counts?

Estimated demux cell counts are provided for multiplexed libraries and refer to the estimated cell counts for each sample that is present in the library.
In order to provide an estimate of the number of cells or nuclei that are present in a given sample before download, we use the estimated number of cells per sample identified by one of the tested demultiplexing methods.
However, these estimated demux cell counts should only be considered a guide; we encourage users to investigate the data on their own and make their own decisions on the best demultiplexing method to use for their research purposes.

Note that not all cells in a library are included in the estimated demux cell count, as some cells may not have been assigned to a sample.
Estimated demux cell counts are only reported for multiplexed samples and are not reported for single-cell or single-nuclei samples that are not multiplexed.
For more about demultiplexing, see the section on {ref}`processing multiplexed libraries <processing_information:multiplexed libraries>`.

## What genes are included in the reference transcriptome?

The {ref}`reference transcriptome index <processing_information:reference transcriptome index>` that was used for alignment was constructed by extracting both spliced cDNA and intronic regions from the primary genome assembly GRCh38, Ensembl database version 104 ([see the code used to generate the reference transcriptome](https://github.com/AlexsLemonade/scpca-nf/blob/main/bin/make_splici_fasta.R)).
The resulting reference transcriptome index contains 60,319 genes.
In addition to protein-coding genes, this list of genes includes pseudogenes and non-coding RNA.
The gene expression data files available for download report all possible genes present in the reference transcriptome, even if not detected in a given library.

## Where can I see the code for generating QC reports?

A QC report for every processed library is included with all downloads, generated from the unfiltered and {ref}`filtered <processing_information:filtering cells>` {ref}`Single-cell gene expression files <sce_file_contents:Single-cell gene expression file contents>`.
You can find the [function for generating a QC report](https://github.com/AlexsLemonade/scpcaTools/blob/main/R/generate_qc_report.R) and the [QC report template documents](https://github.com/AlexsLemonade/scpcaTools/tree/main/inst/rmd) in the package we developed for working with processed ScPCA data, [`scpcaTools`](https://github.com/AlexsLemonade/scpcaTools).

## What if I want to use Seurat instead of Bioconductor?

The files available for download contain [`SingleCellExperiment` objects](http://bioconductor.org/books/3.13/OSCA.intro/the-singlecellexperiment-class.html).
If desired, these can be converted into Seurat objects.

You will need to [install and load the `Seurat` package](https://satijalab.org/seurat/articles/install.html) to work with Seurat objects.

For libraries that only contain RNA-sequencing data (i.e. do not have a CITE-seq library found in the `altExp` of the `SingleCellExperiment` object), you can use the following commands:

```r
library(Seurat)
library(SingleCellExperiment)

# read in RDS file
sce <- readRDS("SCPCL000000_filtered.rds")

# create seurat object from the SCE counts matrix
seurat_object <- CreateSeuratObject(counts = counts(sce),
                                    assay = "RNA",
                                    project = "SCPCL000000")
```
The above code will only maintain information found in the original counts matrix from the `SingleCellExperiment`.
Optionally, if you would like to keep the included cell and gene associated metadata during conversion to the Seurat object you can perform the below additional steps:

```r
# convert colData and rowData to data.frame for use in the Seurat object
cell_metadata <- as.data.frame(colData(sce))
row_metadata <- as.data.frame(rowData(sce))

# add cell metadata (colData) from SingleCellExperiment to Seurat
seurat_object@meta.data <- cell_metadata

# add row metadata (rowData) from SingleCellExperiment to Seurat
seurat_object[["RNA"]]@meta.features <- row_metadata

# add metadata from SingleCellExperiment to Seurat
seurat_object@misc <- metadata(sce)
```

For `SingleCellExperiment` objects from libraries with both RNA-seq and CITE-seq data, you can use the following additional commands to add a second assay containing the CITE-seq counts and associated feature data:

```r
# create assay object in Seurat from CITE-seq counts found in altExp(SingleCellExperiment)
cite_assay <- CreateAssayObject(counts = counts(altExp(sce)))

# optional: add row metadata (rowData) from altExp to assay
cite_row_metadata <- as.data.frame(rowData(altExp(sce)))
cite_assay@meta.features <- cite_row_metadata

# add altExp from SingleCellExperiment as second assay to Seurat
seurat_object[["CITE"]] <- cite_assay
```

## What if I want to use Python instead of R?

We provide single-cell and single-nuclei gene expression data as RDS files, which must be opened in R to view the contents.
If you prefer to work in Python, there are a variety of ways of converting the count data to Python-compatible formats.
We have found that one of the more efficient is conversion via the 10X format using [`DropletUtils::write10xCounts()`](https://rdrr.io/bioc/DropletUtils/man/write10xCounts.html).
Note that you will need to install the [`DropletUtils` package](https://www.bioconductor.org/packages/devel/bioc/html/DropletUtils.html) to use this function.

When used as described below, `DropletUtils::write10xCounts()` will output three files to a new directory, following the format used by Cell Ranger 3.0 (and later):
- the counts matrix in sparse matrix format - `matrix.mtx.gz`
- the row names, or gene names, saved as a TSV - `features.tsv.gz`
- the column names, or cell barcodes, saved as a TSV - `barcodes.tsv.gz`

```r
library(SingleCellExperiment)

# read in the RDS file to be converted
sce <- readRDS("SCPCL000000_filtered.rds")

# write counts to 10X format and save to a folder named "SCPCL000000-rna"
DropletUtils::write10xCounts("SCPCL000000-rna", counts(sce),
                             barcodes = colnames(sce),
                             gene.id = rownames(sce),
                             gene.symbol = rowData(sce)$gene_symbol,
                             version = "3")

```

If a library has associated CITE-seq that exists, you will have to save that separately.

```r
# write CITE-seq counts to 10X format
DropletUtils::write10xCounts("SCPCL000000-cite", counts(altExp(sce)),
                             barcodes = colnames(altExp(sce)),
                             gene.id = rownames(altExp(sce)),
                             gene.type = "Antibody Capture",
                             version = "3")

```

These files can then be directly read into Python using the [`scanpy` package](https://scanpy.readthedocs.io/en/stable/), creating an [`AnnData` object](https://anndata.readthedocs.io/en/latest/index.html).
Note that you will need to [install the `scanpy` package](https://scanpy.readthedocs.io/en/stable/installation.html).

```python
import scanpy as sc

#read in 10X formatted files
rna_file_directory = "SCPCL000000-rna"
anndata_object = sc.read_10x_mtx(rna_file_directory)

cite_file_directory = "SCPCL000000-cite"
cite_anndata = sc.read_10x_mtx(cite_file_directory)

# append CITE-seq anndata to RNA-seq anndata
anndata_object['CITE'] = cite_anndata.to_df()
```

It should be noted that in this conversion the `colData`, `rowData`, and metadata that are found in the original `SingleCellExperiment` objects will not be retained.
If you would like to include this data, you could write out each table separately and load them manually in Python.
Alternatively, you might be interested in this [reference from the authors of `scanpy`](https://theislab.github.io/scanpy-in-R/#converting-from-r-to-python) discussing a different approach to  conversion using Rmarkdown notebooks and the `reticulate` package to directly convert `SingleCellExperiment` object components to `AnnData` object components without writing files locally.
