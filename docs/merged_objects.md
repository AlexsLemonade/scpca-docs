# Merged objects

These objects contain processed gene expression data, cell and gene metrics, associated metadata, and, in the case of multimodal data like ADTs from CITE-seq experiments, data from additional cell-based assays for all libraries in a given ScPCA project.
See {ref}`the section on merged object processing <processing_information:merged objects` for more information on how these objects were prepared.

Merged objects are provided in two formats:
  - As an RDS file containing a [`SingleCellExperiment` object](http://bioconductor.org/books/3.17/OSCA.intro/the-singlecellexperiment-class.html) for use in R.
  - As an HDF5 file containing an [`AnnData` object](https://anndata.readthedocs.io/en/latest/index.html) for use in Python.


Below we present some details about the specific contents of the objects we provide.

## Components of a SingleCellExperiment merged object

To begin, you will need to load the `SingleCellExperiment` package and read the RDS file:

```r
library(SingleCellExperiment)
merged_sce <- readRDS("SCPCP000000_merged.rds")
```

### SingleCellExperiment expression counts

The `counts` and `logcounts` assays of the `SingleCellExperiment` object for single-cell and single-nuclei experiments contains the primary RNA-seq expression data.
The `counts` assay contains raw counts represented as integers, and the `logcounts` assay contains normalized (on a per-library basis) counts as described in {ref}`the data post-processing section <processing_information:processed gene expression data>`.

Both assays include reads aligned to both spliced and unspliced cDNA (see the section on {ref}`Post Alevin-fry processing <processing_information:post alevin-fry processing>`).
The data is stored as a sparse matrix, where each column represents a cell or droplet, and each row represents a gene.
The `counts` and `logcounts` assays can be accessed with the following R code:

```r
counts(merged_sce) # counts matrix
logcounts(merged_sce) # logcounts matrix
```

Column names are cell barcode sequences prefixed with the originating library id, e.g. `SCPCL000000-{barcode}`, and row names are Ensembl gene IDs.
These names can be accessed with the following R code:

```r
colnames(merged_sce) # matrix column names
rownames(merged_sce) # matrix row names
```

There is an additional `spliced` assay which contains the counts matrix with only reads from spliced cDNA:

```r
assay(merged_sce, "spliced") # spliced-only counts matrix
```





## Components of an AnnData merged object

Before getting started, we highly encourage you to familiarize yourself with the general `AnnData` object structure and functions available as part of the [`AnnData` package](https://anndata.readthedocs.io/en/latest/index.html).
For the most part, the `AnnData` objects that we provide are formatted to match the expected data format for [`CELLxGENE`](https://cellxgene.cziscience.com/) following [schema version `3.0.0`](https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/3.0.0/schema.md).

To begin, you will need to load the `AnnData` package and read the HDF5 file:

```python
import anndata
merged_adata_object = anndata.read_h5ad("SCPCP000000_merged_rna.hdf5")
```

### AnnData expression counts

The data matrix `raw.X` of the merged `AnnData` object contains the RNA-seq expression data as primary integer counts, and the data matrix `X` contains the RNA-seq expression data as normalized counts.
The data is stored as a sparse matrix, where each column represents a cell or droplet, and each row represents a gene.
The `raw.X` and `X` matrices can be accessed with the following python code:

```python
merged_adata_object.raw.X # raw count matrix
merged_adata_object.X # normalized count matrix
```

Column names are cell barcode sequences prefixed with the originating library id, e.g. `SCPCL000000-{barcode}`, and row names are Ensembl gene IDs.
These names can be accessed as with the following python code:


```python
merged_adata_object.obs_names # matrix column names
merged_adata_object.var_names # matrix row names
```

