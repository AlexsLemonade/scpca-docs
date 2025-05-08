# Download options

You can obtain Portal data either by downloading a single project, creating a custom dataset with a selection of projects and/or samples, or by choosing one of the Portal-wide download options.
For full information about the files you can download from the Portal, see {ref}`the Downloadable Files page<downloadable_files:Downloadable files>`.

This page describes the different options available to you when downloading Portal data.

## Data format

We provide all single-cell and single-nuclei expression data in both `SingleCellExperiment` objects (`.rds` files) for use in R, and as `AnnData` objects (`.h5ad` files) for use in Python.
You can learn more about using these object types from our FAQ sections on {ref}`using the provided RDS files<faq:How do I use the provided RDS files in R?>` and {ref}`using the provided H5AD files<faq:How do I use the provided H5AD files in Python?>`.

Note that only one data format is currently supported for a single download, including when {ref}`downloading custom datasets<downloadable_files:Custom datasets>`.
To obtain data in both `SingleCellExperiment` and `AnnData` formats, you will need to download these file formats separately.

In addition, note that `AnnData` files are not available for projects with multiplexed libraries all data on the Portal, {ref}`as described here:<download_files:Multiplexed sample libraries>`.

## Modalities

Besides single-cell/nuclei expression, many samples in the Portal have additional sequencing modalities including CITE-seq, bulk RNA-seq, and spatial transcriptomics.

When downloading a project with spatial transcriptomics data, you will have the option to select either the "Single-cell" or "Spatial" modality.
Selecting "Single-cell" will provide you with the single-cell/nuclei expression data only, and selecting "Spatial" will provide you with the spatial data only.
If available, CITE-seq data will be automatically included with the "Single-cell" modality only.
Only one modality can be selected at a time.

If the project also has bulk RNA-seq data, you can choose to include it in your download for either modality selection.
