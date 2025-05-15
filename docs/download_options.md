# Download options

You can obtain Portal data either by downloading a single project, creating a custom dataset with a selection of projects and/or samples, or by choosing one of the Portal-wide download options.
For full information about the files you can download from the Portal, see {ref}`the Downloadable Files page<download_files:Downloadable files>`.

This page describes the different options available to you when downloading Portal data.

## Data format

We provide all single-cell and single-nuclei expression data in both `SingleCellExperiment` objects (`.rds` files) for use in R, or as `AnnData` objects (`.h5ad` files) for use in Python.
The default format for all samples on the Portal with single-cell and single-nuclei expression is set to `SingleCellExperiment (R)`.
You can learn more about using these object types from our FAQ sections on {ref}`using the provided RDS files<faq:How do I use the provided RDS files in R?>` and {ref}`using the provided H5AD files<faq:How do I use the provided H5AD files in Python?>`.

Only {ref}`one data format is currently supported for a single download<faq:Why can't I change the data format in My Dataset?>`, including when {ref}`downloading custom datasets<download_files:Custom datasets>`.
To obtain data in both `SingleCellExperiment` and `AnnData` formats, you will need to download these file formats separately.

In addition, note that expression data for multiplexed libraries is only available in `SingleCellExperiment` format, {ref}`as described here<download_files:Multiplexed sample libraries>`.

## Modalities

Besides single-cell/nuclei expression, many samples in the Portal have additional sequencing modalities including CITE-seq, spatial transcriptomics, and bulk RNA-seq.

In particular, there are two modality options that you may see when {ref}`creating a custom dataset to download<download_files:Custom datasets>` or when {ref}`downloading a full project with the "Download Now" button<download_files:Project downloads>`: `Single-cell` and `Spatial`.

By default, the `Single-cell` modality will be selected for all single-cell and single-nuclei RNA-seq samples and/or projects.
Selecting this download option will provide you with the gene expression data from single-cell or single-nuclei samples and/or projects.
If available, CITE-seq expression data will also be included.

If a sample or project has spatial transcriptomic data, you will also have the option to select the `Spatial` modality for download.
Selecting `Spatial` will provide you with the spatial transcriptomic data only.

If you are creating a custom dataset that contains samples and/or projects with bulk RNA-seq data, you will have the option to include this data in your download as well.
Note that the bulk RNA-seq expression file will always include all samples from the given project with bulk expression, even if you are only downloading a subset of that project's samples.
If you are using the "Download now" button to download a full project that contains bulk RNA-seq expression, it will automatically be included with the download.
<!-- TODO: Confirm if there are any Spatial considerations here we need to add:
https://github.com/AlexsLemonade/scpca-docs/pull/413#issuecomment-2867497722 -->

For more information about the expected file download structure for `Single-cell` and `Spatial` modalities, refer to our {ref}`Downloadable files<download_files:STUB LINK TO SECTION WITH SINGLE-CELL/SPATIAL DOWNLOAD FOLDERS>`.

## Merged objects

When downloading a project, either by using `Download Now` or `Add to Dataset`, you will have the option to either receive the data as objects for individual libraries, or as {ref}`a single merged object with data from all samples in the given project<merged_objects:Merged objects>`.
Please be aware that merged objects have _not_ been integrated or batch-corrected.
Refer to {ref}`this documentation<download_files:Merged object downloads` for the contents of a merged object download specifically.
Note that this applies only to `Single-cell` modality downloads, not `Spatial`.

When {ref}`creating a custom dataset to download<download_files:Custom datasets>`, you will be able to select the option to merge all samples only if you have included all samples from the given project in `My Dataset`.
Merging a subset of samples in a project {ref}`is not currently supported<faq:STUB for https://github.com/AlexsLemonade/scpca-docs/issues/399>`.
In addition, merged objects are not available for all samples or projects, {ref}`as described here <faq:Which projects can I download as merged objects?>`.

Note that even when {ref}`downloading data for all single-cell and single-nuclei samples on the Portal<download_files:Portal-wide downloads>`, merged objects will still be provided per-project.
There will not be a merged object with all samples from all projects, but a single merged object for each project.

## Multiplexed sample libraries

When downloading a project that contains multiplexed samples (see {ref}`What is a multiplexed sample? <faq:What is a multiplexed sample?>`), you will have the option to exclude multiplexed samples from the download.
If selected, the download will contain expression data for only non-multiplexed samples.
Note that, {ref}`as described in our FAQ<faq:Which samples can I download as AnnData objects?>`, `AnnData` objects (`.h5ad` files) are not available for multiplexed samples.
In addition, you will not be able to select the option to merge samples into a single file {ref}`if the project contains multiplexed samples<faq:Which projects can I download as merged objects?>`.
