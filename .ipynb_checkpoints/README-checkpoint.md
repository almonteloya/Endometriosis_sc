## Altered Endometrial Cellular Features in Endometriosis Patients Identified by Single-Cell Analysis Across the Menstrual Cycle

### Getting Started

**Data**: 10x Single-cell RNA-seq endometrium tissue samples. From patients with endometriosis, patients with fibroids and pristine controls. 

This repository included the instruction for the prepossessing of the data and code for downstream analysis. Each pool containg a single library. 


## Preprocessing: 


#### Generation of count matrix:

Thirteen pools were independently processed using the Cell Ranger software pipeline (version 3.1.0). Reads were mapped to the transcriptome using the human hg38 version, and UMI counts were calculated using default parameters.

We processed each UMI count matrix using the R package Seurat (v.4.0.3). A single Seurat object was created for each pool.

#### Demultiplexing:

The samples were demultiplexed using Freemuxlet using the number of samples per pool and a list of valid barcodes. The list of valid barcodes was obtained from each object.  

Using Freemuxlet cells were classified according to their UMI content and doublet information.


#### Filtering: 

Cells classified as singlets by freemuxlet were retained. Further filtering per each pool included removing cells with  high mitochondrial content (>20%), high ribosomal content (>50%), high RNA count (50,000>), low nfeature  (<200) and low RNA count (<300).  

#### Normalization and dimensionality reduction:

We  performed normalization and scaling using default parameters. We calculated variable features (n=3000), principal component analysis (PCA), and SCT Normalization (regressing for cell cycle, ribosomal and mitochondrial genes). Followed by RunUMAP (dims=1:30), FindNeighbors and FindClusters (Leiden algorithm). 



#### Sample assignation:

To assign demultiplexed cells to the sample of origin, genetic profiles between demultiplexed cells and bulk RNA profiles were compared using bcftools (v.1.10). 
Based on this metadata was assigned: Patient_id, disease_condition and menstrual phase. 


#### Result:

13 Seurat objects, each one with: raw counts, normalization, scaling. PCA and UMAP reductions. 


### Integration:


#### Second Filtering: 

We filtered for cells with high mitochondrial content and samples ambiguously assigned by Freemuxlet

Please see: Code/Integration/Filter_pools/


Note: We conducted an exploration using DoubletFinder; however, the cells identified as doublets were determined to be immune cells of interest, which led to the decision not to employ this method. 
Instead we used FindMarkers and filtered clusters which top markers were mitochondrial genes. 


#### Harmony integration: 

Samples were merged and intefrated via harmony. Dimensionality reduction and clustering at different resolutions was performed. We performed another round of filtering. 

Please see: Code/Integration/Harmony_integration.R

## Downstream analysis:

#### SingleR annotation: 

Broad annotation base on: Single-cell transcriptomic atlas of the human endometrium during the menstrual cycle (2020) by Wanxin et al.

See: Code/Downstream_analysis/SingleR
