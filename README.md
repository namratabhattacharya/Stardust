# Stardust: Simultaneous Visualisation of Cells and Marker Genes from scRNA-seq Studies

**Stardust** is a Python package for cell-gene co-embedding and downstream analysis of single-cell transcriptomics data. It enables simultaneous representation of cells and marker genes, preserving both transcriptional similarity and inferred spatial proximity.

For thorough details, refer to our paper: [bioRxiv link](https://www.biorxiv.org/content/10.1101/2022.12.27.521966v2)

---

## Table of Contents
- [Environment Requirements](#environment-requirements)
- [Installation Guide](#installation-guide)
- [Supported Dataset Formats](#supported-dataset-formats)
- [Usage](#usage)
- [Tutorial](#vignette-tutorial)
  - [Co-embedding](#cell-gene-co-embedding)
  - [Analysis](#analysis)
- [Output Directory Structure](#output-directory-structure)
- [Contributing](#contributing)

---

## Environment Requirements

- Python 3.7 or above
- `setuptools` must be installed

---

## Installation Guide
  
**Step 1:** Download the package directly from GitHub. For Linux users
```bash
https://github.com/namratabhattacharya/Stardust.git
```
**Step 2:** Download the binary openOrd.jar from the release tag v1.0 and place it under Stardust_package/stardust/run_stardust/

**Step 3:** Change the directory to the package directory
```bash
cd Stardust/
```
**Step 4:** create the environment (Conda ≥4.9).
```bash
conda env create -f environment.yml
```
**step 5:** activate environment
```bash
conda activate stardust37
```

## Supported dataset format
Stardust accepts single-cell expression data in the following formats:

- **CSV Format**  
  A gene expression matrix file named `expression.csv`, with **cells as rows** and **genes as columns**.  
  - *An example dataset is provided in the `example_data/csv/` folder.*

- **10x Genomics Format**  
  Standard `.10x` output files from the 10x Genomics Cell Ranger pipeline.  
  - *An example dataset is provided in the `example_data/10x/` folder.*

## Usage

```python
import stardust

# Generate co-embedding plots
stardust.run_stardust.run()

# Further downstream analysis (Optional)
stardust.analysis.silhouette() #This function calculates the silhouette coefficient for the embedded cells.
stardust.analysis.heatmap() #This function generates a marker gene heatmap.
stardust.analysis.alluvial() #This function creates an alluvial (flow) plot.
```
Outputs are saved under:
**Stardust_results/visualization_output/4_pass/**

## Vignette tutorial
This vignette uses a melanoma data set from the website [here](https://singlecell.broadinstitute.org/single_cell/study/SCP11/melanoma-intra-tumor-heterogeneity) to demonstrate a standard pipeline. This vignette can be used as a tutorial as well.

### cell-gene co-embedding
To reproduce the cell-gene co-embedding visualisation,you need to run the following commands

```python
import stardust

stardust.run_stardust.run() #generated the coembedding plots under Stardust_resuts/visualization_output/4_pass.
```
On run command, you need to provide the expression matrix path and the file format (i.e either "csv" or "10x" ). 
NOTE: For CSV file format, the expression file name should be "expression.csv"
![alt text](blob/run.PNG)

#### co-embedding on Stardust
<p style="text-align: left;">
  <img src="blob/sd_color_embedding.png" alt="Stardust embedding" width="500"/>
</p>


#### co-embedding on UMAP
<p style="text-align: left;">
  <img src="blob/umap_color_embedding.png" alt="Stardust embedding" width="500"/>
</p>

### analysis
NOTE: To run the analysis functions, you need to run the cell-gene co-embedding. 
#### silhouette co-efficient
```python
stardust.analysis.silhouette()
```
<p style="text-align: left;">
  <img src="blob/silhouette.png" alt="Stardust embedding" width="300"/>
</p>

#### marker heatmap
```python
stardust.analysis.heatmap()
```
<p style="text-align: left;">
  <img src="blob/heatmap.png" alt="Stardust embedding" width="500"/>
</p>

#### clustering_metrics.py — quantify clustering quality (NMI, ARI)
Compare Stardust clusters with Scanpy Louvain clusters
```bash
python clustering_metrics.py
```
## Output Directory structure
- **Stardust_results/visualization_output/4_pass**
  - Contains co-embedding plots, CSV outputs, and layout files
- **Stardust_results/analysis/**
  - Contains silhouette scores, marker gene heatmaps

## Contributing
Souce code: [Github](https://github.com/namratabhattacharya/Stardust.git)  
For further information, contact debarka@iiitd.ac.in or n.bhattacharya@qut.eu.au. 
