Markdown has an overview on the basics of scRNA-seq analysis and the code is very generic as a means to explain with example. The actual code is present in each file mentioned at the respective stages of the analysis.

Highly recommend: [Seurat Guided Tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial#assigning-cell-type-identity-to-clusters)

Great Bioinformatics YouTube channel that has content beyond Seurat: [Collection of Online Tutorials](https://www.youtube.com/@Collection_of_online_tutorials/featured)

# 1. Introduction to Single-Cell RNA Sequencing

Single-cell RNA sequencing (scRNA-seq) is a technique used to measure the gene expression profiles of individual cells. Unlike bulk RNA-seq, which provides an average signal across many cells, scRNA-seq captures the cellular heterogeneity within a tissue or sample. This enables the discovery of distinct cell types, rare cell populations, and dynamic changes in cell states.

The scRNA-seq workflow involves several key steps: data generation, quality control, normalization, dimensionality reduction, clustering, differential expression, and downstream biological interpretation.


# 2. Overview of Seurat

Seurat is an R package designed for the analysis and exploration of single-cell RNA sequencing data. The core of Seurat is the Seurat object, which stores raw data, processed data, metadata, and results. Understanding the Seurat object is essential for navigating the workflow efficiently.

Scanpy is a comparable Python framework that uses the same principles but stores data in an `anndata` object. Both frameworks implement similar workflows; choice depends on the preferred programming language.


# 3. Data Import and Quality Control (`01_non_integrated.Rmd`)

## 3.1. Reading Data

Single-cell data is typically generated using droplet-based technologies like 10x Genomics. The output includes a gene-cell count matrix. This can be read into R using functions such as `Read10X()` and then converted into a Seurat object using `CreateSeuratObject()`. The Seurat object serves as the container for all downstream analyses.

```r
data <- Read10X(data.dir = "path/to/data/")
seurat_object <- CreateSeuratObject(counts = data)
```

## 3.2. Quality Metrics

Quality control ensures that only high-quality cells are retained for analysis. Typical quality metrics include:

- The number of detected genes per cell (nFeature_RNA): very low counts may indicate empty droplets or damaged cells.
- The total number of molecules per cell (nCount_RNA): extremely high counts may suggest doublets (two cells in one droplet).
- The proportion of reads mapping to mitochondrial genes: high mitochondrial content can be a sign of stressed or dying cells.

```{r}
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")
```

## 3.3. Filtering Cells
Cells that do not meet quality thresholds are removed using the `subset()` function. Thresholds are adjusted based on plots such as violin plots or scatter plots that show the distribution of these metrics. There is no rule as to what cutoff to use, generally a percent.mt cutoff of 10 is suggested, though many people go as high as 30-35.

```{r}
seurat_object <- subset(seurat_object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
```

# 4. Normalization and Identification of Variable Features
## 4.1. Normalization
Normalization accounts for differences in sequencing depth between cells. A common method is log-normalization, performed using `NormalizeData()`. This scales gene expression so that comparisons between cells are meaningful.

```{r}
seurat_object <- NormalizeData(seurat_object)
```

An alternative method, `SCTransform()`, models technical noise and performs variance stabilization, which can provide more robust results for complex datasets. This however takes a long time to run so is not recommended for large datasets. Both methods perform well based on the scenario, but the downstream processing for `SCTransform()` includes additional steps, which makes a lot of users lean towards `NormalizeData()` instead.

```{r}
seurat_object <- SCTransform(seurat_object)
```

## 4.2. Identifying Highly Variable Genes
Not all genes carry useful information for distinguishing cell types. `FindVariableFeatures()` identifies genes with high variance across cells, which are most informative for dimensionality reduction and clustering.

Visual inspection of variable genes helps verify that highly expressed housekeeping genes do not dominate the list. 2000 variable genes is the general norm.

```{r}
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)
```

# 5. Scaling and Dimensionality Reduction
## 5.1. Scaling Data
`ScaleData()` centers and scales the expression values for each gene. This ensures that all genes contribute equally to principal component analysis (PCA) and other methods. It is possible to regress out unwanted sources of variation, such as cell cycle effects or batch effects, during scaling, using `vars.to.regress` argument in the function. The default is `NULL`.

```{r}
seurat_object <- ScaleData(seurat_object)
```

## 5.2. Principal Component Analysis (PCA)
PCA reduces the dimensionality of the dataset by capturing the major sources of variation in a smaller number of principal components (PCs). `RunPCA()` computes the PCs, which can then be visualized with `ElbowPlot()` or `DimHeatmap()` to determine how many PCs capture meaningful biological signal. Wherever the elbow flattens out is the cutoff point, hence the name, since the PCs are no longer variable and retaining them can introduce noise in the data.

```{r}
seurat_object <- RunPCA(seurat_object)
ElbowPlot(seurat_object)
```

## 5.3. Nonlinear Embedding
Techniques such as UMAP (`RunUMAP()`) or t-SNE (`RunTSNE()`) project high-dimensional data into two or three dimensions for visualization. These plots show how cells cluster together based on their transcriptomes. tSNE are generally considered to be old-fashioned now, and UMAPs are widely used. The `dims` argument would depend on the cutoff point we get from the `ElbowPlot()`.
`DimPlot()` (dimensionality reduction plot) function helps to visual the generated UMAP.

```{r}
seurat_object <- RunUMAP(seurat_object, dims = 1:10)
DimPlot(seurat_object, reduction = "umap")
```

# 6. Integration (`02_Integration.Rmd`) 

Integration is the process of removing batch effects that exist in the data which could have been introduced by variability during processing, such as library prep, processing on different days, or samples sequenced at separate time points. A poorly integrated dataset has minimal overlap between samples which makes the data hard to compare. A well integrated dataset has cells intermixed irrespective of sample origin. This can be visualised using `DimPlot()`. There are exceptions where some clusters come only from some samples which can happen due to true biological differences between samples, like if the patient sample was significantly different from healthy controls, like is the case in extreme conditions.
There are many integration packages that exist, including Seurat's anchor based integration. There are pros and cons to every approach but Harmony is used industry wide due to being really fast and accurate.

```{r}
library(harmony)
seurat_object <- RunHarmony(seurat_object, group.by.vars = "sample_names")
```


# 7. Clustering Cells (`03_Clustering.Rmd`)
## 7.1. Constructing a Graph
`FindNeighbors()` constructs a nearest-neighbor graph based on the selected PCs. This graph represents the similarity between cells in the reduced feature space, and is useful to cluster similar cells together which is the following step.

```{r}
seurat_object <- FindNeighbors(seurat_object, dims = 1:10)
```

## 7.2. Community Detection
`FindClusters()` applies a community detection algorithm to the graph, grouping cells into clusters that ideally correspond to distinct cell types or states. The `resolution` parameter controls the granularity of clustering; higher resolution results in more clusters.

```{r}
seurat_object <- FindClusters(seurat_object, resolution = 0.5)
```

## 7.3. Visualizing Clusters
Clusters are visualized on UMAP or t-SNE plots using `DimPlot()`. Cells are colored by their cluster assignments to assess the separation and composition of clusters.

```{r}
DimPlot(seurat_object, reduction = "umap", label = TRUE)
```

# 8. Identifying Marker Genes and Differential Expression 
## 8.1. Finding Markers
`FindAllMarkers()` identifies genes that are significantly enriched in each cluster compared to others. These marker genes help interpret the biological identity of each cluster. Markers can also be found to compare conditions. This helps to identify what markers might be enriched in the condition as compared to healthy controls. Condition markers are commonly plotted on a volcano plot for pictorial comparison.

```{r}
markers <- FindAllMarkers(seurat_object)
```

## 8.2. Validating Markers
Marker expression can be visualized with `FeaturePlot()`, which overlays gene expression on UMAP or t-SNE plots, or with `VlnPlot()`, which shows expression distributions within clusters. Cross-referencing known marker genes from the literature assists in assigning cell type identities.

```{r}
FeaturePlot(seurat_object, features = c("CD3D", "MS4A1"))
VlnPlot(seurat_object, features = c("CD3D"))
```

# 9. Visualization and Interpretation
## 9.1. Plot Types
Common plot types include:

1. UMAP/t-SNE: Shows cluster structure and relationships.
2. Feature Plot: Displays gene expression spatially.
3. Violin Plot: Shows expression distributions within clusters.
4. Dot Plot: Summarizes gene expression across clusters.
5. Heatmap: Sumarizes gene expression similar to a dotplot but in a coloured tile format.

Visual inspection at each stage helps confirm whether filtering, clustering, and annotation align with biological expectations.

```{r}
DimPlot(seurat_object, reduction = "umap", label = TRUE)
FeaturePlot(seurat_object, features = c("CD3D", "MS4A1"))
VlnPlot(seurat_object, features = c("CD3D"))
DotPlot(seurat_object, features = c("CD3D", "MS4A1", "LYZ"))
Heatmap()
```

## 9.2. Iterative Refinement
Interpretation is iterative. Clusters may be split further or merged based on marker genes and biological context. Cells with ambiguous identity can be excluded. Biological hypotheses can be developed from the marker genes and cluster relationships.

# 10. Subset Analysis (`04_Subset_Clustering.Rmd`)
Subsetting the data allows to focus on the cell type(s) of interest by isolating it from the rest of the dataset. This is done using the same `subset()` function but uses the `seurat_clusters` as an identifier or a new annotated cell type column.

```{r}
glia_subset <- subset(seurat_object, idents = "2")
tcell_subset <- subset(seurat_object, idents = c("CD4 T", "CD8 T", "gDT"))
```

- The subsetted data needs to go through a lot of the same processing as the main object. QC is not required since we are subsetting the data post QC, but the data needs to be re-normalized and re-scaled. This is because both normalization and scaling are relative to the entire dataset, and once we have reduced that to a small fraction of the original, we need to adjust gene expression relative to the smaller subset.
- Integration is not required since the dataset has already been integrated, but sometimes the data still shows batch effect. This can be checked by plotting a `DimPlot()` on the samples, and if they are not inter-mixed well, then the subset might need re-integration.
- After this the subset can be reclustered using the same steps as the main object for cellular subtypes (e.g, myeloid cells to classical monocytes, inflammatory macrophages, dendritic cells, etc.)
- Condition-wise comparisons are also interesting here since the difference between condition and healthy control can be studied in detail in each subcluster.
- All the further analysis like GSEA, trajectory analysis, cellchat can be done on the subsetted data as well.

# 11. Gene Set Enrichment Analysis - GSEA (`05_GSEA_Clusterprofiler.Rmd`)
Gene Set Enrichment Analysis (GSEA) is used to identify whether specific biological pathways or gene sets are overrepresented in differentially expressed genes. This helps to interpret clusters in terms of known biological functions rather than individual genes.

- In the workflow, `clusterProfiler` is a commonly used R package for running GSEA on marker genes. Significant pathways can be visualized with bar plots or dot plots to show enrichment scores and adjusted p-values.
- Typical workflow steps include preparing a ranked gene list, performing enrichment, and visualizing results:

```r
library(clusterProfiler)
gsea_result <- gseGO(geneList = ranked_genes, OrgDb = org.Hs.eg.db, ont = "BP")
```
This ranks genes by log fold-change or a similar statistic and tests for enrichment in Gene Ontology terms or KEGG pathways.

---

Beyond the core workflow, additional analyses include:

# 12. Trajectory Inference (`06_Trajectory_Analysis.Rmd`)
Trajectory inference arranges single cells along a pseudotime axis to infer dynamic biological processes such as cell differentiation or activation. This analysis is particularly useful for studying developmental lineages or transitional cell states that are not captured by discrete clusters alone.

- Tools such as `Monocle3` or `Slingshot` perform trajectory inference by learning a graph of cells and fitting trajectories through the embedding (such as UMAP).
- A typical trajectory analysis involves creating a `cell_data_set` object, reducing dimensions, clustering, learning a graph, and ordering cells.
- The inferred pseudotime values can be visualized on the UMAP to show how cells transition from an origin state to a target state.

# 13. Cell-Cell Communication (`07_CellChat.Rmd`)

Strongly recommend following along the YouTube tutorial : [CellChat Tutorial](https://youtu.be/mOSG31aaj14?feature=shared)

Cell-cell communication analysis predicts potential ligand-receptor interactions between cell types or clusters. This provides insight into how cells may influence each other through secreted factors and signaling pathways.

- Tools such as `CellPhoneDB` and `CellChat` compare known ligand-receptor pairs with the expression data to identify enriched interactions.
- Results can be visualized using network plots or chord diagrams to highlight which cell types send or receive signals.
- In a typical CellChat workflow, the object is created, interactions are inferred, and results are visualized.
- This analysis extends the interpretation of clusters by linking them through predicted signaling interactions, revealing potential mechanisms of coordination in the tissue microenvironment.
