import scanpy as sc
import numpy as np

adata = sc.read_h5ad("seurat_object_new.h5ad")

# Check metadata 
adata.obs

# SeuratDisk and sceasy (amongst other packages), commonly convert cluster/label information into numeric values so those need to be fixed

# Converting seurat clustering information into string from numeric continuous values
adata.obs['seurat_clusters'] = adata.obs['seurat_clusters'].astype(str)

# Cluster information gets lost in conversion and only the index numeric values are retained. So annotation information needs to be added back to new anndata
cluster_mapping = { 
                        "0" : "Plasma",
                        "1" : "Plasma",
                        "2" : "T cells"
                  }

# Apply mapping to seurat_clusters indexing
adata.obs['cluster_label'] = adata.obs['seurat_clusters'].map(cluster_mapping)

# Change data type of column to categorical
adata.obs["cluster_label"] = adata.obs["cluster_label"].astype("category")

# Check UMAP to see changes are in place and are correct
sc.pl.umap(adata, color = "seurat_clusters", title='')
sc.pl.umap(adata, color = "cluster_label", title='')

# Check normalised counts. The output should reflect a Compressed Sparse Row format, ideally in numpy.float32 for efficiency 
adata.raw.X
