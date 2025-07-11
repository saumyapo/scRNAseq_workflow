[![Analysis Tools:](https://img.shields.io/badge/Analysis%20Tools:-Seurat%20v5&#46;1&#46;0-orange.svg)](https://satijalab.org/seurat/)
[![Integration Tools:](https://img.shields.io/badge/Integration%20Tools:-Harmony%20v1&#46;2&#46;3-FFD700)](https://github.com/immunogenomics/harmony)
[![Pathway Analysis Tools:](https://img.shields.io/badge/Package%20Analysis%20Tools:-grey)](https://github.com/saumyapo/scRNAseq_workflow/) [![clusterProfiler](https://img.shields.io/badge/clusterProfiler%20v4&#46;12&#46;6-darkred.svg)](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html) [![fgsea](https://img.shields.io/badge/fgsea%20v1&#46;24&#46;0-darkgreen.svg)](https://github.com/alserglab/fgsea)

[![Cell-Cell Communication:](https://img.shields.io/badge/Cell%E2%80%94Cell%20Communication:-grey)](https://github.com/saumyapo/scRNAseq_workflow/) [![Liana](https://img.shields.io/badge/Liana%20v0&#46;1&#46;14-blue.svg)](https://saezlab.github.io/liana/index.html) [![CellChat](https://img.shields.io/badge/CellChat%20v2&#46;1&#46;2-purple.svg)](https://github.com/sqjin/CellChat)
[![RShiny:](https://img.shields.io/badge/RShiny:-ShinyCell%20v2&#46;1&#46;0-pink)](https://github.com/SGDDNB/ShinyCell)

# Single cell RNA-seq Workflow

Highly recommend: [Seurat Guided Tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial#assigning-cell-type-identity-to-clusters)

Great bioinformatics tutorial YouTube channel: [Collection of Online Tutorials](https://www.youtube.com/@Collection_of_online_tutorials/featured)

## Folder breakdown
1. `analysis` contains scripts for Seurat V5 scRNA-seq workflow from loading the files in to annotation, GSEA and optional analysis like trajectory analysis and cell-cell communication.
2. `RShiny` contains RShiny code that allows the user to interactively view Seurat or CellChat data.
3. `scripts` contains miscellaneous scripts related to analysis, like Seurat Object to AnnData conversion, or making some plots of interest.


## Citations
1. Seurat v5: <br>
Hao, Y., Stuart, T., Kowalski, M.H. et al. Dictionary learning for integrative, multimodal and scalable single-cell analysis. Nat Biotechnol 42, 293–304 (2024). https://doi.org/10.1038/s41587-023-01767-y
2. Harmony: <br>
Korsunsky, I., Millard, N., Fan, J. et al. Fast, sensitive and accurate integration of single-cell data with Harmony. Nat Methods 16, 1289–1296 (2019). https://doi.org/10.1038/s41592-019-0619-0
3. clusterProfiler: <br>
Yu G (2024). “Thirteen years of clusterProfiler.” The Innovation, 5(6), 100722. doi:10.1016/j.xinn.2024.100722, https://doi.org/10.1016/j.xinn.2024.100722.
4. fgsea: <br>
Korotkevich G, Sukhov V, Sergushichev A (2019). “Fast gene set enrichment analysis.” bioRxiv. doi:10.1101/060012, http://biorxiv.org/content/early/2016/06/20/060012.
5. CellChat: <br>
Jin, S., Guerrero-Juarez, C.F., Zhang, L. et al. Inference and analysis of cell-cell communication using CellChat. Nat Commun 12, 1088 (2021). https://doi.org/10.1038/s41467-021-21246-9
6. Liana: <br>
Dimitrov, D., Schäfer, P.S.L., Farr, E. et al. LIANA+ provides an all-in-one framework for cell–cell communication inference. Nat Cell Biol 26, 1613–1622 (2024). https://doi.org/10.1038/s41556-024-01469-w
7. ShinyCell: <br>
John F Ouyang, Uma S Kamaraj, Elaine Y Cao, Owen J L Rackham, ShinyCell: simple and sharable visualization of single-cell gene expression data, Bioinformatics, Volume 37, Issue 19, October 2021, Pages 3374–3376, https://doi.org/10.1093/bioinformatics/btab209
8. CellChat Shiny: <br>
https://github.com/sqjin/CellChatShiny
