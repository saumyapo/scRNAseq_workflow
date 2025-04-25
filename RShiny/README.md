Directory contains scripts for RShiny applications that can be used to visualise single-cell RNAseq data or cellchat data.

## Visualising any Seurat object interactively using `ShinyCell`:
* As long as all the packages are installed, the generate RShiny App is plug and play.
* Only needs the Seurat object RDS file to generate the entire RShiny application. Can run function to visualise just one RDS per application, or also input multiple RDS files to be visualised in a single RShiny application.
* The output of running the `ShinyCell` code is a directory called `shinyApp` which contains `server.R` and `ui.R` as well as the RDS file split into appropriate segments for easy visualisations.
* Can directly `Run App` through either `server.R` or `ui.R` and do not need to make any modifications.

