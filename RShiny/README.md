Directory contains scripts for RShiny applications that can be used to visualise single-cell RNAseq data or cellchat data.

`ShinyCell` is a great package that only needs the Seurat object RDS file to generate a RShiny application that can help to visualise data real time. The output of running the code is a directory called `shinyApp` which contains `server.R` and `ui.R` as well as the RDS file split into appropriate segments for easy visualisations. As long as all the packages are installed, the RShiny App is plug and play.
