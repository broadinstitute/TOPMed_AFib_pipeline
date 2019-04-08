if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager",repos="https://cloud.r-project.org")

BiocManager::install(c('SeqArray'), dependencies=TRUE, clean=TRUE, ask=FALSE, INSTALL_opts='--no-docs --no-demo --byte-compile',version="3.8")

sessionInfo()

quit("no")
