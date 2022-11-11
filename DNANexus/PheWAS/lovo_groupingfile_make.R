#!/usr/bin/env Rscript

#### binary traits
args=(commandArgs(TRUE))
group1=as.character(args[1])
group2=as.character(args[2])
gene=as.character(args[3])

library(data.table)
# Group file 1 (LOF and LOF+missense)
g1 <- get(load(group1))
g1$varid <- paste0(g1$chr, ":", g1$pos, ":", g1$ref, ":", g1$alt)
g1 <- g1[which(grepl(gene, g1$group_id)), ]

# Group file 2 (ultra-rare missense)
g2 <- get(load(group2))
g2$varid <- paste0(g2$chr, ":", g2$pos, ":", g2$ref, ":", g2$alt)
g2 <- g2[which(grepl(gene, g2$group_id)), ]

# Find unique variants across annotations
varz <- unique(c(g1$varid, g2$varid))

# Make LOVO groupings for the gene
lovo_g1 <- g1
lovo_g2 <- g2
for(var in varz){
    
    inter1 <- g1
    if(var %in% inter1$varid){
        inter1 <- inter1[-(which(inter1$varid==var)), ]
    }
    inter1$group_id <- paste0(inter1$group_id, ":LOVO", var)
    lovo_g1 <- rbind(lovo_g1, inter1)

    inter2 <- g2
    if(var %in% inter2$varid){
        inter2 <- inter2[-(which(inter2$varid==var)), ]
    }
    inter2$group_id <- paste0(inter2$group_id, ":LOVO", var)
    lovo_g2 <- rbind(lovo_g2, inter2)

}

save(lovo_g1, file='LOVO_group1.RData')
save(lovo_g2, file='LOVO_group2.RData')
