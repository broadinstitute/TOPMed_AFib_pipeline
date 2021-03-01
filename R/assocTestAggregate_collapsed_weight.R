library(Matrix)
setGeneric("assocTestAggregate_v2", function(gdsobj, ...) standardGeneric("assocTestAggregate_v2"))


.match.arg <- function(test) {
    if (length(test) > 1) test <- NULL
    match.arg(test, choices=c("Collapse","Burden", "SKAT", "fastSKAT", "SMMAT", "fastSMMAT", "SKATO"))
}


setMethod("assocTestAggregate_v2",
          "SeqVarIterator",
          function(gdsobj, null.model, AF.max=1,
                   weight.beta=c(1,1), weight.user=NULL,
                   test=c("Collapse","Burden", "SKAT", "fastSKAT", "SMMAT", "SKATO"),
                   # burden.test=c("Score"),
                   # pval.method=c("davies", "kuonen", "liu"),
                   neig = 200, ntrace = 500,
                   rho = seq(from = 0, to = 1, by = 0.1),
                   sparse=TRUE, imputed=FALSE,
                   male.diploid=TRUE, genome.build=c("hg19", "hg38"),
                   verbose=TRUE) {

              # check argument values
              test <- .match.arg(test)
              # burden.test <- match.arg(burden.test)
              # pval.method <- match.arg(pval.method)

              # don't use sparse matrices for imputed dosages
              if (imputed) sparse <- FALSE

              # coerce null.model if necessary
              if (sparse) null.model <- GENESIS:::.nullModelAsMatrix(null.model)

              # filter samples to match null model
              sample.index <- GENESIS:::.setFilterNullModel(gdsobj, null.model, verbose=verbose)

              # do we need to match on alleles?
              match.alleles <- any(c("ref", "alt") %in% names(mcols(currentRanges(gdsobj))))

              # check ploidy
              if (SeqVarTools:::.ploidy(gdsobj) == 1) male.diploid <- FALSE

              # results
              res <- list()
              res.var <- list()
              i <- 1
              n.iter <- length(variantFilter(gdsobj))
              set.messages <- ceiling(n.iter / 100) # max messages = 100
              iterate <- TRUE
              while (iterate) {
                  var.info <- variantInfo(gdsobj, alleles=match.alleles, expanded=TRUE)

                  if (!imputed) {
                      geno <- expandedAltDosage(gdsobj, use.names=FALSE, sparse=sparse)[sample.index,,drop=FALSE]
                  } else {
                      geno <- imputedDosage(gdsobj, use.names=FALSE)[sample.index,,drop=FALSE]
                  }

                  if (match.alleles) {
                      index <- GENESIS:::.matchAlleles(gdsobj, var.info)
                      var.info <- var.info[index,,drop=FALSE]
                      geno <- geno[,index,drop=FALSE]
                  } else {
                      index <- NULL
                  }

                  # number of non-missing samples
                  # n.obs <- colSums(!is.na(geno))
                  n.obs <- GENESIS:::.countNonMissing(geno, MARGIN = 2)

                  # allele frequency
                  freq <- GENESIS:::.alleleFreq(gdsobj, geno, variant.index=index, sample.index=sample.index,
                                      male.diploid=male.diploid, genome.build=genome.build)

                  # filter monomorphic variants
                  keep <- GENESIS:::.filterMonomorphic(geno, count=n.obs, freq=freq$freq, imputed=imputed)

                  # exclude variants with freq > max
                  keep <-  keep & freq$freq <= AF.max
                  if (!all(keep)) {
                      var.info <- var.info[keep,,drop=FALSE]
                      geno <- geno[,keep,drop=FALSE]
                      n.obs <- n.obs[keep]
                      freq <- freq[keep,,drop=FALSE]
                  }

                  # weights
                  if (is.null(weight.user)) {
                      # Beta weights
                      weight <- GENESIS:::.weightFromFreq(freq$freq, weight.beta)
                  } else {
                      # user supplied weights
                      weight <- currentVariants(gdsobj)[[weight.user]][expandedVariantIndex(gdsobj)]
                      if (!is.null(index)) weight <- weight[index]
                      weight <- weight[keep]

                      weight0 <- is.na(weight) | weight == 0
                      if (any(weight0)) {
                          keep <- !weight0
                          var.info <- var.info[keep,,drop=FALSE]
                          geno <- geno[,keep,drop=FALSE]
                          n.obs <- n.obs[keep]
                          freq <- freq[keep,,drop=FALSE]
                          weight <- weight[keep]
                      }
                  }

                  # number of variant sites
                  n.site <- length(unique(var.info$variant.id))

                  # number of alternate alleles
                  n.alt <- sum(geno, na.rm=TRUE)

                  # number of samples with observed alternate alleles > 0
                  n.sample.alt <- sum(rowSums(geno, na.rm=TRUE) > 0)

                  res[[i]] <- data.frame(n.site, n.alt, n.sample.alt)
                  res.var[[i]] <- cbind(var.info, n.obs, freq, weight)

                  if (n.site > 0) {
                      # zero impute missing values
                    if (any(n.obs < nrow(geno)) & test %in% c("Collapse","Burden")) {
                        geno <- zeroImpute_Sean(geno, freq$freq)
                    }else if(any(n.obs < nrow(geno)) & !test %in% c("Collapse","Burden")) {
                      # mean impute missing values
                        geno <- GENESIS:::.meanImpute(geno, freq$freq)
                      }

                      # do the test
                      assoc <- testVariantSet(null.model, G=geno, weights=weight,
                                              test=test, # burden.test=burden.test,
                                              neig = neig, ntrace = ntrace,
                                              rho=rho)
                                              # pval.method=pval.method)
                      res[[i]] <- cbind(res[[i]], assoc, stringsAsFactors=FALSE)
                  }

                  if (verbose & n.iter > 1 & i %% set.messages == 0) {
                      message(paste("Iteration", i , "of", n.iter, "completed"))
                  }
                  i <- i + 1
                  iterate <- SeqVarTools::iterateFilter(gdsobj,verbose=FALSE)
              }

              res <- list(results=dplyr::bind_rows(res), variantInfo=res.var)
              GENESIS:::.annotateAssoc(gdsobj, res)
          })


#
zeroImpute_Sean <- function(geno, freq) {
    miss.idx <- Matrix::which(is.na(geno))
    miss.var.idx <- ceiling(miss.idx/nrow(geno))
    geno[miss.idx] <- 0 * freq[miss.var.idx]
    geno
}


testVariantSet <- function( nullmod, G, weights,
                                      test = c("Collapse","Burden", "SKAT", "fastSKAT", "SMMAT", "fastSMMAT", "SKATO"),
                                      # burden.test = c("Score"),
                                      neig = 200, ntrace = 500,
                                      rho = seq(from = 0, to = 1, by = 0.1)){
                                     # pval.method = c("davies", "kuonen", "liu"),
                                     # return.scores = FALSE, return.scores.cov = FALSE){

              test <- match.arg(test)
              # burden.test <- match.arg(burden.test)
              # pval.method <- match.arg(pval.method)

              G <- GENESIS:::.genoAsMatrix(nullmod, G)
              if (test == "Collapse") {
                  out <- .testVariantSetCollapse(nullmod, G, weights, burden.test = "Score")
              }
              if (test == "Burden") {
                  out <- GENESIS:::.testVariantSetBurden(nullmod, G, weights, burden.test = "Score")
              }
              if (test == "SKAT") {
                  out <- GENESIS:::.testVariantSetSKAT(nullmod, G, weights, neig = Inf, ntrace = Inf)
                                             # return.scores, return.scores.cov)
              }
              if(test == "fastSKAT"){
                  out <- GENESIS:::.testVariantSetSKAT(nullmod, G, weights, neig, ntrace)
              }
              if (test == "SMMAT") {
                  out <- GENESIS:::.testVariantSetSMMAT(nullmod, G, weights, neig = Inf, ntrace = Inf)
              }
              if(test == "fastSMMAT"){
                  out <- GENESIS:::.testVariantSetSMMAT(nullmod, G, weights, neig, ntrace)
              }
              if(test == "SKATO"){
                  out <- GENESIS:::.testVariantSetSKATO(nullmod, G, weights, rho)
              }
              return(out)
          }



          ## create the burden score, than calls the appropriate single variant test function.
          ## can easily implement GxE interaction with the burden score... later!
.testVariantSetCollapse <- function(nullmod, G, weights, burden.test){
              # multiply G by weights and compute burden
              if(is(G, "Matrix")){
                  burden <- rowSums(G %*% Matrix::Diagonal(x = weights))
              }else{
                  burden <- colSums(t(G) * weights)
              }
              burden<-ifelse(burden>=1,1,0)
              # adjust burden for covariates and random effects
              Gtilde <- GENESIS:::calcGtilde(nullmod, burden)
              if(is.null(nullmod$RSS0)){
                  nullmod$RSS0 <- as.numeric(crossprod(nullmod$Ytilde))
              }

              if (burden.test == "Score") {
                  out <- GENESIS:::.testGenoSingleVarScore(Gtilde, G = burden, resid = nullmod$resid, RSS0 = nullmod$RSS0)
              }
              # if (burden.test == "Wald"){
              #     out <- .testGenoSingleVarWald(Gtilde, Ytilde = nullmod$Ytilde,
              #                                   n = length(nullmod$Ytilde), k = ncol(nullmod$model.matrix))
              # }
              return(out)
          }

####
##### function start
aggreated_test<-function(num,gdsfile,varfile,groupfile,phenfile,nullfile,cutoff,stat,outfile){

##### samples
phen1<-fread(phenfile,header=T,data.table=F,sep="\t")
names(phen1)[1]<-"sample.id"
samid0<-phen1$sample.id

######
###### QCed variants
vardata<-fread(varfile,header=F,sep="\t",data.table=F)
varid0<-vardata$V1

######
###### read gds file
gds <- seqOpen(gdsfile, allow.duplicate=T)
samples <- seqGetData(gds, "sample.id")
missamples<-samples[!samples %in% samid0]
misphen<-data.frame(matrix(NA,nrow=length(missamples),ncol=ncol(phen1)))
colnames(misphen)<-names(phen1)
misphen$sample.id<-missamples
combphen<-rbind(phen1,misphen)
rownames(combphen)<-combphen$sample.id
combphen2<-combphen[samples,]

######
###### # construct a SeqVarData object
seqData <- SeqVarData(gds, sampleData=AnnotatedDataFrame(combphen2))

######
###### filter the gdsfile
seqSetFilter(seqData, sample.id=samid0, variant.id=varid0)

######
###### annotation file
annot<-get(load(groupfile))
annot$pos<-as.numeric(annot$pos)
######
###### grouping file
gr<-aggregateGRangesList(annot)

######
###### create the iterator
iterator <- SeqVarListIterator(seqData, variantRanges=gr)

######
###### load null model
nullmod<-get(load(nullfile))

###### perfrom assocation test
assoc <- assocTestAggregate_v2(iterator, nullmod, AF.max=cutoff, test=stat, weight.user="weight",verbose=TRUE)

#####
save(assoc,file=outfile)
#####
seqClose(gds)
}
