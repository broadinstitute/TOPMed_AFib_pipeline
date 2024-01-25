
testVariantSet_ExtractKernelStatistics_ScoresAndCovarianceMatrices_Sean <- function(nullmod, G, weights, var.info, neig = Inf, ntrace = Inf, Use.SPA=F, freq,
                                                                                       SAIGEGENEplus_collapse_threshold=1){

        # Check for use.SPA, which is not yet supported
        if(Use.SPA){
                stop("SPA not yet implemented for ExtractKernelStatistics function. Stopping.")
        }

	# Modify var.info so output is in chr:pos:ref:alt format and can be compared across studies
        var.id.name <- paste0(var.info$chr, ":", var.info$pos, ":", var.info$ref, ":", var.info$alt)
        colnames(G) <- var.id.name


        if(is(G, "Matrix")){
                burden <- rowSums(G %*% Diagonal(x = weights))
                G <- G %*% Diagonal(x = weights)
        }else{
              	burden <- colSums(t(G) * weights)
                G <- t(t(G) * weights)
        }


	if(is.null(nullmod$RSS0)){
                nullmod$RSS0 <- as.numeric(crossprod(nullmod$Ytilde))
        }

	# Calculate SKAT statistic
        U <- as.vector(crossprod(G, nullmod$resid)) # WGPY
        # SKAT test statistic
        Q <- sum(U^2)

        # adjust G for covariates and random effects
        burdentilde <- GENESIS:::calcGtilde(nullmod, burden)
        Gtilde <- GENESIS:::calcGtilde(nullmod, G) # P^{1/2}GW

        # Compute SKAT Variance
        ncolGtilde <- ncol(Gtilde)
        nrowGtilde <- nrow(Gtilde)

        if (ncolGtilde <= nrowGtilde) {
                V <- crossprod(Gtilde)
        } else {
                V <- tcrossprod(Gtilde)
        }

	burden_out <- GENESIS:::.testGenoSingleVarScore(burdentilde, G = burden, resid = nullmod$resid, RSS0 = nullmod$RSS0)
        colnames(burden_out) <- paste0("Burden_", colnames(burden_out))
        single_var_out <- GENESIS:::.testGenoSingleVarScore(Gtilde, G = G, resid = nullmod$resid, RSS0 = nullmod$RSS0)
        colnames(V)<-rownames(V)<-rownames(single_var_out)<-var.id.name # added

        out <- list(NULL)
        out[['burden_out']] <- burden_out
        out[['single_var_out']] <- single_var_out
        out[['covariance_matrix']] <- V
        return(out)
}



####
#### added spa Jan/25/2024
testVariantSet_ExtractKernelStatistics_ScoresAndCovarianceMatrices_Sean <- function(nullmod, G, weights, var.info, neig = Inf, ntrace = Inf, Use.SPA=F, freq,
                                                                                       SAIGEGENEplus_collapse_threshold=1){

       
                
        
	# Modify var.info so output is in chr:pos:ref:alt format and can be compared across studies
        var.id.name <- paste0(var.info$chr, ":", var.info$pos, ":", var.info$ref, ":", var.info$alt)
        colnames(G) <- var.id.name


        if(is(G, "Matrix")){
                burden <- rowSums(G %*% Diagonal(x = weights))
                G <- G %*% Diagonal(x = weights)
        }else{
              	burden <- colSums(t(G) * weights)
                G <- t(t(G) * weights)
        }


	if(is.null(nullmod$RSS0)){
                nullmod$RSS0 <- as.numeric(crossprod(nullmod$Ytilde))
        }

	# Calculate SKAT statistic
        U <- as.vector(crossprod(G, nullmod$resid)) # WGPY
        # SKAT test statistic
        Q <- sum(U^2)

        # adjust G for covariates and random effects
        burdentilde <- GENESIS:::calcGtilde(nullmod, burden)
        Gtilde <- GENESIS:::calcGtilde(nullmod, G) # P^{1/2}GW

        # Compute SKAT Variance
        ncolGtilde <- ncol(Gtilde)
        nrowGtilde <- nrow(Gtilde)

        if (ncolGtilde <= nrowGtilde) {
                V <- crossprod(Gtilde)
        } else {
                V <- tcrossprod(Gtilde)
        }


        # Check for use.SPA
        if(Use.SPA){
        burden_out <- GENESIS:::.testGenoSingleVarScore(burdentilde, G = burden, resid = nullmod$resid, RSS0 = nullmod$RSS0)
        colnames(burden_out) <- paste0("Burden_", colnames(burden_out))
        single_var_out <- GENESIS:::.testGenoSingleVarScore(Gtilde, G = G, resid = nullmod$resid, RSS0 = nullmod$RSS0)
        colnames(V)<-rownames(V)<-rownames(single_var_out)<-var.id.name # added

        out <- list(NULL)
        out[['burden_out']] <- burden_out
        out[['single_var_out']] <- single_var_out
        out[['covariance_matrix']] <- V

        }else{
        
	# We will start with single marker tests for each of the markers
        out <- GENESIS:::.testGenoSingleVarScore(Gtilde, G = G, resid = nullmod$resid, RSS0 = nullmod$RSS0)
        outnull <- out
        # Run SPA for each marker and the combined burden of very rare variants
        out <- SPA_pval_Sean(score.result = out, nullmod = nullmod, G = as.matrix(G), pval.thresh = 0.05)
        # Compute SPA adjusted variance
        out$SPA.Score.Variance <- (out$Score^2) / qchisq(out$SPA.pval, lower.tail=F, df=1)
        out[out$SPA.Score.Variance==0,'SPA.Score.Variance'] <- sqrt(outnull[which(out$SPA.Score.Variance==0),'Score.SE'])
        out[out$SPA.Score.Variance==0,'SPA.pval'] <- outnull[which(out$SPA.Score.Variance==0),'Score.pval']
        #out <- out[,c("Score", "SPA.Score.Variance", "SPA.pval", "Est", "Est.SE")]
        colnames(out)[c(5,6)] <- c("Raw.Est", "Raw.Est.SE")
        single_var_out <- out
        colnames(V)<-rownames(V)<-rownames(single_var_out)<-var.id.name # added

        # Compute SPA adjusted Sum of Variances (SAIGE-GENE, AJHG), we will use this for SKAT test later
        V_tilde <- single_var_out$SPA.Score.Variance
        Vsum_tilde <- sum(single_var_out$SPA.Score.Variance)
        
        # We will also compute a burden test for all markers
        out <- GENESIS:::.testGenoSingleVarScore(burdentilde, G = burden, resid = nullmod$resid, RSS0 = nullmod$RSS0)
        outnull <- out
                      
        #Run SPA for the burden
        out <- SPA_pval_Sean(score.result = out, nullmod = nullmod, G = as.matrix(burden), pval.thresh = 0.05)
        # Compute SPA adjusted variance
        out$SPA.Score.Variance <- (out$Score^2) / qchisq(out$SPA.pval, lower.tail=F, df=1)
        out[out$SPA.Score.Variance==0,'SPA.Score.Variance'] <- sqrt(outnull[which(out$SPA.Score.Variance==0),'Score.SE'])
        out[out$SPA.Score.Variance==0,'SPA.pval'] <- outnull[which(out$SPA.Score.Variance==0),'Score.pval']
        #out <- out[,c("Score", "SPA.Score.Variance", "SPA.pval", "Est", "Est.SE")]
        colnames(out)[c(5,6)] <- c("Raw.Est", "Raw.Est.SE")
        colnames(out) <- paste0("Burden_", colnames(out))
        burden_out<-out
        
        # Compute SPA adjusted Burden Variance (SAIGE-GENE, AJHG)
        #Vsum_downwardhat <- out$Burden_SPA.Score.Variance

        # Compute ratio to find more conservative variance (SAIGE-GENE, AJHG)
        #r <- Vsum_tilde / Vsum_downwardhat
        #r_tilde <- min(1, r)

        # Adjust SKAT variance (SAIGE-GENE, AJHG)
        #diag(V) <- V_tilde
        #V <- V / r_tilde

        out <- list(NULL)
        out[['burden_out']] <- burden_out
        out[['single_var_out']] <- single_var_out
        out[['covariance_matrix']] <- V
        }

        return(out)
}






testVariantSet_Sean <- function( nullmod, G, weights, freq, use.weights=F, var.info,
                            test = c("Burden", "SKAT", "fastSKAT", "SMMAT", "fastSMMAT", "SKATO", "SKAT_SAIGEGENEplus", "ExtractKernelStatistics"),
                            burden.test = c("Score","Score.SPA"), collapse = FALSE, recessive=FALSE, recessive.model = c("strict", "putative"),
                            vc.type = "regular weighted", vc.test=c("Score","Score.SPA"), SAIGEGENEplus_collapse_threshold=10,
                            neig = 200, ntrace = 500,
                            rho = seq(from = 0, to = 1, by = 0.1)){
                           # pval.method = c("davies", "kuonen", "liu"),
                           # return.scores = FALSE, return.scores.cov = FALSE){

    test <- match.arg(test)
    burden.test <- match.arg(burden.test)
    vc.type <- match.arg(vc.type)
    # pval.method <- match.arg(pval.method)

    G <- GENESIS:::.genoAsMatrix(nullmod, G)
    if (test == "Burden") {
        if(collapse){
                burden.type <- "collapsing test"
        }else if(!use.weights){
                burden.type <- "regular burden"
        }else{
              	burden.type <- "externally weighted burden"
        }
	#cat('Running Burden test type', burden.type, 'using Pval method ', burden.test, '...\n')
        out <- testVariantSetBurden_Sean(nullmod, G, weights, burden.test = burden.test, collapse = collapse, recessive = recessive)
    }
    if (test == "SKAT") {
        if(vc.test=="Score.SPA"){
                stop('SPA is not yet implemented for', test, '...\n')
        }
	#cat('Running variance component-based test type', test, 'type', vc.type, '...\n')
        out <- testVariantSetSKAT_Sean(nullmod, G, weights, neig = Inf, ntrace = Inf)
                                   # return.scores, return.scores.cov)
    }
    if(test == "fastSKAT"){
        if(vc.test=="Score.SPA"){
                stop('SPA is not yet implemented for', test, '...\n')
        }
	#cat('Running variance component-based test type', test, 'type', vc.type, '...\n')
        out <- testVariantSetSKAT_Sean(nullmod, G, weights, neig, ntrace)
    }
    if (test == "SMMAT") {
        if(vc.test=="Score.SPA"){
                stop('SPA is not yet implemented for', test, '...\n')
        }
	#cat('Running variance component-based test type', test, 'type', vc.type, '...\n')
        out <- testVariantSetSMMAT_Sean(nullmod, G, weights, neig = Inf, ntrace = Inf)
    }
    if(test == "fastSMMAT"){
        if(vc.test=="Score.SPA"){
                stop('SPA is not yet implemented for', test, '...\n')
        }
	#cat('Running variance component-based test type', test, 'type', vc.type, '...\n')
        out <- testVariantSetSMMAT_Sean(nullmod, G, weights, neig, ntrace)
    }
    if(test == "SKATO"){
        if(vc.test=="Score.SPA"){
                stop('SPA is not yet implemented for', test, '...\n')
        }
	#cat('Running variance component-based test type', test, 'type', vc.type, '...\n')
        out <- testVariantSetSKATO_Sean(nullmod, G, weights, rho)
    }
    if(test == "SKAT_SAIGEGENEplus"){
        #cat('Running variance component-based test type', test, 'type', vc.type, 'using pvalue method', vc.test, '...\n')
        Use.SPA <- F
        if(vc.test=="Score.SPA"){
                Use.SPA <- T
        }
	out <- testVariantSetSKAT_SAIGEGENEplus_Sean(nullmod, G, weights, neig = Inf, ntrace = Inf, Use.SPA=Use.SPA, freq=freq,
                                                     SAIGEGENEplus_collapse_threshold=SAIGEGENEplus_collapse_threshold)
    }
    if(test == "ExtractKernelStatistics"){
        #cat('Extracting Kernel Statistics...\n')
        Use.SPA <- F
        if(vc.test=="Score.SPA"){
        	Use.SPA <- T
        }
	# SPA not yet supported.... Will implement later
        # SAIGEGENEplus_collapse not yet implemented ... Will work on this later.
        out <- testVariantSet_ExtractKernelStatistics_ScoresAndCovarianceMatrices_Sean(nullmod, G, weights, var.info, neig = Inf, ntrace = Inf, Use.SPA=Use.SPA, freq=freq,
                                                                                       SAIGEGENEplus_collapse_threshold=SAIGEGENEplus_collapse_threshold)
    }
    return(out)
}
