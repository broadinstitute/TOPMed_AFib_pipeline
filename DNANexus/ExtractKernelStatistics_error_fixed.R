
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
        rownames(single_var_out)<-var.id.name # added
        out <- list(NULL)
        out[['burden_out']] <- burden_out
        out[['single_var_out']] <- single_var_out
        out[['covariance_matrix']] <- V
        return(out)
}
