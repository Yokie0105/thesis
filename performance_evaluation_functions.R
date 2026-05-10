library(voomCLR)
library(limma)  
library(modeest)

setwd("C:/Users/sarah/OneDrive - KU Leuven/Desktop/Master's Thesis")

################################################################################

# VOOMCLR FUNCTION

run.voomCLR <- function(sim.data, design.formula, exposure.coef="GroupACPA+") {
  
  # Simulated data
  Y <- sim.data$Y
  metadata <- sim.data$metadata
  celltypes <- rownames(Y)
  
  # Truth
  true.beta <- sim.data$betas$beta 
  is.true.signal <- abs(true.beta) > 0  
  
  ## VoomCLR 
  
  # Analysis
  design <- model.matrix(as.formula(design.formula), data=metadata)
  
  v <- voomCLR(counts=Y, design=design, varCalc="analytical", varDistribution="NB", plot=FALSE, span=0.8)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  
  # Output
  raw.beta <- fit$coefficients[, exposure.coef]
  raw.se <- fit$stdev.unscaled[, exposure.coef] * sqrt(fit$s2.post)
  df.tot <- fit$df.total
  
  ## Mode for bias correction
  
  # Bias correction 
  tt <- topTableBC(fit, coef=exposure.coef, n=Inf)
  tt <- tt[celltypes, ] # Keep original order of celltypes
  
  # Bias-corrected beta
  mode.beta <- tt$logFC
  
  # Estimated bias
  mode.bias <- raw.beta - mode.beta
  
  ## Voting for bias correction
  
  # Bias correction
  votes <- matrix(NA, nrow=length(celltypes), ncol=length(celltypes))
  rownames(votes) <- celltypes
  colnames(votes) <- celltypes
  
  for (p in 1:length(celltypes)) {
    votes[p, ] <- raw.beta - raw.beta[p]
  }
  
  # Bias-corrected beta
  vote.beta <- apply(votes, 2, function(x) modeest::mlv(x, method="meanshift", kernel="gaussian"))
  
  # Estimated bias
  vote.bias <- raw.beta - vote.beta
  
  # P values
  vote.p <- 2 * pt(-abs(vote.beta / raw.se), df=df.tot)
  vote.adj.p <- p.adjust(vote.p, method="BH")
  
  ## Output
  
  # Final output
  results <- data.frame(celltype = celltypes,
                        true.beta = true.beta,
                        is.true.signal = is.true.signal,
                        df = df.tot,
                        se = raw.se, 
                        
                        raw.beta = raw.beta,
                        raw.error = raw.beta - true.beta,
                        
                        mode.beta = mode.beta,
                        mode.bias = mode.bias,
                        mode.error = mode.beta - true.beta,
                        mode.p = tt$P.Value,
                        mode.adj.p = tt$adj.P.Val,
                        
                        vote.beta = vote.beta,
                        vote.bias = vote.bias,
                        vote.error = vote.beta - true.beta,
                        vote.p = vote.p,
                        vote.adj.p = vote.adj.p)
  
  return(results)
  
}

################################################################################

# PERFORMANCE EVALUATION FUNCTION

evaluate.performance <- function(results, alphas=c(0.01, 0.05, 0.10)) {
  
  # True (number of) signal cell types
  true.k <- sum(results$is.true.signal)
  true.signal.celltypes <- results$celltype[results$is.true.signal]
  
  # Predicted signal cell types
  top.k.mode <- results$celltype[order(results$mode.adj.p)][1:true.k]
  top.k.vote <- results$celltype[order(results$vote.adj.p)][1:true.k]
  
  # Coverage function
  calculate.coverage <- function(est.beta, true.beta, se, df, target.cov) {
    t.crit <- qt(1 - (1 - target.cov) / 2, df)
    lower <- est.beta - t.crit * se
    upper <- est.beta + t.crit * se
    return(mean(true.beta >= lower & true.beta <= upper))
  }
  
  # Rates function
  calculate.rates <- function(adj.p, is.true.signal, alpha) {
    pred.pos <- adj.p < alpha
    
    TP <- sum(pred.pos & is.true.signal)
    FP <- sum(pred.pos & !is.true.signal)
    TN <- sum(!pred.pos & !is.true.signal)
    FN <- sum(!pred.pos & is.true.signal)
    
    TPR <- ifelse((TP + FN)==0, 0, TP / (TP + FN)) 
    FPR <- ifelse((FP + TN)==0, 0, FP / (FP + TN)) 
    FDR <- ifelse((TP + FP)==0, 0, FP / (TP + FP)) 
    
    return(c(TPR=TPR, FPR=FPR, FDR=FDR))
  }
  
  metrics <- c(
    # Coverages
    mode.cov.90 = calculate.coverage(results$mode.beta, results$true.beta, results$se, results$df, 0.90),
    vote.cov.90 = calculate.coverage(results$vote.beta, results$true.beta, results$se, results$df, 0.90),
    
    mode.cov.95 = calculate.coverage(results$mode.beta, results$true.beta, results$se, results$df, 0.95),
    vote.cov.95 = calculate.coverage(results$vote.beta, results$true.beta, results$se, results$df, 0.95),
    
    mode.cov.99 = calculate.coverage(results$mode.beta, results$true.beta, results$se, results$df, 0.99),
    vote.cov.99 = calculate.coverage(results$vote.beta, results$true.beta, results$se, results$df, 0.99),
    
    # Top k accuracy
    mode.k.accuracy = as.numeric(sum(true.signal.celltypes %in% top.k.mode) == true.k),
    vote.k.accuracy = as.numeric(sum(true.signal.celltypes %in% top.k.vote) == true.k)
  )
  
  for (a in alphas) {
    # Rates (TPR, FPR, FDR)
    mode.rates <- calculate.rates(results$mode.adj.p, results$is.true.signal, a)
    vote.rates <- calculate.rates(results$vote.adj.p, results$is.true.signal, a)
    
    metrics[paste0("mode_TPR_", a)] <- mode.rates["TPR"]
    metrics[paste0("mode_FPR_", a)] <- mode.rates["FPR"]
    metrics[paste0("mode_FDR_", a)] <- mode.rates["FDR"]
    
    metrics[paste0("vote_TPR_", a)] <- vote.rates["TPR"]
    metrics[paste0("vote_FPR_", a)] <- vote.rates["FPR"]
    metrics[paste0("vote_FDR_", a)] <- vote.rates["FDR"]
  }
  
  return(metrics)
  
}