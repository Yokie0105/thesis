library(dplyr)
library(MCMCpack)
library(ggplot2)

setwd("C:/Users/sarah/OneDrive - KU Leuven/Desktop/Master's Thesis")

################################################################################

# PARAMETRIC SIMULATION FUNCTION

simulate.RA.counts <- function(n0=31, n1=31, P=29, lambda=14074, 
                               conf.level="Moderate", prop.signal=0.25, sd.conf=0.2,
                               prop.Sex=0.2, prop.Age=0.2, prop.BMI=0.2, prop.Ethnicity=0.2, prop.Race=0.2,
                               plot.covariates=FALSE) {
  
  # Sample size and indices of exposure groups
  n <- n0 + n1
  idx0 <- 1:n0
  idx1 <- (n0+1):n
  
  # Number of signal cell types: default = sparsity
  k <- round(prop.signal * P)
  
  ## Covariate data: fixed confounders
  
  # Sex ~ binomial
  Sex <- c(sample(c("Female","Male"), size=n0, prob=c(0.871,0.129), replace=TRUE),
           sample(c("Female","Male"), size=n1, prob=c(0.645,0.355), replace=TRUE))
  
  # BMI ~ lognormal
  # Set realistic cut-offs
  p0.min.BMI <- plnorm(13, meanlog=3.284, sdlog=0.198)
  p0.max.BMI <- plnorm(42, meanlog=3.284, sdlog=0.198)
  
  p1.min.BMI <- plnorm(13, meanlog=3.293, sdlog=0.177)  
  p1.max.BMI <- plnorm(42, meanlog=3.293, sdlog=0.177)
  
  # Parametric simulation
  BMI <- c(qlnorm(runif(n0, p0.min.BMI, p0.max.BMI), meanlog=3.284, sdlog=0.198),
           qlnorm(runif(n1, p1.min.BMI, p1.max.BMI), meanlog=3.293, sdlog=0.177))
  
  # Race ~ multinomial
  Race <- rep("Caucasian", n)
  while(
    sum(Race[idx0] == "African American") < 1 | sum(Race[idx1] == "African American") < 1 |
    sum(Race[idx0] == "Asian") < 1 | sum(Race[idx1] == "Asian") < 1 |
    sum(Race[idx0] == "Other") < 1 | sum(Race[idx1] == "Other") < 1) {
    
    Race <- c(sample(c("Caucasian","African American","Asian","Other"), size=n0, prob=c(0.851, 0.032, 0.020,0.097), replace=TRUE),
              sample(c("Caucasian","African American","Asian","Other"), size=n1, prob=c(0.742,0.161,0.065,0.032), replace=TRUE))
  }                                                                                   
  
  ## Covariate data: tuneable confounders
  
  # Levels of confounding
  conf.params <- list("No" = 
                        list(shape0=4.639, shape1=4.639, scale0=67, scale1=67, eff.Age=0, 
                             p0=0.161, p1=0.161, eff.Ethnicity=0),
                      "Moderate" = 
                        list(shape0=4.154, shape1=5.370, scale0=61, scale1=68, eff.Age=1,
                             p0=0.258, p1=0.065, eff.Ethnicity=1),
                      "Strong" = 
                        list(shape0=3.800, shape1=6.200, scale0=55, scale1=70, eff.Age=2,
                             p0=0.452, p1=0.050, eff.Ethnicity=2))
  
  # Check the provided conf.level
  if (!(conf.level %in% names(conf.params))) { stop("conf.level must be 'No', 'Moderate', or 'Strong'") }
  
  # Age parameters for provided conf.level
  scale0 <- conf.params[[conf.level]]$scale0
  scale1 <- conf.params[[conf.level]]$scale1
  shape0 <- conf.params[[conf.level]]$shape0
  shape1 <- conf.params[[conf.level]]$shape1
  eff.Age <- conf.params[[conf.level]]$eff.Age
  
  # Ethnicity parameters for provided conf.level
  p0 <- conf.params[[conf.level]]$p0
  p1 <- conf.params[[conf.level]]$p1
  eff.Ethnicity <- conf.params[[conf.level]]$eff.Ethnicity
  
  # Age ~ weibull
  # Set realistic cut-offs
  p0.min.Age <- pweibull(20, shape=shape0, scale=scale0)
  p0.max.Age <- pweibull(80, shape=shape0, scale=scale0)
  
  p1.min.Age <- pweibull(20, shape=shape1, scale=scale1)
  p1.max.Age <- pweibull(80, shape=shape1, scale=scale1)
  
  # Parametric simulation
  if (conf.level == "No") {
    
    # Matched data across exposure (conf.level = "No")
    Age.0 <- round(qweibull(runif(n0, p0.min.Age, p0.max.Age), shape=shape0, scale=scale0))
    if (n0 == n1) { Age.1 <- Age.0 } else { Age.1 <- sample(Age.0, size=n1, replace=TRUE) }
    Age <- c(Age.0, Age.1)
    
  } else {
    
    # Different data across exposure (conf.level = "Moderate"/"Strong")
    Age <- c(round(qweibull(runif(n0, p0.min.Age, p0.max.Age), shape=shape0, scale=scale0)),
             round(qweibull(runif(n1, p1.min.Age, p1.max.Age), shape=shape1, scale=scale1)))
  }
  
  # Ethnicity ~ binomial 
  if (conf.level == "No") {
    
    # Matched data across exposure (conf.level = "No")
    Ethnicity.0 <- rep("Non-Hispanic", n0)
    while(sum(Ethnicity.0 == "Hispanic") < 1) {
      Ethnicity.0 <- sample(c("Non-Hispanic","Hispanic"), size=n0, prob=c(1-p0, p0), replace=TRUE)
    }
    if (n0 == n1) { Ethnicity.1 <- Ethnicity.0 } else { Ethnicity.1 <- sample(Ethnicity.0, size=n1, replace=TRUE) }
    Ethnicity <- c(Ethnicity.0, Ethnicity.1)
    
  } else {
    
    # Different data across exposure (conf.level = "Moderate"/"Strong")
    Ethnicity <- rep("Non-Hispanic", n) 
    while(sum(Ethnicity[idx0] == "Hispanic") < 1 | sum(Ethnicity[idx1] == "Hispanic") < 1) {
      Ethnicity <- c(sample(c("Non-Hispanic","Hispanic"), size=n0, prob=c(1-p0, p0), replace=TRUE),
                     sample(c("Non-Hispanic","Hispanic"), size=n1, prob=c(1-p1, p1), replace=TRUE))
    }
  }
  
  ## Count data: parameter set-up
  
  # Baseline cell counts (ni) ~ negative binomial (size=overdispersion)
  mu0 <- rep(0, P)
  while (!all(mu0!=0)) {                               
    mu0 <- rnbinom(n=P, size=1/2, mu=lambda/P)         
  }
  
  # Effect sizes
  beta <- sample(c(rep(1,k), rep(0, P-k)), replace=F) * 
    rlnorm(n=P, meanlog=-0.9, sdlog=0.2) *            
    sample(c(-1,1), size=P, replace=TRUE)   
  
  if (prop.Sex == 0) {
    beta.Sex <- rep(0, P)
  } else {
    beta.Sex <- rep(0, P)
    while (all(beta.Sex == 0)) {                             
      beta.Sex <- rbinom(n=P, size=1, prob=prop.Sex) * 
        rlnorm(n=P, meanlog=-0.8, sdlog=0.2) * 
        sample(c(-1,1), size=P, replace=TRUE)            
    }
  }
  
  if (prop.Age == 0 | eff.Age == 0) {
    beta.Age <- rep(0, P)
  } else {
    beta.Age <- rep(0, P)                                 
    while (all(beta.Age == 0)) { 
      # Baseline effect size (n=1)
      base.beta.Age <- rlnorm(n=1, meanlog=-3.5, sdlog=0.2) * 
        sample(c(-1,1), size=1, replace=TRUE) *
        # Add confounding effect for conf.level
        eff.Age 
      
      # Effect sizes per cell type (n=P)
      beta.Age <- rbinom(n=P, size=1, prob=prop.Age) * 
        base.beta.Age *
        # Add bias: default = constant across all cell types, zero
        rlnorm(n=P, meanlog=0, sdlog=sd.conf)
    }
  }
  
  if (prop.BMI == 0) {
    beta.BMI <- rep(0, P)
  } else {
    beta.BMI <- rep(0, P)
    while (all(beta.BMI == 0)) {                             
      beta.BMI <- rbinom(n=P, size=1, prob=prop.BMI) * 
        rlnorm(n=P, meanlog=-3.3, sdlog=0.2) * 
        sample(c(-1,1), size=P, replace=TRUE)            
    }
  }
  
  if (prop.Ethnicity == 0 | eff.Ethnicity == 0) {
    beta.Ethnicity <- rep(0, P)
  } else {
    beta.Ethnicity <- rep(0, P)                  
    while (all(beta.Ethnicity == 0)) {   
      # Baseline effect size (n=1)
      base.beta.Ethnicity <- rlnorm(n=1, meanlog=-0.5, sdlog=0.2) * 
        sample(c(-1,1), size=1, replace=TRUE) *
        # Add confounding effect for conf.level
        eff.Ethnicity  
      
      # Effect sizes per cell type (n=P)
      beta.Ethnicity <- rbinom(n=P, size=1, prob=prop.Ethnicity) * 
        base.beta.Ethnicity *
        # Add bias: default = constant across all cell types, zero
        rlnorm(n=P, meanlog=0, sdlog=sd.conf)
    }
  }
  
  betas.Race <- matrix(0, nrow=P, ncol=3)
  colnames(betas.Race) <- c("beta1.Race", "beta2.Race", "beta3.Race")
  if (prop.Race != 0) {
    while (all(betas.Race == 0)) {                         
      affected.cells.Race <- rbinom(n=P, size=1, prob=prop.Race)    
      for (j in 1:3) {
        betas.Race[, j] <- affected.cells.Race * rlnorm(n=P, meanlog=-0.3, sdlog=0.2) * sample(c(-1,1), size=P, replace=TRUE)          
      }
    }
  }
  
  # Total cell counts (Ni) ~ poisson
  N <- rpois(n=n, lambda=lambda)
  
  ## Count data: simulation
  
  # Multinomial model
  Y <- matrix(NA, nrow=P, ncol=n)
  for (i in 1:n) {
    beta.Race <- if(Race[i] == "Caucasian") 0 
    else betas.Race[, match(Race[i], c("African American","Asian","Other")), drop=FALSE]
    
    mu.i <- mu0 * exp(beta * (ifelse(i > n0, 1, 0)) +
                        beta.Sex * (ifelse(Sex[i]=="Male", 1, 0)) + 
                        beta.Age * Age[i] +
                        beta.BMI * BMI[i] +
                        beta.Ethnicity * (ifelse(Ethnicity[i]=="Hispanic", 1, 0)) +
                        beta.Race * 1)
    # -> Log-fold-change = log(mu1/mu0) = beta
    
    # Relative abundance probabilities (πi) = multinomial logit link
    pi.i <- mu.i/sum(mu.i)
    
    # Observed cell counts ~ multinomial
    Y[,i] <- rmultinom(n=1, size=N[i], prob=pi.i)
  }
  
  # Output
  rownames(Y) <- paste0("Celltype", 1:P)
  colnames(Y) <- paste0("Sample", 1:n)
  Group <- factor(c(rep("ACPA-", n0), rep("ACPA+", n1)))
  
  metadata <- data.frame(Group, Sex, Age, BMI, Race, Ethnicity)
  signal.celltypes <- rownames(Y)[abs(beta)>0]
  
  betas <- data.frame(beta, beta.Sex, beta.Age, beta.BMI, beta.Ethnicity)
  betas <- cbind(betas, betas.Race)
  
  ## Visualization of tuneable confounders
  
  # Dataframe for plotting
  df.plot <- data.frame(Exposure = Group, 
                        Age = Age, 
                        Ethnicity = factor(Ethnicity, levels = c("Non-Hispanic", "Hispanic")))
  
  # Mean Age by exposure
  mean.Age.by.Exposure <- df.plot %>%
    group_by(Exposure) %>%
    summarize(mean.Age = mean(Age, na.rm=TRUE), .groups="drop")
  
  # Age plot
  plot.Age <- ggplot(df.plot, aes(x=Age)) +
    geom_density() +
    facet_grid(~Exposure) +
    geom_vline(data=mean.Age.by.Exposure, aes(xintercept=mean.Age, color=Exposure)) +
    scale_color_manual(values=c("ACPA+"="#00BFC4", "ACPA-"="#F8766D")) +
    labs(x="Age (years)", y="Density") +
    theme_bw(base_size=16)
  
  # Ethnicity plot
  plot.Ethnicity <- ggplot(df.plot, aes(x=Exposure, fill=Ethnicity)) +
    geom_bar(position="fill") +
    labs(x="Exposure", y="Proportion") +
    theme_bw(base_size=16)
  
  # Only print if requested
  if (plot.covariates) {
    print(plot.Age)
    print(plot.Ethnicity)
  }
  
  return(list(Y=Y, signal.celltypes=signal.celltypes, betas=betas, metadata=metadata,
              plot.Age=plot.Age, plot.Ethnicity=plot.Ethnicity))
  
}