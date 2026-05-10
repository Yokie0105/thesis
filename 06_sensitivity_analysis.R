library(dplyr)
library(tidyr)
library(ggplot2)

source("parametric_simulation_function.R")
source("performance_evaluation_functions.R")

################################################################################

# DATA SIMULATION

set.seed(123)

iters <- 250
sparsity.levels <- c(0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00)

# Toggle TRUE/FALSE whether data simulation needs to be run again
run.simulation.sparsity <- TRUE

if (run.simulation.sparsity) {
  
  sim.results.sparsity <- list()
  cat("Starting data simulation...\n")
  
  for (level in as.character(sparsity.levels)) {
    
    cat(paste("Generating", iters, "datasets with", level, "proportion of signal cell types...\n"))
    sim.results.sparsity[[level]] <- list()
    
    for (i in 1:iters) {
      # Simulation function
      sim.data <- simulate.RA.counts(prop.signal=as.numeric(level))
      
      # Store simulation results
      sim.results.sparsity[[level]][[i]] <- sim.data
    }
  }
  
  cat("Data simulation complete!\n")
  
  # Save as RDS
  saveRDS(sim.results.sparsity, file="simulated_RA_data_sparsity.rds")
  cat("Data saved to simulated_RA_data_sparsity.rds\n")
  
} else {
  
  # Load previously saved RDS file
  cat("Loading previously simulated data testing sparsity assumption...\n")
  sim.results.sparsity <- readRDS("simulated_RA_data_sparsity.rds")
  cat("Data loaded successfully!\n")
  
}

################################################################################

# SENSITIVITY ANALYSIS

## Metrics

sparsity.levels <- names(sim.results.sparsity)
iters <- length(sim.results.sparsity[[1]])

all.examples <- list()
all.est.biases <- list()
all.rates <- list()
all.cov <- list()
all.k.acc <- list()

cat("Starting performance evaluation for sparsity assumption...\n")

for (level in sparsity.levels) {
  
  cat(paste("Evaluating models for", level, "proportion of signal cell types...\n"))
  
  # Get simulated data
  sim.data.level <- sim.results.sparsity[[level]]
  
  iter.voom.unadj <- list()
  iter.perf.unadj <- list()
  
  iter.voom.adj <- list()
  iter.perf.adj <- list()
  
  for (i in 1:iters) {
    
    # Model unadjusted for confounding
    voom.unadj <- run.voomCLR(sim.data.level[[i]], design.formula="~ Group + Sex + Race")
    iter.voom.unadj[[i]] <- voom.unadj
    iter.perf.unadj[[i]] <- evaluate.performance(voom.unadj)
    
    # Model adjusted for confounding
    voom.adj <- run.voomCLR(sim.data.level[[i]], design.formula="~ Group + Sex + Age + Race + Ethnicity")
    iter.voom.adj[[i]] <- voom.adj
    iter.perf.adj[[i]] <- evaluate.performance(voom.adj)
  }
  
  cat(paste("Finished evaluation for", level, "proportion of signal cell types.\n"))
  
  # 1. Output per iteration
  
  # Extract one iteration for illustration
  iter.example <- bind_rows(iter.voom.unadj[[1]] %>% mutate(Model="Unadjusted", prop.signal=as.numeric(level)),
                            iter.voom.adj[[1]] %>% mutate(Model="Adjusted", prop.signal=as.numeric(level)))
  
  all.examples[[level]] <- bind_rows(iter.example %>% transmute(celltype, Model, prop.signal, Correction="Mode", bias=mode.bias, raw.error=raw.error, corrected.error=mode.error),
                                     iter.example %>% transmute(celltype, Model, prop.signal, Correction="Vote", bias=vote.bias, raw.error=raw.error, corrected.error=vote.error))
  
  # 2. Output per level
  
  # VoomCLR results
  level.voom.unadj <- do.call(rbind, iter.voom.unadj)
  level.voom.adj   <- do.call(rbind, iter.voom.adj)
  
  # Average performance metrics
  level.avg.perf.unadj <- colMeans(do.call(rbind, iter.perf.unadj))
  level.avg.perf.adj   <- colMeans(do.call(rbind, iter.perf.adj))
  
  # Estimated biases
  level.est.bias.unadj <- level.voom.unadj %>%
    dplyr::select(mode.bias, vote.bias) %>%
    pivot_longer(cols=everything(), names_to="Correction", values_to="est.bias") %>%
    mutate(Correction = ifelse(grepl("mode", Correction), "Mode", "Vote"), Model="Unadjusted", prop.signal=as.numeric(level)) %>%
    dplyr::select(Model, Correction, prop.signal, est.bias)
  
  level.est.bias.adj <- level.voom.adj %>%
    dplyr::select(mode.bias, vote.bias) %>%
    pivot_longer(cols=everything(), names_to="Correction", values_to="est.bias") %>%
    mutate(Correction = ifelse(grepl("mode", Correction), "Mode", "Vote"), Model="Adjusted", prop.signal=as.numeric(level))%>%
    dplyr::select(Model, Correction, prop.signal, est.bias)
  
  all.est.biases[[level]] <- bind_rows(level.est.bias.unadj, level.est.bias.adj)
  
  # Performance metrics
  level.perf <- data.frame(metric=names(level.avg.perf.unadj),
                           Unadjusted=level.avg.perf.unadj,
                           Adjusted=level.avg.perf.adj) %>%
    pivot_longer(cols=-metric, names_to="Model", values_to="value") %>%
    mutate(Model = factor(Model, levels=c("Unadjusted", "Adjusted")),
           Correction = ifelse(grepl("mode", metric), "Mode", "Vote"),
           prop.signal = as.numeric(level))
  
  # Rates 
  level.rates <- level.perf %>%
    filter(grepl("_0.05", metric)) %>%
    mutate(rate.type=gsub("mode_|vote_|_.*", "", metric)) %>%
    dplyr::select(Model, Correction, prop.signal, rate.type, value) %>%
    pivot_wider(names_from=rate.type, values_from=value) %>%
    mutate(Alpha=0.05) %>%
    dplyr::select(Model, Correction, prop.signal, Alpha, TPR, FPR, FDR)
  
  all.rates[[level]] <- level.rates
  
  # Coverages
  level.cov <- level.perf %>%
    filter(grepl("cov.95", metric)) %>%
    mutate(target.cov=0.95) %>%
    dplyr::select(Model, Correction, prop.signal, target.cov, cov=value)
  
  all.cov[[level]] <- level.cov
  
  # Top k accuracy
  level.k.acc <- level.perf %>%
    filter(grepl("k.accuracy", metric)) %>%
    dplyr::select(Model, Correction, prop.signal, k.accuracy=value)
  
  all.k.acc[[level]] <- level.k.acc
  
}

cat("Performance evaluation for sparsity assumption complete!\n")

# 3. Final output

cat("Aggregating sparsity results...\n")

# Aggregated output
final.examples.sparsity <- do.call(rbind, all.examples) %>%
  mutate(Model = factor(Model, levels=c("Unadjusted", "Adjusted")),
         Correction = factor(Correction, levels=c("Mode", "Vote")))

final.est.biases.sparsity <- do.call(rbind, all.est.biases) %>%
  mutate(Model = factor(Model, levels=c("Unadjusted", "Adjusted")))

final.rates.sparsity <- do.call(rbind, all.rates) %>%
  mutate(Model = factor(Model, levels=c("Unadjusted", "Adjusted")))

final.cov.sparsity <- do.call(rbind, all.cov) %>%
  mutate(Model = factor(Model, levels=c("Unadjusted", "Adjusted")))

final.k.acc.sparsity <- do.call(rbind, all.k.acc) %>%
  mutate(Model = factor(Model, levels=c("Unadjusted", "Adjusted")))

# Aggregated mean biases
final.mean.biases.sparsity <- final.est.biases.sparsity %>%
  group_by(Model, Correction, prop.signal) %>%
  summarize(mean.bias=mean(est.bias), sd.bias=sd(est.bias), .groups="drop")

cat("Sparsity results aggregated successfully!\n")

## Visualization

cat("Visualizing sparsity results...\n")

# Estimation error (per iteration)
ggplot(final.examples.sparsity %>% 
         filter(Correction == "Mode") %>% # Filtered to omit correction overlap visually
         pivot_longer(cols=c(raw.error, corrected.error), names_to="error.calc", values_to="Error") %>%
         mutate(error.calc = factor(error.calc, levels=c("raw.error", "corrected.error"))),
       aes(x=as.factor(prop.signal), y=Error, fill=error.calc)) +
  geom_point(aes(color=error.calc), position=position_jitterdodge(jitter.width=0.1), alpha=0.2, size=1, show.legend=FALSE) +
  geom_boxplot(outlier.shape=NA) +
  geom_hline(yintercept=0, linetype="dashed", color="black") +
  facet_wrap(~Model, labeller=labeller(Model = c("Unadjusted"="Unadjusted model", "Adjusted"="Adjusted model"))) +
  labs(x="Proportion of signal cell types", y="Estimation error on logFC", fill="Error calculation") +
  scale_fill_manual(values=c("raw.error"="#F8766D", "corrected.error"="#00BFC4"), 
                    labels=c("Raw coefficients", "Corrected coefficients"),
                    breaks=c("raw.error", "corrected.error")) +
  coord_cartesian(ylim=c(-0.6, 0.6)) +
  theme_bw(base_size=16)

# Bias (per level)
ggplot(final.est.biases.sparsity %>% filter(Correction == "Mode"), aes(x=as.factor(prop.signal), y=est.bias, fill=Model)) +
  geom_violin(position=position_dodge(0.5), color=NA, alpha=0.6) +
  geom_boxplot(width=0.2, position=position_dodge(0.5), alpha=0.8, outlier.size=0.5) +
  geom_hline(yintercept=0, linetype="dashed", color="black") +
  labs(x="Proportion of signal cell types", y="Estimated bias") +
  scale_fill_manual(values=c("Unadjusted"="#F8766D", "Adjusted"="#00BFC4")) +
  coord_cartesian(ylim=c(-0.6, 0.6)) +
  theme_bw(base_size=16)

# FDR
ggplot(final.rates.sparsity, aes(x=as.factor(prop.signal), y=FDR, color=Model, linetype=Correction)) +
  geom_line(aes(group=interaction(Model, Correction)), linewidth=0.6) +
  geom_point(size=2.5) +
  geom_hline(yintercept=0.05, color="black", linetype="dotted") +
  labs(x="Proportion of signal cell types", y="False Discovery Rate (FDR)") +
  coord_cartesian(ylim=c(0, 1)) +
  theme_bw(base_size=16)

# TPR 
ggplot(final.rates.sparsity, aes(x=as.factor(prop.signal), y=TPR, color=Model, linetype=Correction)) +
  geom_line(aes(group=interaction(Model, Correction)), linewidth=0.6) +
  geom_point(size=2.5) +
  labs(x="Proportion of signal cell types", y="True Positive Rate (TPR)") +
  coord_cartesian(ylim=c(0, 1)) +
  theme_bw(base_size=16)

# Coverage 
ggplot(final.cov.sparsity, aes(x=as.factor(prop.signal), y=cov, color=Model, linetype=Correction)) +
  geom_line(aes(group=interaction(Model, Correction)), linewidth=0.6) +
  geom_point(size=2.5) +
  geom_hline(yintercept=0.95, color="black", linetype="dotted") + 
  labs(x="Proportion of signal cell types", y="Coverage percentage") +
  coord_cartesian(ylim=c(0, 1)) +
  theme_bw(base_size=16)

# Top k accuracy
ggplot(final.k.acc.sparsity, aes(x=as.factor(prop.signal), y=k.accuracy, color=Model, linetype=Correction)) +
  geom_line(aes(group=interaction(Model, Correction)), linewidth=0.6) +
  geom_point(size=2.5) +
  labs(x="Proportion of signal cell types", y="Top-k accuracy") +
  coord_cartesian(ylim=c(0, 1)) +
  theme_bw(base_size=16)

cat("Sparsity visualization complete!...\n")