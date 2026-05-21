library(dplyr)
library(tidyr)
library(ggplot2)
library(iCOBRA)

source("parametric_simulation_function.R")
source("performance_evaluation_functions.R")

################################################################################

# DATA SIMULATION

set.seed(123)

iters <- 1000
conf.levels <- c("No", "Moderate", "Strong")

# Toggle TRUE/FALSE whether data simulation needs to be run again
run.simulation <- TRUE

if (run.simulation) {
  
  sim.results <- list()
  cat("Starting data simulation...\n")
  
  for (level in conf.levels) {
    
    cat(paste("Generating", iters, "datasets with", level, "confounding for Age and Ethnicity...\n"))
    sim.results[[level]] <- list()
    
    for (i in 1:iters) {
      # Simulation function
      sim.data <- simulate.RA.counts(conf.level=level)
      
      # Store simulation results
      sim.results[[level]][[i]] <- sim.data
    }
  }
  
  cat("Data simulation complete!\n")
  
  # Save as RDS
  saveRDS(sim.results, file="simulated_RA_data.rds")
  cat("Data saved to simulated_RA_data.rds\n")
  
} else {
  
  # Load previously saved RDS file
  cat("Loading previously simulated data...\n")
  sim.results <- readRDS("simulated_RA_data.rds")
  cat("Data loaded successfully!\n")
  
}

################################################################################

# SIMULATION STUDY

## Metrics

conf.levels <- names(sim.results)
iters <- length(sim.results[["Strong"]])

all.examples <- list()
all.est.biases <- list()
all.rates <- list()
all.cov <- list()
all.ROC <- list()
all.k.acc <- list()

cat("Starting performance evaluation...\n")

for (level in conf.levels) {
  
  cat(paste("Evaluating models for", level, "confounding...\n"))
  
  # Get simulated data
  sim.data.level <- sim.results[[level]]
  
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
  
  cat(paste("Finished evaluation for", level, "confounding.\n"))
  
  # 1. Output per iteration
  
  # Extract one iteration for illustration
  iter.example <- bind_rows(iter.voom.unadj[[1]] %>% mutate(Model="Unadjusted", Confounding=level),
                            iter.voom.adj[[1]] %>% mutate(Model="Adjusted", Confounding=level))
  
  all.examples[[level]] <- bind_rows(iter.example %>% transmute(celltype, Model, Confounding, Correction="Mode", bias=mode.bias, raw.error=raw.error, corrected.error=mode.error),
                                     iter.example %>% transmute(celltype, Model, Confounding, Correction="Vote", bias=vote.bias, raw.error=raw.error, corrected.error=vote.error))
  
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
    mutate(Correction = ifelse(grepl("mode", Correction), "Mode", "Vote"), Model="Unadjusted", Confounding=level, est.bias) %>%
    dplyr::select(Model, Correction, Confounding, est.bias)
  
  level.est.bias.adj <- level.voom.adj %>%
    dplyr::select(mode.bias, vote.bias) %>%
    pivot_longer(cols=everything(), names_to="Correction", values_to="est.bias") %>%
    mutate(Correction = ifelse(grepl("mode", Correction), "Mode", "Vote"), Model="Adjusted", Confounding=level, est.bias)%>%
    dplyr::select(Model, Correction, Confounding, est.bias)
  
  all.est.biases[[level]] <- bind_rows(level.est.bias.unadj, level.est.bias.adj)
  
  # Performance metrics
  level.perf <- data.frame(metric=names(level.avg.perf.unadj),
                           Unadjusted=level.avg.perf.unadj,
                           Adjusted=level.avg.perf.adj) %>%
    pivot_longer(cols=-metric, names_to="Model", values_to="value") %>%
    mutate(Model = factor(Model, levels=c("Unadjusted", "Adjusted")),
           Correction = ifelse(grepl("mode", metric), "Mode", "Vote"),
           Confounding = level)
  
  # Rates 
  level.rates <- level.perf %>%
    filter(grepl("_0.", metric)) %>%
    mutate(Alpha=as.numeric(gsub(".*_", "", metric)),
           rate.type=gsub("mode_|vote_|_.*", "", metric)) %>%
    dplyr::select(-metric) %>%
    pivot_wider(names_from=rate.type, values_from=value) %>%
    dplyr::select(Model, Correction, Confounding, Alpha, TPR, FPR, FDR)
  
  all.rates[[level]] <- level.rates
  
  # Coverages
  level.cov <- level.perf %>%
    filter(grepl("cov", metric)) %>%
    mutate(target.cov=as.numeric(gsub(".*cov\\.", "", metric)) / 100) %>%
    dplyr::select(Model, Correction, Confounding, target.cov, cov=value)
  
  all.cov[[level]] <- level.cov
  
  # Cobra data
  level.pval <- data.frame(unadj.mode=level.voom.unadj$mode.p,
                           unadj.vote=level.voom.unadj$vote.p,
                           adj.mode=level.voom.adj$mode.p,
                           adj.vote=level.voom.adj$vote.p)
  level.pval[is.na(level.pval)] <- 1
  
  level.padj <- data.frame(unadj.mode=level.voom.unadj$mode.adj.p,
                           unadj.vote=level.voom.unadj$vote.adj.p,
                           adj.mode=level.voom.adj$mode.adj.p,
                           adj.vote=level.voom.adj$vote.adj.p)
  level.padj[is.na(level.padj)] <- 1
  
  level.score <- 1 - level.pval
  level.truth <- data.frame(status=as.numeric(level.voom.adj$is.true.signal))
  
  rownames(level.pval) <- paste0("sim_", 1:nrow(level.pval))
  rownames(level.padj) <- rownames(level.pval)
  rownames(level.score) <- rownames(level.pval)
  rownames(level.truth) <- rownames(level.pval)
  
  cobra.data <- COBRAData(pval=level.pval, padj=level.padj, score=level.score, truth=level.truth)
  cobra.perf <- calculate_performance(cobra.data, binary_truth="status", aspects=c("roc", "fdrtpr", "fdrtprcurve"))
  cobra.plot <- prepare_data_for_plot(cobra.perf, facetted=FALSE)
  
  level.ROC <- cobra.plot@roc
  level.ROC$Confounding <- level
  
  all.ROC[[level]] <- level.ROC
  
  # Top k accuracy
  level.k.acc <- level.perf %>%
    filter(grepl("k.accuracy", metric)) %>%
    dplyr::select(Model, Correction, Confounding, k.accuracy=value)
  
  all.k.acc[[level]] <- level.k.acc
  
}

cat("Performance evaluation complete!\n")

# 3. Final output

cat("Aggregating results...\n")

# Aggregated examples
final.examples <- do.call(rbind, all.examples) %>%
  mutate(Confounding = factor(Confounding, levels=c("No", "Moderate", "Strong")),
         Model = factor(Model, levels=c("Unadjusted", "Adjusted")),
         Correction = factor(Correction, levels=c("Mode", "Vote")))

final.examples.long <- final.examples %>%
  pivot_longer(cols=c(raw.error, corrected.error), names_to="error.calc", values_to="Error")

final.examples.mean <- final.examples.long %>%
  group_by(Model, Confounding, Correction, error.calc) %>%
  summarize(mean.error = mean(Error), .groups="drop")

# Aggregated estimated and mean biases
final.est.biases <- do.call(rbind, all.est.biases) %>%
  mutate(Confounding = factor(Confounding, levels=c("No", "Moderate", "Strong")),
         Model = factor(Model, levels=c("Unadjusted", "Adjusted")))

final.mean.biases <- final.est.biases %>%
  group_by(Model, Correction, Confounding) %>%
  summarize(mean.bias=mean(est.bias), sd.bias=sd(est.bias),.groups = "drop")

# Aggregated rates
final.rates <- do.call(rbind, all.rates) %>%
  mutate(Confounding = factor(Confounding, levels=c("No", "Moderate", "Strong")),
         Model = factor(Model, levels=c("Unadjusted", "Adjusted")))

# Aggregated coverages
final.cov <- do.call(rbind, all.cov) %>%
  mutate(Confounding = factor(Confounding, levels=c("No", "Moderate", "Strong")),
         Model = factor(Model, levels=c("Unadjusted", "Adjusted")))

# Aggregated ROC
final.ROC <- do.call(rbind, all.ROC) %>%
  mutate(Confounding = factor(Confounding, levels=c("No", "Moderate", "Strong")),
         Model = factor(ifelse(grepl("unadj", method), "Unadjusted", "Adjusted"), levels=c("Unadjusted", "Adjusted")),
         Correction = ifelse(grepl("mode", method), "Mode", "Vote"))

# Aggregated top k accuracies
final.k.acc <- do.call(rbind, all.k.acc) %>%
  mutate(Confounding = factor(Confounding, levels=c("No", "Moderate", "Strong")),
         Model = factor(Model, levels=c("Unadjusted", "Adjusted")))

# Covariate examples 
covariate.examples <- bind_rows(
  sim.results[["No"]][[1]]$metadata %>% mutate(Confounding="No"),
  sim.results[["Moderate"]][[1]]$metadata %>% mutate(Confounding="Moderate"),
  sim.results[["Strong"]][[1]]$metadata %>% mutate(Confounding="Strong")) %>%
  mutate(Confounding = factor(Confounding, levels=c("No", "Moderate", "Strong")),
         Ethnicity = factor(Ethnicity, levels=c("Non-Hispanic", "Hispanic")))

mean.Age.by.Exposure <- covariate.examples %>%
  group_by(Confounding, Group) %>%
  summarize(mean.Age = mean(Age, na.rm=TRUE), .groups="drop")

cat("Results aggregated successfully!\n")

## Visualization

cat("Visualizing results...\n")

# FDR
ggplot(final.rates, aes(x=Alpha, y=FDR, color=Confounding, shape=Model, linetype=Correction)) +
  geom_line(linewidth=0.6) +
  geom_point(size=2.5) +
  geom_abline(intercept=0, slope=1, color="black", linetype="dotted") +
  labs(x="Significance level", y="False Discovery Rate (FDR)") +
  scale_x_continuous(breaks=c(0.01, 0.05, 0.10)) +
  scale_color_manual(values=c("No"="#7CAE00", "Moderate"="#00BFC4", "Strong"="#F8766D")) +
  coord_cartesian(ylim=c(0, 1)) +
  theme_bw(base_size=16)

# TPR 
ggplot(final.rates, aes(x=Alpha, y=TPR, color=Confounding, shape=Model, linetype=Correction)) +
  geom_line(linewidth=0.6) +
  geom_point(size=2.5) +
  labs(x="Significance level", y="True Positive Rate (TPR)") +
  scale_x_continuous(breaks=c(0.01, 0.05, 0.10)) +
  scale_color_manual(values=c("No"="#7CAE00", "Moderate"="#00BFC4", "Strong"="#F8766D")) +
  coord_cartesian(ylim=c(0, 1)) +
  theme_bw(base_size=16)

# Coverage 
ggplot(final.cov, aes(x=target.cov, y=cov, color=Confounding, shape=Model, linetype=Correction)) +
  geom_line(linewidth=0.6) +
  geom_point(size=2.5) +
  geom_abline(intercept=0, slope=1, color="black", linetype="dotted") + 
  labs(x="Target coverage", y="Coverage percentage") +
  scale_x_continuous(breaks=c(0.90, 0.95, 0.99)) +
  scale_color_manual(values=c("No"="#7CAE00", "Moderate"="#00BFC4", "Strong"="#F8766D")) +
  coord_cartesian(ylim=c(0, 1)) +
  theme_bw(base_size=16)

# ROC curve
ggplot(final.ROC, aes(x=FPR, y=TPR, color=Confounding, alpha=Model, linewidth=Model, linetype=Correction)) +
  geom_path() +
  geom_abline(intercept=0, slope=1, color="black", linetype="dotted") +
  labs(x="False Positive Rate (FPR)", y="True Positive Rate (TPR)") +
  scale_color_manual(values=c("No"="#7CAE00", "Moderate"="#00BFC4", "Strong"="#F8766D")) +
  scale_alpha_manual(values=c("Adjusted"=1, "Unadjusted"=0.4)) +
  scale_linewidth_manual(values=c("Adjusted"=0.7, "Unadjusted"=0.3)) +
  coord_cartesian(xlim=c(0, 1), ylim=c(0, 1)) +
  theme_bw(base_size=16)

# Estimation error (per iteration): mode
ggplot(final.examples.long %>% filter(Correction=="Mode"),
       aes(x=as.factor(celltype), y=Error, color=error.calc)) +
  geom_hline(aes(yintercept=0, color="black"), linetype="dashed") +
  geom_point(size=2) +
  geom_hline(data=final.examples.mean %>% filter(Correction=="Mode") %>% arrange(desc(error.calc)),
             aes(yintercept=mean.error, color=error.calc), linetype="dashed") +
  facet_grid(Model ~ Confounding, labeller=labeller(Model = c("Unadjusted"="Unadjusted model", "Adjusted"="Adjusted model"))) +
  labs(x="Cell types (1-29)", y="Estimation error on logFC", color="Error calculation") +
  scale_color_manual(values=c("raw.error"="#F8766D", "corrected.error"="#00BFC4"),
                     labels=c("Raw coefficients", "Corrected coefficients"),
                     breaks=c("raw.error", "corrected.error")) +
  coord_cartesian(ylim=c(-0.5, 0.5)) +
  theme_bw(base_size=16) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) 

# Estimation error (per iteration): vote
ggplot(final.examples.long %>% filter(Correction=="Vote"),
       aes(x=as.factor(celltype), y=Error, color=error.calc)) +
  geom_hline(aes(yintercept=0, color="black"), linetype="dashed") +
  geom_point(size=2) +
  geom_hline(data=final.examples.mean %>% filter(Correction=="Vote") %>% arrange(desc(error.calc)),
             aes(yintercept=mean.error, color=error.calc), linetype="dashed") +
  facet_grid(Model ~ Confounding, labeller=labeller(Model = c("Unadjusted"="Unadjusted model", "Adjusted"="Adjusted model"))) +
  labs(x="Cell types (1-29)", y="Estimation error on logFC", color="Error calculation") +
  scale_color_manual(values=c("raw.error"="#F8766D", "corrected.error"="#00BFC4"),
                     labels=c("Raw coefficients", "Corrected coefficients"),
                     breaks=c("raw.error", "corrected.error")) +
  coord_cartesian(ylim=c(-0.5, 0.5)) +
  theme_bw(base_size=16) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())

# Bias (per level)
ggplot(final.est.biases %>% filter(Correction == "Mode"), aes(x=Confounding, y=est.bias, fill=Model)) +
  geom_violin(position=position_dodge(0.5), color=NA, alpha=0.6) +
  geom_boxplot(width=0.2, position=position_dodge(0.5), alpha=0.8, outlier.size=0.5) +
  geom_hline(yintercept=0, linetype="dashed", color="black") +
  labs(x="Confounding level", y="Estimated bias") +
  scale_fill_manual(values=c("Unadjusted"="#F8766D", "Adjusted"="#00BFC4")) +
  coord_cartesian(ylim=c(-0.75, 0.75)) +
  theme_bw(base_size=16)

# Top k accuracy
ggplot(final.k.acc, aes(x=Confounding, y=k.accuracy, shape=Model, linetype=Correction)) +
  geom_line(aes(group=interaction(Model, Correction)), linewidth=0.6) +
  geom_point(size=2.5) +
  labs(x="Confounding level", y="Top-k accuracy") +
  coord_cartesian(ylim=c(0, 1)) +
  theme_bw(base_size=16)

# Age plot
ggplot(covariate.examples, aes(x=Age)) +
  geom_density() +
  facet_grid(Confounding ~ Group) +
  geom_vline(data=mean.Age.by.Exposure, aes(xintercept=mean.Age, color=Group)) +
  scale_color_manual(values=c("ACPA+"="#00BFC4", "ACPA-"="#F8766D")) +
  labs(x="Age (years)", y="Density") +
  theme_bw(base_size=24)

# Ethnicity plot
ggplot(covariate.examples, aes(x=Group, fill=Ethnicity)) +
  geom_bar(position="fill") +
  facet_wrap(~Confounding) +
  labs(x="Exposure group", y="Proportion") +
  theme_bw(base_size=24)

cat("Visualization complete!...\n")