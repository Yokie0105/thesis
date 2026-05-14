library(voomCLR)
library(dplyr)
library(tibble)
library(ggplot2)
library(modeest)

# Load original data
ori.data <- readRDS("original_RA_data.rds")
ori.data <- na.omit(ori.data) 
head(ori.data)
summary(ori.data)

# Load imbalanced data
imb.data <- readRDS("imbalanced_RA_data.rds")
imb.data <- na.omit(imb.data) 
head(imb.data)
summary(imb.data)

# Cell types
celltypes <- c("ASDC", "CD14_monocyte", "CD16_monocyte", "CD56bright_NK", 
               "CD56dim_NK", "CD8aa", "DN_T_cell", "Effector_B_cell", 
               "Erythrocyte", "ILC", "Intermediate_monocyte", "MAIT", 
               "Memory_B_cell", "Memory_CD4_T_cell", "Memory_CD8_T_cell", 
               "Naive_B_cell", "Naive_CD4_T_cell", "Naive_CD8_T_cell", 
               "Plasma_cell", "Platelet", "Progenitor_cell", 
               "Proliferating_NK_cell", "Proliferating_T_cell", 
               "Transitional_B_cell", "Treg", "cDC1", "cDC2", "gdT", "pDC")

################################################################################

# CELL COMPOSITION ANALYSIS: ORIGINAL DATA

## Mode for bias correction

# Observed counts 
Y <- t(as.matrix(ori.data[, celltypes]))

# Relative abundances 
ori.data.RA <- ori.data
ori.data.RA[, celltypes] <- t(microbiome::transform(Y, "compositional"))

# Mean relative abundances
df.mean.RA <- data.frame(Celltype=names(colMeans(ori.data.RA[, celltypes])), mean.RA=colMeans(ori.data.RA[, celltypes]))

# Model analysis
# Adjustment: Race + Ethnicity
design.adj <- model.matrix(~ Exposure + Sex + Age + Race + Ethnicity, data=ori.data)

v.adj <- voomCLR(counts=Y, design=design.adj, varCalc="analytical", varDistribution="NB", plot=FALSE, span=0.8)
fit.adj <- lmFit(v.adj, design.adj)
fit.adj <- eBayes(fit.adj)

# Adjusted model results
fit.adj$coefficients
fit.adj$p.value
fit.adj$df.residual
fit.adj$df.prior

# Exposure outcome
tt.adj.exposure <- topTableBC(fit.adj, coef="ExposureACPA+", n=Inf)
tt.adj.exposure 

# Arrange results
results.adj.mode <- tt.adj.exposure %>%
  rownames_to_column(var="Celltype") %>%
  left_join(df.mean.RA, by="Celltype") %>%
  mutate(Celltype = gsub("_", " ", Celltype),
         Celltype = ifelse(Celltype %in% c("CD56bright NK", "CD56dim NK"), 
                           paste0(Celltype, " cell"), Celltype)) %>%
  arrange(adj.P.Val) %>% 
  mutate(Celltype = factor(Celltype, levels = rev(unique(Celltype))),
         Category = case_when(grepl("B cell|Plasma", Celltype) ~ "B cell",
                              grepl("T cell|CD8|CD4|Treg|gdT|MAIT", Celltype) ~ "T cell",
                              grepl("NK", Celltype) ~ "NK cell",
                              grepl("monocyte", Celltype) ~ "Monocyte",
                              grepl("DC|ASDC", Celltype) ~ "Dendritic cell",
                              TRUE ~ "Other"),
         Significance = ifelse(adj.P.Val < 0.05, "Significant", "Not significant"),
         Asterix = ifelse(Significance=="Significant", "*", ""))

# Visualize results
ggplot(results.adj.mode, aes(x=logFC, y=Celltype)) +
  geom_vline(xintercept=0, linetype="dashed") +
  geom_segment(aes(x=0, xend=logFC, y=Celltype, yend=Celltype, color=Category)) +
  geom_point(aes(color=Category, size=mean.RA)) +
  geom_text(aes(label=Asterix, hjust=ifelse(logFC >= 0, -1.0, 1.5)), vjust=0.75, size=7, show.legend=FALSE) +
  guides(size="none", color=guide_legend(nrow=1)) +
  scale_x_continuous(limits=c(-0.75, 0.75), breaks=seq(-0.75, 0.75, by=0.25)) +
  labs(x="LogFC", y="Cell type") +
  theme_bw(base_size=16) +
  theme(legend.position="bottom")

## Voting for bias correction

# Previous adjusted model results
beta.adj <- fit.adj$coefficients[, "ExposureACPA+"]
se.adj <- fit.adj$stdev.unscaled[, "ExposureACPA+"] * sqrt(fit.adj$s2.post)

# Voting procedure
votes.adj <- matrix(NA, nrow=length(celltypes), ncol=length(celltypes))
colnames(votes.adj) <- names(beta.adj)
rownames(votes.adj) <- names(beta.adj)

for (p in 1:length(celltypes)) {
  votes.adj[p, ] <- beta.adj - beta.adj[p]
}

vote.mode.beta.adj <- apply(votes.adj, 2, function(x) mlv(x, method="meanshift", kernel="gaussian"))
names(vote.mode.beta.adj) <- names(beta.adj)

# Ensure similar output as previously
vote.mode.t.adj <- vote.mode.beta.adj / se.adj
vote.mode.pval.adj <- 2 * pt(-abs(vote.mode.t.adj), df=fit.adj$df.total)
vote.mode.adj.pval.adj <- p.adjust(vote.mode.pval.adj, method="BH")

tt.adj.vote.exposure <- tt.adj.exposure
tt.adj.vote.exposure$logFC <- vote.mode.beta.adj[rownames(tt.adj.exposure)]
tt.adj.vote.exposure$t <- vote.mode.t.adj[rownames(tt.adj.exposure)]
tt.adj.vote.exposure$P.Value <- vote.mode.pval.adj[rownames(tt.adj.exposure)]
tt.adj.vote.exposure$adj.P.Val <- vote.mode.adj.pval.adj[rownames(tt.adj.exposure)]

# Arrange results
results.adj.vote <- tt.adj.vote.exposure %>%
  rownames_to_column(var="Celltype") %>%
  left_join(df.mean.RA, by="Celltype") %>%
  mutate(Celltype = gsub("_", " ", Celltype),
         Celltype = ifelse(Celltype %in% c("CD56bright NK", "CD56dim NK"), 
                           paste0(Celltype, " cell"), Celltype)) %>%
  arrange(adj.P.Val) %>% 
  mutate(Celltype = factor(Celltype, levels = rev(unique(Celltype))),
         Category = case_when(grepl("B cell|Plasma", Celltype) ~ "B cell",
                              grepl("T cell|CD8|CD4|Treg|gdT|MAIT", Celltype) ~ "T cell",
                              grepl("NK", Celltype) ~ "NK cell",
                              grepl("monocyte", Celltype) ~ "Monocyte",
                              grepl("DC|ASDC", Celltype) ~ "Dendritic cell",
                              TRUE ~ "Other"),
         Significance = ifelse(adj.P.Val < 0.05, "Significant", "Not significant"),
         Asterix = ifelse(Significance=="Significant", "*", ""))

# Visualize results
ggplot(results.adj.vote, aes(x=logFC, y=Celltype)) +
  geom_vline(xintercept=0, linetype="dashed") +
  geom_segment(aes(x=0, xend=logFC, y=Celltype, yend=Celltype, color=Category)) +
  geom_point(aes(color=Category, size=mean.RA)) +
  geom_text(aes(label=Asterix, hjust=ifelse(logFC >= 0, -1.0, 1.5)), vjust=0.75, size=7, show.legend=FALSE) +
  guides(size="none", color=guide_legend(nrow=1)) +
  scale_x_continuous(limits=c(-0.75, 0.75), breaks=seq(-0.75, 0.75, by=0.25)) +
  labs(x="LogFC", y="Cell type") +
  theme_bw(base_size=16) +
  theme(legend.position="bottom")

################################################################################

# CELL COMPOSITION ANALYSIS: IMBALANCED DATA

## Mode for bias correction

# Observed counts 
Y <- t(as.matrix(imb.data[, celltypes]))

# Relative abundances 
imb.data.RA <- imb.data
imb.data.RA[, celltypes] <- t(microbiome::transform(Y, "compositional"))

# Mean relative abundances 
df.mean.RA <- data.frame(Celltype=names(colMeans(imb.data.RA[, celltypes])), mean.RA=colMeans(imb.data.RA[, celltypes]))

# Model analysis
# Adjustment: Sex + Age + Race + Ethnicity
design.adj <- model.matrix(~ Exposure + Sex + Age + Race + Ethnicity, data=imb.data)

v.adj <- voomCLR(counts=Y, design=design.adj, varCalc="analytical", varDistribution="NB", plot=FALSE, span=0.8)
fit.adj <- lmFit(v.adj, design.adj)
fit.adj <- eBayes(fit.adj)

# Adjusted model results
fit.adj$coefficients
fit.adj$p.value
fit.adj$df.residual
fit.adj$df.prior

# Exposure outcome
tt.adj.exposure <- topTableBC(fit.adj, coef="ExposureACPA+", n=Inf)
tt.adj.exposure 

# Arrange results
results.adj.mode <- tt.adj.exposure %>%
  rownames_to_column(var="Celltype") %>%
  left_join(df.mean.RA, by="Celltype") %>%
  mutate(Celltype = gsub("_", " ", Celltype),
         Celltype = ifelse(Celltype %in% c("CD56bright NK", "CD56dim NK"), 
                           paste0(Celltype, " cell"), Celltype)) %>%
  arrange(adj.P.Val) %>% 
  mutate(Celltype = factor(Celltype, levels = rev(unique(Celltype))),
         Category = case_when(grepl("B cell|Plasma", Celltype) ~ "B cell",
                              grepl("T cell|CD8|CD4|Treg|gdT|MAIT", Celltype) ~ "T cell",
                              grepl("NK", Celltype) ~ "NK cell",
                              grepl("monocyte", Celltype) ~ "Monocyte",
                              grepl("DC|ASDC", Celltype) ~ "Dendritic cell",
                              TRUE ~ "Other"),
         Significance = ifelse(adj.P.Val < 0.05, "Significant", "Not significant"),
         Asterix = ifelse(Significance=="Significant", "*", ""))

# Visualize results
ggplot(results.adj.mode, aes(x=logFC, y=Celltype)) +
  geom_vline(xintercept=0, linetype="dashed") +
  geom_segment(aes(x=0, xend=logFC, y=Celltype, yend=Celltype, color=Category)) +
  geom_point(aes(color=Category, size=mean.RA)) +
  geom_text(aes(label=Asterix, hjust=ifelse(logFC >= 0, -1.0, 1.5)), vjust=0.75, size=7, show.legend=FALSE) +
  guides(size="none", color=guide_legend(nrow=1)) +
  scale_x_continuous(limits=c(-0.75, 0.75), breaks=seq(-0.75, 0.75, by=0.25)) +
  labs(x="LogFC", y="Cell type") +
  theme_bw(base_size=16) +
  theme(legend.position="bottom")

## Voting for bias correction

# Previous adjusted model results
beta.adj <- fit.adj$coefficients[, "ExposureACPA+"]
se.adj <- fit.adj$stdev.unscaled[, "ExposureACPA+"] * sqrt(fit.adj$s2.post)

# Voting procedure
votes.adj <- matrix(NA, nrow=length(celltypes), ncol=length(celltypes))
colnames(votes.adj) <- names(beta.adj)
rownames(votes.adj) <- names(beta.adj)

for (p in 1:length(celltypes)) {
  votes.adj[p, ] <- beta.adj - beta.adj[p]
}

vote.mode.beta.adj <- apply(votes.adj, 2, function(x) mlv(x, method="meanshift", kernel="gaussian"))
names(vote.mode.beta.adj) <- names(beta.adj)

# Ensure similar output as previously
vote.mode.t.adj <- vote.mode.beta.adj / se.adj
vote.mode.pval.adj <- 2 * pt(-abs(vote.mode.t.adj), df=fit.adj$df.total)
vote.mode.adj.pval.adj <- p.adjust(vote.mode.pval.adj, method="BH")

tt.adj.vote.exposure <- tt.adj.exposure
tt.adj.vote.exposure$logFC <- vote.mode.beta.adj[rownames(tt.adj.exposure)]
tt.adj.vote.exposure$t <- vote.mode.t.adj[rownames(tt.adj.exposure)]
tt.adj.vote.exposure$P.Value <- vote.mode.pval.adj[rownames(tt.adj.exposure)]
tt.adj.vote.exposure$adj.P.Val <- vote.mode.adj.pval.adj[rownames(tt.adj.exposure)]

# Arrange results
results.adj.vote <- tt.adj.vote.exposure %>%
  rownames_to_column(var="Celltype") %>%
  left_join(df.mean.RA, by="Celltype") %>%
  mutate(Celltype = gsub("_", " ", Celltype),
         Celltype = ifelse(Celltype %in% c("CD56bright NK", "CD56dim NK"), 
                           paste0(Celltype, " cell"), Celltype)) %>%
  arrange(adj.P.Val) %>% 
  mutate(Celltype = factor(Celltype, levels = rev(unique(Celltype))),
         Category = case_when(grepl("B cell|Plasma", Celltype) ~ "B cell",
                              grepl("T cell|CD8|CD4|Treg|gdT|MAIT", Celltype) ~ "T cell",
                              grepl("NK", Celltype) ~ "NK cell",
                              grepl("monocyte", Celltype) ~ "Monocyte",
                              grepl("DC|ASDC", Celltype) ~ "Dendritic cell",
                              TRUE ~ "Other"),
         Significance = ifelse(adj.P.Val < 0.05, "Significant", "Not significant"),
         Asterix = ifelse(Significance=="Significant", "*", ""))

# Visualize results
ggplot(results.adj.vote, aes(x=logFC, y=Celltype)) +
  geom_vline(xintercept=0, linetype="dashed") +
  geom_segment(aes(x=0, xend=logFC, y=Celltype, yend=Celltype, color=Category)) +
  geom_point(aes(color=Category, size=mean.RA)) +
  geom_text(aes(label=Asterix, hjust=ifelse(logFC >= 0, -1.0, 1.5)), vjust=0.75, size=7, show.legend=FALSE) +
  guides(size="none", color=guide_legend(nrow=1)) +
  scale_x_continuous(limits=c(-0.75, 0.75), breaks=seq(-0.75, 0.75, by=0.25)) +
  labs(x="LogFC", y="Cell type") +
  theme_bw(base_size=16) +
  theme(legend.position="bottom")