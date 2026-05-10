library(dplyr)
library(tidyr)
library(ggplot2)
library(cobalt)
library(compositions)
library(factoextra)

setwd("C:/Users/sarah/OneDrive - KU Leuven/Desktop/Master's Thesis")

# Load data
data <- readRDS("original_RA_data.rds")
data <- na.omit(data)
head(data)
summary(data)

################################################################################

# DOWNSAMPLING

# Separate exposure groups
data.neg <- data %>% filter(Exposure == "ACPA-")
data.pos <- data %>% filter(Exposure == "ACPA+")

# Define sampling weights
data.pos.weighted <- data.pos %>% 
  mutate(
    # Sex:
    # Large weight for male, small weight for female
    w.sex = ifelse(Sex == "Male", 1.0, 0.1),
    
    # Age:
    # Large weight for older age, small weight for younger age
    w.age = (Age - min(Age, na.rm=TRUE)) / (max(Age, na.rm=TRUE) - min(Age, na.rm=TRUE)),
    
    # BMI:
    # Large weight for lower BMI, small weight for higher BMI
    w.BMI = (max(BMI, na.rm=TRUE) - BMI) / (max(BMI, na.rm=TRUE) - min(BMI, na.rm=TRUE)),
    
    # Race:
    # Large weight for African American, Asian, other; small weight for Caucasian
    w.race = case_when(Race == "African American" ~ 1.0,
                       Race == "Asian" ~ 1.0,       
                       Race == "Other" ~ 1.0,       
                       Race == "Caucasian" ~ 0.1),       
    
    # Ethnicity
    # Large weight for non-Hispanic, small weight for Hispanic
    w.ethnicity = ifelse(Ethnicity == "Non-Hispanic", 1.0, 0.1),
    
    # Calculate final weight 
    final.weight = (w.sex * w.age * w.BMI * w.race * w.ethnicity) + 0.001
  )

# Downsample to 31 individuals
set.seed(123)
data.pos.conf <- data.pos.weighted %>%
  slice_sample(n = 31, weight_by=final.weight) %>%
  dplyr::select(-starts_with("w."), -final.weight) # Clean up the weight columns

# Recombine exposure groups
data.conf <- bind_rows(data.neg, data.pos.conf)

# Save downsampled RA data
saveRDS(data.conf, file="imbalanced_RA_data.rds")

################################################################################

# DATA VISUALIZATION

## Propensity scores

# Logistic regression
data <- data %>% mutate(Exposure.binary = ifelse(Exposure=="ACPA+", 1, 0))
ps.model.noconf <- glm(Exposure.binary ~ Sex + Age + BMI + Race + Ethnicity, 
                       family=binomial(link="logit"), data=data)

data.conf <- data.conf %>% mutate(Exposure.binary = ifelse(Exposure=="ACPA+", 1, 0))
ps.model.conf <- glm(Exposure.binary ~ Sex + Age + BMI + Race + Ethnicity, 
                     family=binomial(link="logit"), data=data.conf)

# Propensity scores
data$ps <- predict(ps.model.noconf, type="response")
data.conf$ps <- predict(ps.model.conf, type="response")

ggplot(data, aes(x=ps, fill=Exposure)) + 
  geom_density(alpha=0.6) + 
  labs(x="Propensity score", y="Density") +
  theme_bw(base_size=16)
ggplot(data.conf, aes(x=ps, fill=Exposure)) + 
  geom_density(alpha=0.6) + 
  labs(x="Propensity score", y="Density") +
  theme_bw(base_size=16)

## Standardized Mean Difference (SMD)

# Balance statistics
balance.noconf <- bal.tab(Exposure ~ Sex + Age + BMI + Race + Ethnicity, data=data)
balance.conf <- bal.tab(Exposure ~ Sex + Age + BMI + Race + Ethnicity, data=data.conf)

# SMD
SMD.noconf <- balance.noconf$Balance %>%
  mutate(Covariate=rownames(.)) %>%
  dplyr::select(Covariate, Diff.Un) %>%
  rename("Original"=Diff.Un)
SMD.conf <- balance.conf$Balance %>%
  mutate(Covariate=rownames(.)) %>%
  dplyr::select(Covariate, Diff.Un) %>%
  rename("Imbalanced"=Diff.Un)

SMD.joined <- SMD.noconf %>%
  left_join(SMD.conf, by="Covariate") %>%
  pivot_longer(cols=c("Original", "Imbalanced"), names_to="Data", values_to="SMD") %>%
  mutate(Covariate=factor(Covariate, levels=rev(unique(Covariate))))

# Directional SMD
ggplot(SMD.joined, aes(x=SMD, y=Covariate, color=Data)) +
  geom_vline(xintercept=0, linetype="solid") +
  geom_vline(xintercept=c(-0.1, 0.1), linetype="dashed") +
  geom_point(size=3.5) +
  labs(x="Standardized Mean Difference (SMD)", y="") +
  scale_y_discrete(labels=c("Sex_Male" = "Sex: male",
                            "Age" = "Age",
                            "BMI" = "BMI",
                            "Race_Caucasian" = "Race: Caucasian",
                            "Race_African American" = "Race: African American",
                            "Race_Asian" = "Race: Asian",
                            "Race_Other" = "Race: other",
                            "Ethnicity_Hispanic" = "Ethnicity: Hispanic")) +
  scale_color_manual(values = c("Original"="#F8766D", "Imbalanced"="#00BFC4"),
                     breaks = c("Original", "Imbalanced")) +
  theme_bw(base_size=16)

# Absolute values
ggplot(SMD.joined, aes(x=abs(SMD), y=Covariate, color=Data)) +
  geom_vline(xintercept=0, linetype="solid") +
  geom_vline(xintercept=0.1, linetype="dashed") +
  geom_point(size=3.5) +
  labs(x="Absolute Standardized Mean Difference (SMD)", y="") +
  scale_y_discrete(labels=c("Sex_Male" = "Sex: male",
                            "Age" = "Age",
                            "BMI" = "BMI",
                            "Race_Caucasian" = "Race: Caucasian",
                            "Race_African American" = "Race: African American",
                            "Race_Asian" = "Race: Asian",
                            "Race_Other" = "Race: other",
                            "Ethnicity_Hispanic" = "Ethnicity: Hispanic")) +
  scale_color_manual(values = c("Original"="#F8766D", "Imbalanced"="#00BFC4"),
                     breaks = c("Original", "Imbalanced")) +
  theme_bw(base_size=16)

## Principal Component Analysis (PCA)

# Cell types
celltypes <- c("ASDC", "CD14_monocyte", "CD16_monocyte", "CD56bright_NK", 
               "CD56dim_NK", "CD8aa", "DN_T_cell", "Effector_B_cell", 
               "Erythrocyte", "ILC", "Intermediate_monocyte", "MAIT", 
               "Memory_B_cell", "Memory_CD4_T_cell", "Memory_CD8_T_cell", 
               "Naive_B_cell", "Naive_CD4_T_cell", "Naive_CD8_T_cell", 
               "Plasma_cell", "Platelet", "Progenitor_cell", 
               "Proliferating_NK_cell", "Proliferating_T_cell", 
               "Transitional_B_cell", "Treg", "cDC1", "cDC2", "gdT", "pDC")

# Compositional PCA on cell type counts
comps <- acomp(data.conf[, celltypes]) 
pc.comps <- princomp(comps) 

# Add covariate data
pc.scores <- data.frame(pc.comps$scores)
pc.scores$Exposure <- data.conf$Exposure 
pc.scores$Sex <- data.conf$Sex
pc.scores$Age <- data.conf$Age
pc.scores$BMI <- data.conf$BMI
pc.scores$Race <- data.conf$Race 
pc.scores$Ethnicity <- data.conf$Ethnicity

# Scree plot
fviz_eig(pc.comps, addlabels=TRUE) 

# Biplots
ggplot(pc.scores, aes(x=Comp.1, y=Comp.2, col=Exposure)) + 
  geom_point(size=3) + 
  labs(x="PC1 (21.3%)", y="PC2 (12.8%)") +
  scale_color_discrete(breaks=c("ACPA-", "ACPA+"), 
                       labels=c("ACPA-        ", "ACPA+        ")) +
  theme_bw(base_size=16) +
  theme(aspect.ratio=1)

ggplot(pc.scores, aes(x=Comp.1, y=Comp.2, col=Sex, shape=Exposure)) + 
  geom_point(size=3) +
  labs(x="PC1 (21.3%)", y="PC2 (12.8%)") +
  scale_shape_discrete(breaks=c("ACPA-", "ACPA+"), 
                       labels=c("ACPA-        ", "ACPA+        ")) +
  theme_bw(base_size=16) +
  theme(aspect.ratio=1)

ggplot(pc.scores, aes(x=Comp.1, y=Comp.2, col=Age, shape=Exposure)) + 
  geom_point(size=3, alpha=0.8) + 
  scale_color_gradient(low="blue", high="red") +
  labs(x="PC1 (21.3%)", y="PC2 (12.8%)") +
  scale_shape_discrete(breaks=c("ACPA-", "ACPA+"), 
                       labels=c("ACPA-        ", "ACPA+        ")) +
  theme_bw(base_size=16) +
  theme(aspect.ratio=1) +
  guides(col=guide_colorbar(order=2), shape = guide_legend(order=1))

ggplot(pc.scores, aes(x=Comp.1, y=Comp.2, col=BMI, shape=Exposure)) + 
  geom_point(size=3, alpha=0.8) + 
  scale_color_gradient(low="blue", high="red") +
  labs(x="PC1 (21.3%)", y="PC2 (12.8%)") +
  scale_shape_discrete(breaks=c("ACPA-", "ACPA+"), 
                       labels=c("ACPA-        ", "ACPA+        ")) +
  theme_bw(base_size=16) +
  theme(aspect.ratio=1) +
  guides(col=guide_colorbar(order=2), shape = guide_legend(order=1))

ggplot(pc.scores, aes(x=Comp.1, y=Comp.2, col=Race, shape=Exposure)) + 
  geom_point(size=3) + 
  labs(x="PC1 (21.3%)", y="PC2 (12.8%)") +
  scale_color_discrete(breaks=c("Caucasian", "African American", "Asian", "Other"), 
                       labels=c("Caucasian", "A. American", "Asian", "Other")) +
  theme_bw(base_size=16) +
  theme(aspect.ratio=1)

ggplot(pc.scores, aes(x=Comp.1, y=Comp.2, col=Ethnicity, shape=Exposure)) + 
  geom_point(size=3) + 
  labs(x="PC1 (21.3%)", y="PC2 (12.8%)") +
  theme_bw(base_size=16) +
  theme(aspect.ratio=1)