library(dplyr)
library(tidyr)
library(ggplot2)
library(robustbase)
library(microbiome)

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

# CLR transformation
CLR <- function(Y=Y) { 
  geo.mean <- exp(rowMeans(log(Y))) 
  Z <- t(log(t(Y) / geo.mean))
}

################################################################################

# EXPLORATORY DATA ANALYSIS: ORIGINAL DATA

## Covariate data exploration

# Check for duplicates
length(unique(ori.data$Subject)) == nrow(ori.data)

# Categorical variable counts
table(ori.data$Group)
table(ori.data$Exposure)

table(ori.data$Sex)
table(ori.data$Exposure, ori.data$Sex)

table(ori.data$Race)
table(ori.data$Exposure, ori.data$Race)

table(ori.data$Ethnicity)
table(ori.data$Exposure, ori.data$Ethnicity)

# Categorical variable proportions
prop.table(table(ori.data$Group))
prop.table(table(ori.data$Exposure))

prop.table(table(ori.data$Sex))
prop.table(table(ori.data$Exposure, ori.data$Sex), margin=1)

prop.table(table(ori.data$Race))
prop.table(table(ori.data$Exposure, ori.data$Race), margin=1)

prop.table(table(ori.data$Ethnicity))
prop.table(table(ori.data$Exposure, ori.data$Ethnicity), margin=1)

# Age
summary(ori.data$Age)
by(ori.data$Age, ori.data$Exposure, summary)

# BMI
summary(ori.data$BMI)
by(ori.data$BMI, ori.data$Exposure, summary)

## Cell count data exploration

# Observed counts
Y <- as.matrix(ori.data[, celltypes]) 

# CLR-transformed counts
Z <- CLR(Y + 0.5) 

ori.data.CLR <- ori.data
ori.data.CLR[, celltypes] <- Z

# Relative abundances
P <- t(microbiome::transform(t(Y), "compositional"))

ori.data.RA <- ori.data
ori.data.RA[, celltypes] <- P

# Sorted cell types
sorted.celltypes <- names(sort(colMedians(Y)))
sorted.celltypes.CLR <- names(sort(colMedians(Z)))
sorted.celltypes.RA <- names(sort(colMedians(P)))

# Delete underscores in sorted cell type labels
sorted.celltypes <- gsub("_", " ", names(sort(colMedians(Y))))
sorted.celltypes.CLR <- gsub("_", " ", names(sort(colMedians(Z))))
sorted.celltypes.RA <- gsub("_", " ", names(sort(colMedians(P))))

# Check if sorted cell types are equal
all(sorted.celltypes==sorted.celltypes.CLR, sorted.celltypes==sorted.celltypes.RA)

# Pivot to long format
ori.data.long <- ori.data %>%
  pivot_longer(cols=all_of(celltypes), names_to="Celltype", values_to="Y") %>%
  mutate(Celltype = gsub("_", " ", Celltype),
         Celltype = factor(Celltype, levels=sorted.celltypes))

ori.data.CLR.long <- ori.data.CLR %>%
  pivot_longer(cols=all_of(celltypes), names_to="Celltype", values_to="Z") %>%
  mutate(Celltype = gsub("_", " ", Celltype),
         Celltype = factor(Celltype, levels=sorted.celltypes.CLR))

# Observed count barchart
ggplot(ori.data.long, aes(x=Exposure, y=Y, fill=Celltype)) + 
  geom_bar(position="fill", stat="identity") +
  labs(y="Mean relative abundance", fill="Cell type") +
  theme_bw(base_size=16)

# CLR-transformed count boxplots
ggplot(ori.data.CLR.long, aes(x=Celltype, y=Z, fill=Celltype)) +
  geom_boxplot(show.legend=FALSE) +
  facet_wrap(~Exposure) +
  labs(x="Cell type", y="CLR-transformed count") +
  theme_bw(base_size=16) +
  theme(axis.text.x = element_text(angle=45, hjust=1, size=10))

ggplot(ori.data.CLR.long, aes(x=Celltype, y=Z, fill=Exposure)) +
  geom_boxplot() +
  labs(x="Cell type", y="CLR-transformed count") +
  theme_bw(base_size=16) +
  theme(axis.text.x = element_text(angle=45, hjust=1))

ggplot(ori.data.CLR.long, aes(x=Celltype, y=Z, fill=Exposure)) +
  geom_boxplot() +
  facet_wrap(~Sex) +
  labs(x="Cell type", y="CLR-transformed count") +
  theme_bw(base_size=16) +
  theme(axis.text.x = element_text(angle=45, hjust=1, size=10))

ggplot(ori.data.CLR.long, aes(x=Celltype, y=Z, fill=Exposure)) +
  geom_boxplot() +
  facet_wrap(~Ethnicity) +
  labs(x="Cell type", y="CLR-transformed count") +
  theme_bw(base_size=16) +
  theme(axis.text.x = element_text(angle=45, hjust=1, size=10))

ggplot(ori.data.CLR.long, aes(x=Celltype, y=Z, fill=Exposure)) +
  geom_boxplot() +
  facet_wrap(~Race) +
  labs(x="Cell type", y="CLR-transformed count") +
  theme_bw(base_size=16) +
  theme(axis.text.x = element_text(angle=45, hjust=1, size=10))

# Shannon index
div <- microbiome::diversity(t(Y))
div.df <- data.frame(Shannon = div$shannon,
                     Exposure = ori.data$Exposure,
                     Age = ori.data$Age,
                     BMI = ori.data$BMI)

ggplot(div.df, aes(x=Age, y=Shannon, color=Exposure, fill=Exposure)) + 
  geom_point() + 
  geom_smooth(method="loess", se=TRUE) + 
  labs(x="Age (year)", y="Shannon", color="Exposure") +
  theme_bw(base_size=16)

ggplot(div.df, aes(x=BMI, y=Shannon, color=Exposure, fill=Exposure)) + 
  geom_point() + 
  geom_smooth(method="loess", se=TRUE) + 
  labs(x="BMI (kg/m²)", y="Shannon", color="Exposure") +
  theme_bw(base_size=16)

################################################################################

# EXPLORATORY DATA ANALYSIS: IMBALANCED DATA

## Covariate data exploration

# Check for duplicates
length(unique(imb.data$Subject)) == nrow(imb.data)

# Categorical variable counts
table(imb.data$Group)
table(imb.data$Exposure)

table(imb.data$Sex)
table(imb.data$Exposure, imb.data$Sex)

table(imb.data$Race)
table(imb.data$Exposure, imb.data$Race)

table(imb.data$Ethnicity)
table(imb.data$Exposure, imb.data$Ethnicity)

# Categorical variable proportions
prop.table(table(imb.data$Group))
prop.table(table(imb.data$Exposure))

prop.table(table(imb.data$Sex))
prop.table(table(imb.data$Exposure, imb.data$Sex), margin=1)

prop.table(table(imb.data$Race))
prop.table(table(imb.data$Exposure, imb.data$Race), margin=1)

prop.table(table(imb.data$Ethnicity))
prop.table(table(imb.data$Exposure, imb.data$Ethnicity), margin=1)

# Age
summary(imb.data$Age)
by(imb.data$Age, imb.data$Exposure, summary)

# BMI
summary(imb.data$BMI)
by(imb.data$BMI, imb.data$Exposure, summary)

## Cell count data exploration

# Observed counts
Y <- as.matrix(imb.data[, celltypes]) 

# CLR-transformed counts
Z <- CLR(Y + 0.5) 

imb.data.CLR <- imb.data
imb.data.CLR[, celltypes] <- Z

# Relative abundances
P <- t(microbiome::transform(t(Y), "compositional"))

imb.data.RA <- imb.data
imb.data.RA[, celltypes] <- P

# Sorted cell types
sorted.celltypes <- names(sort(colMedians(Y)))
sorted.celltypes.CLR <- names(sort(colMedians(Z)))
sorted.celltypes.RA <- names(sort(colMedians(P)))

# Delete underscores in sorted cell type labels
sorted.celltypes <- gsub("_", " ", names(sort(colMedians(Y))))
sorted.celltypes.CLR <- gsub("_", " ", names(sort(colMedians(Z))))
sorted.celltypes.RA <- gsub("_", " ", names(sort(colMedians(P))))

# Check if sorted cell types are equal
all(sorted.celltypes==sorted.celltypes.CLR, sorted.celltypes==sorted.celltypes.RA)

# Pivot to long format
imb.data.long <- imb.data %>%
  pivot_longer(cols=all_of(celltypes), names_to="Celltype", values_to="Y") %>%
  mutate(Celltype = gsub("_", " ", Celltype),
         Celltype = factor(Celltype, levels=sorted.celltypes))

imb.data.CLR.long <- imb.data.CLR %>%
  pivot_longer(cols=all_of(celltypes), names_to="Celltype", values_to="Z") %>%
  mutate(Celltype = gsub("_", " ", Celltype),
         Celltype = factor(Celltype, levels=sorted.celltypes.CLR))

# Observed count barchart
ggplot(imb.data.long, aes(x=Exposure, y=Y, fill=Celltype)) + 
  geom_bar(position="fill", stat="identity") +
  labs(y="Mean relative abundance", fill="Cell type") +
  theme_bw(base_size=16)

# CLR-transformed count boxplots
ggplot(imb.data.CLR.long, aes(x=Celltype, y=Z, fill=Celltype)) +
  geom_boxplot(show.legend=FALSE) +
  facet_wrap(~Exposure) +
  labs(x="Cell type", y="CLR-transformed count") +
  theme_bw(base_size=16) +
  theme(axis.text.x = element_text(angle=45, hjust=1, size=10))

ggplot(imb.data.CLR.long, aes(x=Celltype, y=Z, fill=Exposure)) +
  geom_boxplot() +
  labs(x="Cell type", y="CLR-transformed count") +
  theme_bw(base_size=16) +
  theme(axis.text.x = element_text(angle=45, hjust=1))

ggplot(imb.data.CLR.long, aes(x=Celltype, y=Z, fill=Exposure)) +
  geom_boxplot() +
  facet_wrap(~Sex) +
  labs(x="Cell type", y="CLR-transformed count") +
  theme_bw(base_size=16) +
  theme(axis.text.x = element_text(angle=45, hjust=1, size=10))

ggplot(imb.data.CLR.long, aes(x=Celltype, y=Z, fill=Exposure)) +
  geom_boxplot() +
  facet_wrap(~Ethnicity) +
  labs(x="Cell type", y="CLR-transformed count") +
  theme_bw(base_size=16) +
  theme(axis.text.x = element_text(angle=45, hjust=1, size=10))

ggplot(imb.data.CLR.long, aes(x=Celltype, y=Z, fill=Exposure)) +
  geom_boxplot() +
  facet_wrap(~Race) +
  labs(x="Cell type", y="CLR-transformed count") +
  theme_bw(base_size=16) +
  theme(axis.text.x = element_text(angle=45, hjust=1, size=10))

# Shannon index
div <- microbiome::diversity(t(Y))
div.df <- data.frame(Shannon = div$shannon,
                     Exposure = imb.data$Exposure,
                     Age = imb.data$Age,
                     BMI = imb.data$BMI)


ggplot(div.df, aes(x=Age, y=Shannon, color=Exposure, fill=Exposure)) + 
  geom_point() + 
  geom_smooth(method="loess", se=TRUE) + 
  labs(x="Age (year)", y="Shannon", color="Exposure") +
  theme_bw(base_size=16)


ggplot(div.df, aes(x=BMI, y=Shannon, color=Exposure, fill=Exposure)) + 
  geom_point() + 
  geom_smooth(method="loess", se=TRUE) + 
  labs(x="BMI (kg/m²)", y="Shannon", color="Exposure") +
  theme_bw(base_size=16)