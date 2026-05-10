library(dplyr)
library(tidyr)

setwd("C:/Users/sarah/OneDrive - KU Leuven/Desktop/Master's Thesis")

################################################################################

# PREPROCESSING

# Load raw RA data
frequencies.long <- read.csv("altra_AIFI_L2_frequencies.csv")
metadata <- read.csv("aifi_altra_sound-life_clinical_lab_results.csv")
h5ad <- read.csv("extracted_h5ad_metadata.csv")
race.ethnicity <- read.csv("2026-02-28_ALTRA_BRI_Race_Ethnicity.csv")

# Filter for first visit per subject
unique.subjects <- frequencies.long %>%
  arrange(sample.daysSinceFirstVisit) %>%
  dplyr::select(subject.subjectGuid, sample.sampleKitGuid, sample.daysSinceFirstVisit) %>%
  distinct() %>%
  group_by(subject.subjectGuid) %>%
  slice(1) %>%
  ungroup()

frequencies.long <- frequencies.long %>%
  semi_join(unique.subjects, by = c("subject.subjectGuid", "sample.sampleKitGuid"))

# Select important variables from frequency data
frequencies.long <- frequencies.long %>%
  dplyr::select(subject.subjectGuid, 
                sample.sampleKitGuid, 
                subject.diseaseGroup, 
                subject.biologicalSex, 
                sample.subjectAgeAtDraw,
                total_cells,
                AIFI_L2, 
                AIFI_L2_count)

# Pivot to wide format
frequencies.wide <- frequencies.long %>%
  pivot_wider(names_from=AIFI_L2, values_from=AIFI_L2_count, values_fill=list(AIFI_L2_count=0))

# Select important variables from metadata
metadata <- metadata %>%
  dplyr::select(sample.sampleKitGuid, 
                am.bmi) %>%
  distinct()

# Prepare race/ethnicity data 
race.ethnicity <- race.ethnicity %>%
  dplyr::select(sample.sampleKitGuid, subject.race, subject.ethnicity) %>%
  distinct()

colnames(frequencies.wide)

# Merge everything together
merged.data <- frequencies.wide %>%
  left_join(metadata, by="sample.sampleKitGuid") %>%
  left_join(h5ad, by="sample.sampleKitGuid") %>%
  left_join(race.ethnicity, by="sample.sampleKitGuid") %>%
  
  # Rename columns 
  rename(Subject = subject.subjectGuid,
         Kit = sample.sampleKitGuid,
         Group = subject.diseaseGroup,
         Sex = subject.biologicalSex,
         Age = sample.subjectAgeAtDraw,
         BMI = am.bmi,
         Batch = batch_id,
         Race = subject.race,
         Ethnicity = subject.ethnicity,
         ASDC = ASDC,
         CD14_monocyte = `CD14 monocyte`,
         CD16_monocyte = `CD16 monocyte`,
         CD56bright_NK = `CD56bright NK cell`,
         CD56dim_NK = `CD56dim NK cell`,
         CD8aa = CD8aa,
         DN_T_cell = `DN T cell`,
         Effector_B_cell = `Effector B cell`,
         Erythrocyte = Erythrocyte,
         ILC = ILC,
         Intermediate_monocyte = `Intermediate monocyte`,
         MAIT = MAIT,
         Memory_B_cell = `Memory B cell`,
         Memory_CD4_T_cell = `Memory CD4 T cell`,
         Memory_CD8_T_cell = `Memory CD8 T cell`,
         Naive_B_cell = `Naive B cell`,
         Naive_CD4_T_cell = `Naive CD4 T cell`,
         Naive_CD8_T_cell = `Naive CD8 T cell`,
         Plasma_cell = `Plasma cell`,
         Platelet = Platelet,
         Progenitor_cell = `Progenitor cell`,
         Proliferating_NK_cell = `Proliferating NK cell`,
         Proliferating_T_cell = `Proliferating T cell`,
         Transitional_B_cell = `Transitional B cell`,
         Treg = Treg,
         cDC1 = cDC1,
         cDC2 = cDC2,
         gdT = gdT,
         pDC = pDC,
         Total = total_cells) %>%
  
  # Add and arrange by exposure variable
  mutate(Exposure = if_else(Group=="Control (HC1)", "Healthy", "Diseased"),
         Exposure = factor(Exposure, levels=c("Healthy", "Diseased"), labels=c("ACPA-", "ACPA+"))) %>%
  arrange(Exposure) %>%
  
  # Convert categorical variables to factors
  mutate(Group = factor(Group, levels=c("Control (HC1)", "ACPA+ At-risk (ARI)", "ACPA+ Early RA (ERA)")),
         Sex = factor(Sex, levels=c("Female", "Male")),
         Race = factor(Race, levels=c("Caucasian", "African American", "Asian", "Other")),
         Ethnicity = factor(Ethnicity, levels=c("Non-Hispanic origin", "Hispanic or Latino origin"), labels=c("Non-Hispanic", "Hispanic"))) %>%
  
  # Relocate cell type, batch and exposure variables
  relocate(ASDC,
           Effector_B_cell, Memory_B_cell, Naive_B_cell, Transitional_B_cell, Plasma_cell,
           CD14_monocyte, CD16_monocyte, Intermediate_monocyte,
           cDC1, cDC2, pDC,
           CD56bright_NK, CD56dim_NK, Proliferating_NK_cell,
           CD8aa, DN_T_cell, MAIT, Memory_CD4_T_cell, Memory_CD8_T_cell,
           Naive_CD4_T_cell, Naive_CD8_T_cell, Proliferating_T_cell, Treg, gdT,
           Erythrocyte, ILC, Platelet, Progenitor_cell, Total,
           .after = last_col()) %>%
  relocate(Batch, .after=Kit) %>%
  relocate(Exposure, .after=Group)

# Get final data
data <- merged.data

# Check 
head(data)
summary(data)
dim(data)
levels(data$Exposure)

# Sample sizes
data %>%
  count(Exposure, name="n") %>%
  mutate(percentage = round((n/nrow(data)) * 100, 1))

# Complete case sample sizes
cc <- na.omit(data)
cc %>%
  count(Exposure, name="n") %>%
  mutate(percentage = round((n/nrow(cc)) * 100, 1))

# Save preprocessed RA data
saveRDS(data, file="original_RA_data.rds")