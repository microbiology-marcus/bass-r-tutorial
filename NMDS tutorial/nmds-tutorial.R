# NMDS tutorial

# Author:   Marcus E. Johnson
# Email:    marcus.johnson@ec.gc.ca
# Date:     2024-01-15

# Copyright (c) Environment and Climate Change Canada, 2024

## Purpose:
# This R script accompanies the tutorial "NMDS analysis of sequencing data"
# and uses the example dataset from phyloseq::GlobalPatterns by Caporaso et al 
# available at https://doi.org/10.1073/pnas.1000080107 

## Notes:
# The section headers match the accompanying documnet and may not be ordered
# correctly numerically in the document. Missing header numbers do not indicate
# a missing section.


# Start-up ----------------------------------------------------------------
## 1.1 Required packages ----

library(vegan)
library(tidyverse)

## 1.3 Included example ----

df_otu <- read.table(file = "example_data/counts.txt",
                     sep = "\t", header = T, row.names = 1)
df_meta <- read.table(file = "example_data/metadata.txt",
                      sep = "\t", header = T, row.names = 1)

head(df_meta)


# Analysis ----------------------------------------------------------------
## 2.1 Correctly formatting data ----

# check data type
class(df_meta$Description)
class(df_meta$Location)

# convert test variables to factors
df_meta$Description <- as.factor(df_meta$Description)
df_meta$Location <- as.factor(df_meta$Location)

# example of an ordered list
location_order <- c("Human", "Environment", "Mock") # list out all factors
ordered(df_meta$Location, location_order)

## 2.2 Conducting analysis ----
# Bray-Curtis analysis
bray_curtis_dist <- vegdist(t(df_otu), method = "bray")

# conduct ANOVA using dist. matrix
permanova_result <- adonis2(                 
  formula = bray_curtis_dist ~ Description,          # set formula of variables
  data = df_meta)                    # include the relevant metadata data frame

# example ANOVA with more than one variable
example_permanova <- adonis2(                 
  formula = bray_curtis_dist ~ Description + Location,   # formula of variables
  data = df_meta)                    # include the relevant metadata data frame

print(permanova_result) # collect results

# conduct an NMDS on the Bray-Curtis matrix
nmds_result <- metaMDS(bray_curtis_dist) 

# collect the scores from the nmds
nmds_data <- as.data.frame(scores(nmds_result))
nmds_data <- cbind(nmds_data, df_meta) # bind the columns from the meta data


# Graphing the results ----------------------------------------------------
## 3.1 Generating graphs ----

# create ggplot
nmds_plot <- ggplot(data = nmds_data, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(col = Description, shape = Location), size = 2) +
  theme_bw()

nmds_plot # view the plot

## 3.2 Exporting graphs ----

ggsave(
  filename = "output/NMDS_scatterplot.PNG", # within the folder "output"
  plot = nmds_plot,
  width = 6.5, height = 4, 
  dpi = 300)
