#### Introduction to R - Al-Zghoul Lab, JUST
#### Shadi Shahatit - RA, 2026
# Libraries ---------------------------------------------------------------

# Install packages (run only once if they are not installed)

# install.packages("tidyverse")
# install.packages("ggplot2")
# install.packages("dplyr")
# install.packages("readr")
# install.packages("ggsci")
# install.packages("patchwork")
# install.packages("readxl")
# install.packages("BiocManager")
# or at once:
# install.packages(c("package_name","package_name"))

# Load libraries

library(tidyverse)
library(ggplot2)
library(dplyr)
library(readr)
library(ggsci)
library(patchwork)
library(readxl)

# General Info ------------------------------------------------------------

# <- assigns values to objects (aka, =)
# %>% passes output from one function into the next (aka, pipe |)
# naming convention: use clear short meaningful names, NO spaces, NO start with numbers

# Variables and Structures ------------------------------------------------

# variable types in R

# Atomic types (basic values)
numeric   <- 3.14        # real numbers
integer   <- 5L          # whole numbers
character <- "control"   # text strings
logical   <- TRUE        # TRUE / FALSE

# Data structures (collections of values)
vector  <- c(1, 2, 3)                   # same type, 1D
matrix  <- matrix(1:6, nrow = 2)        # same type, 2D
list    <- list(1, "A", TRUE)           # mixed types, 1D
dataframe <- data.frame(                # mixed types, 2D
  group = c("A", "B"),
  value = c(10, 20))

factor <- factor(c("control", "treated")) # categorical variable with levels

# Common variable: vectors, dataframes, and lists

# dataframes examples

experiment_df <- data.frame(
  group = c("control", "treated_low", "treated_mid", "treated_high"),
  weight = c(62, 58, 54, 64))

str(experiment_df)
head(experiment_df)
mean(experiment_df$weight)
summary(experiment_df)
nrow(experiment_df)
colnames(experiment_df)

# $ selects a column from a dataframe
experiment_df$group
experiment_df[experiment_df$group == "control",]

# Data Structure Example --------------------------------------------

# biological context:
# weight change using different doses of chemotherapy for cancer

# create your own data

set.seed(123)
control <- rnorm(20, mean = 65, sd = 4)          # no treatment
chemo_low <- rnorm(20, mean = 60, sd = 4)        # low-dose chemotherapy
chemo_high <- rnorm(20, mean = 52, sd = 4)       # high-dose chemotherapy
group <- c(
  rep("ctrl", 20),
  rep("low_dose", 20),
  rep("high_does", 20))
weight <- c(control, chemo_low, chemo_high)
example_weight_df <- data.frame(group, weight)

## OR

# load your data from another source (e.g., Excel)

# example_weight_df <- read_excel("file_path.xlsx")

# descriptive statistics

stats <- example_weight_df %>%
  group_by(group) %>%
  summarise(
    mean_weight = mean(weight),
    sd_weight = sd(weight))

# hypothesis testing

# t test

t_test_result <- t.test(
  weight[group == "ctrl"],
  weight[group == "high_does"])

t_test_result

# anova

anova_result <- aov(weight ~ group, data = example_weight_df)
summary(anova_result)

# Visualization ------------------------------------------------------------

# point plot

ggplot(example_weight_df, aes(x = group, y = weight, color = group))+
  geom_point(position = position_jitter(width = 0.15))+
  labs(x = "Treatment", y = "Weight",
       # title = "Individual Measurements Across Treatment Groups"
       )+
  theme_minimal()

# bar plot

ggplot(example_weight_df, aes(x = group, y = weight, fill = group))+
  stat_summary(fun = mean, geom = "bar")+
  labs(x = "Treatment", y = "Weight")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# box plot

ggplot(example_weight_df, aes(x = group, y = weight, fill = group)) +
  geom_boxplot() +
  labs(
    title = "Distribution of weight in cancer treatment groups",
    x = "Treatment",
    y = "Weight (g)")+
  theme_classic()


