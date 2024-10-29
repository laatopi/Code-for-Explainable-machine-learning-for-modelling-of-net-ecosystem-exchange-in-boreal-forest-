# config.R

# Load the `here` package for consistent relative paths
if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
library(here)

# Define paths relative to the project root using `here()`
raw_data_path <- here("rawData")
base_data_path <- here("data")
base_output_path <- here("results")

# Ensure required folders exist
if (!dir.exists(raw_data_path)) dir.create(raw_data_path)
if (!dir.exists(base_output_path)) dir.create(base_output_path)
