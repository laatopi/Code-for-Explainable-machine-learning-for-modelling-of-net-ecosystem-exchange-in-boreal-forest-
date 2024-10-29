# config.R

# Configuration file for paths used in the project
# Users can adjust these paths based on their local environment.

# Base path for data files (relative to project root)
base_data_path <- "./data/"

# Base path for output files (relative to project root)
base_output_path <- "./results/"

# Ensure data and results folders exist (optional but helpful)
if (!dir.exists(base_data_path)) dir.create(base_data_path)
if (!dir.exists(base_output_path)) dir.create(base_output_path)