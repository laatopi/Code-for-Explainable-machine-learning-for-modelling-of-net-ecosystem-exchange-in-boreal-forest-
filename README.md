# Greenhouse Gas Flux and Climate Data Analysis

This repository contains R scripts for processing and modeling SMEAR station data to analyze greenhouse gas flux and climate interactions.
Read the paper for detailed results:
[Paper](https://egusphere.copernicus.org/preprints/2023/egusphere-2023-2559/egusphere-2023-2559.pdf)
## Repository Structure

##### ├── config.R                   # Configures project paths
##### ├── data/                      # Processed data used for training the models
##### ├── rawData/                   # Contains raw input data
##### ├── results/                   # Stores model results and evaluation metrics
##### ├── scripts/                   # Contains processing and modeling scripts
##### │   ├── data_manipulation.R    # Data loading and processing
##### │   ├── single_site_models.R   # Model Training for Single Site Models
##### │   └── combined_site_models.R # Model training for combined Site Models
##### └── README.md                  # Project overview


## Setup and Requirements

- Install R and required packages: `dplyr`, `purrr`, `stringr`, `caret`, `Cubist`, `Metrics`, `here`.
- Set paths in `config.R` as needed.
- Run the model training scripts.
- Plot the results! (Codes not included in repository)
