# General description
Spatiotemporal variability in the association between mental illness and substance use mortality and unemployment in the contiguous US

# Basic information
- This repository includes analysis codes of the published article in _Applied Geography_ [[Link](https://doi.org/10.1016/j.apgeog.2022.102664)]
- The outcome data was from [the Institute of Health Metrics and Evaluation (IHME)](https://www.healthdata.org/)
- The covariate data were from United States Census Bureau and [U.S. Bureau of Labor Statistics](https://www.bls.gov/lau/)
- The main analysis was done with R and R packages INLA, dplyr, tidyr, spdep, and Matrix.

# Notice
- I strongly recommend users to activate [PARDISO solver](https://pardiso-project.org/r-inla/) in [R-INLA](https://www.r-inla.org) for spatiotemporal model fitting to reduce time to fit models.

# Repository structure
- `./Code`: code files for the main analysis
- `./Data`: the data for the analysis (cleaned)

