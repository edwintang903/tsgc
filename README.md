# Time Series Growth Curves (tsgc) Package

<!-- badges: start -->
<!-- badges: end -->


## Overview

The `tsgc` package is designed for forecasting epidemics, including the detection of new waves and turning points, using a dynamic Gompertz model. It is suitable for predicting future values of variables that, when cumulated, are subject to some unknown saturation level. A leading-indicator extension allows a related series that moves ahead of the variable of interest, such as cases anticipating hospitalisations, to improve short-horizon forecasts, with model and lag selection supported by rolling-origin cross-validation. The underlying state space framework applies equally to data observed at daily, weekly, monthly, quarterly, or annual frequency. This approach is not only applicable to epidemics but also to domains like the diffusion of new products, thanks to its flexibility in adapting to changes in social behavior and policy. The `tsgc` package is demonstrated using COVID-19 confirmed cases and hospitalisation data, as well as non-epidemic applications.

## Installation

To install the latest version of the `tsgc` package from GitHub, use the following R command:
  
```r
# Install from GitHub
  install.packages("devtools")
  library(devtools)
  devtools::install_github("edwintang903/tsgc")
```

or install from the locally downloaded package as:
  
```r
devtools::install()
```

## Usage

Here is a basic example of setting up and estimating a model with the `tsgc` package:
  
  ```r
library(tsgc)

# Load example data
data("gauteng", package = "tsgc")

# Initialize and estimate the model
model <- SSModelDynamicGompertz$new(Y = gauteng)
results <- model$estimate()

# View results
print(results)
```

A leading-indicator model can be set up in much the same way, using a
second series that moves ahead of the variable of interest:

```r
library(tsgc)

# Load example data: daily cases and hospitalisations for England
data("england", package = "tsgc")

# Initialize and estimate the model, with cases as a 4-day leading
# indicator for hospitalisations
model <- SSModelLeadingIndicator$new(Y = england[, 1:2], n.lag = 4)
results <- model$estimate()

# View results
print(results)
```

## Features

- **Dynamic Gompertz Model:** Implements time series growth curve methods based on a dynamic Gompertz model.
- **Leading Indicator Model:** Models a target series jointly with a related series that leads it by a specified lag, such as cases anticipating hospitalisations, to improve short-horizon forecasts.
- **Model Comparison and Lag Selection:** Rolling-origin cross-validation (`cross_val()`) supports comparing candidate models and selecting the leading-indicator lag empirically.
- **Forecast Accuracy Tools:** Reports standard scale-dependent and scale-free accuracy measures (MAPE, sMAPE, MAE, RMSE) against a holdout sample.
- **Reproduction Number Estimation:** Translates growth-rate estimates into instantaneous reproduction number ($R_t$) estimates and credible intervals (`estimate_r0()`).
- **Flexible Applications:** While focused on epidemics, `tsgc` is also applicable in other areas, such as marketing, product diffusion, and other domains with a growth-curve-like trajectory.
- **Multiple Data Frequencies:** The same state space framework applies to data observed daily, weekly, monthly, quarterly, or annually.
- **Reinitialization Support:** Offers functionalities for reinitializing data series and model to account for multiple waves in an epidemic.
- **Comprehensive Documentation:** Includes detailed examples and vignettes to guide users through forecasting exercises.

## Dependencies

This package requires R (version 3.5.0 or higher) and depends on several other R packages for handling state space models and time series data, including `KFAS`, `xts`, `zoo`, and `here`.

## Getting Help

For detailed documentation and examples, refer to the package's vignettes. Should you encounter any issues or have questions, please file them in the GitHub Issues section of the `tsgc` repository.

## Contributing

Contributions to `tsgc` are welcome, including bug reports, feature requests, and pull requests. Please see the GitHub repository for contribution guidelines.

## License

This package is released under the GNU General Public License v3.0.

## Citation

If you use the `tsgc` package in your research, please cite it as follows:

Ashby, M., Harvey, A., Kattuman, P., & Thamotheram, C. (2021). Forecasting epidemic trajectories: Time Series Growth Curves package `tsgc`. Cambridge Centre for Health Leadership & Enterprise. URL: [https://www.jbs.cam.ac.uk/wp-content/uploads/2024/03/cchle-tsgc-paper-2024.pdf]

## Acknowledgments

Our gratitude goes to the Cambridge Centre for Health Leadership & Enterprise, University of Cambridge Judge Business School, and Public Health England/UK Health Security Agency for their support. Special thanks to Thilo Klein and Stefan Scholtes for their constructive comments, and to all contributors to the development and documentation of the `tsgc` package.