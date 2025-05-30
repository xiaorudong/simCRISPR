# simCRISPR

**simCRISPR** is an R package for simulating pooled CRISPR screening data under various experimental designs. 
It models biological and technical variability, including sgRNA knockout effects, treatment and interaction effects, PCR amplification, and sequencing noise. 
The framework is designed to support method development and benchmarking in CRISPR analysis.


## Installation

You can install `simCRISPR` from GitHub using:

```{r}
if (!requireNamespace("devtools", quietly=TRUE))
    install.packages("devtools")
devtools::install_github("xiaorudong/simCRISPR")
```
## Example
Here is a basic example of how to use the `simCRISPR` package to generate pooled CRISPR screening data:

```{r}
library(simCRISPR)

raw_obj <- sim_crispr()
initial_data <- raw_obj$sim_raw
sim_data <- seq_add(initial_data)
```
