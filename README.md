# ER status and positive nodes in SEER vs METABRIC: a Bayesian comparison

Course project for STA510 Bayesian statistical methods (University of Gothenburg).

## Question

Are the associations of ER status and of the number of positive lymph nodes with death within
5 years the same in two breast cancer datasets, SEER (US registry) and METABRIC (UK/Canada cohort),
and by how much do they differ?

The difference is `delta = beta_METABRIC - beta_SEER` on the log-odds scale (one for ER, one for nodes);
`exp(delta)` is the ratio of the odds ratios in the two datasets.

## Models

Bayesian logistic regression (`brms`, Stan via `cmdstanr`) with weakly informative normal priors built from
published studies. Four models are compared:

| Model | ER and nodes | Age, size, grade | Pooling |
|---|---|---|---|
| M0 | one shared coefficient | shared | complete |
| MH | `beta_j = beta_0 + u_j`, `u_j ~ N(0, sigma^2)`, `sigma ~ exponential(2)` | shared | partial |
| M1 (main) | separate per dataset | shared | none |
| M2 | separate per dataset | separate per dataset | none |

Checks: prior predictive check, posterior predictive checks by group, convergence diagnostics and
leave-one-out cross-validation.

## Main result

No clear evidence that the ER or node association differs between SEER and METABRIC; a moderate
difference is not excluded. Under M1, the median `delta` for ER is 0.15 (95% interval -0.31 to 0.62),
i.e. `exp(delta)` = 1.16 (0.73 to 1.86); for nodes 0.07 (-0.21 to 0.36). All 95% intervals include 0.
The models predict about equally well (loo differences below 4).

## Repository contents

```
project.R        full analysis: data, priors, models, fitting, checks, comparison
data/            raw data go here (not included; see data/README.md)
figures/         posterior predictive check figures and trace plots
report/          report.pdf (LaTeX source can be added here)
slides/          presentation.pdf
```

## Running the analysis

1. Download the data into `data/` (see `data/README.md`).
2. Install the R packages `brms`, `cmdstanr` (with CmdStan), `ggplot2` and `readxl`.
3. Run `project.R` from the repository root (or set the environment variable `PROJECT` to that folder).
   Fitted models are cached by `brms` as `fit_*.rds` files in the working directory, and figures are
   written to `figures/`.

The script is organised in numbered sections: setup, data, priors and models, prior predictive check, fit,
convergence, estimands (delta), posterior predictive checks, model comparison.

## Possible extensions

- Log tumour size and a flexible age term (the linear terms fit only approximately).
- Partial pooling of the age slope across the two datasets (between M1 and M2), to see whether it improves prediction.
- Prior sensitivity analysis (other rates for the exponential prior on sigma, narrower or wider priors on the coefficients).
- Repeat the analysis with the unrestricted METABRIC sample.

## References

The priors are based on published sources, listed in the report (`report/report.pdf`).
