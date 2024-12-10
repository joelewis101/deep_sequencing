# AMR genes

`../amr-gene-logreg.stan` is the standard random effect logistic regression
model (AMR gene presence/absence)
`../amr-gene-multilevel-gp-cholfactv2.stan` is the final gaussian process
logistic regression model

- `model-amr.R` fits and saves the models  
- `model-amr-diagnostics?(-GP).qmd` pulls and fits model diagnostics for the two
models  
- `building-correlated-random-eff-mod.html` outlines the gaussian process model and
the results from it 
- `resistome.qmd` is a description of the dataset and outputs of the random effect
models, focussing on the Gaussian process model (i.e. the final model).

# Taxonomy

The taxonomy models are

- `negbin-randintercept.stan` negative binomial absolute read count model with random
- intercept by participant
- `negbin_gp_nopreds` negative binomial absolute read count model with gaussian process
- within-participant correlation with no predictor variables
- `negbin_gp` negative binomial absolute read count model with gaussian process
- within-participant correlation
- `poisson_gp.stan` poisson bsolute read count model with gaussian process
- within-participant correlation (v bad fit)

`fit_taxonomy_models.R` will fit the `negbin_gp` models and save the output
`model-taxonomy-diagnostics-GP.qmd` will use the output to plot the model
diagnostics
`taxonomy.qmd` will plot the results

# Other scripts

`simulate_from_fitted_models.R` will run some simulations from the posterior and
generate the plots used in the FIS/ECCMID abstract,



