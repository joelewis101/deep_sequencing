# AMR genes

../amr-gene-logreg.stan is the standard random effect model
../amr-gene-multilevel-gp-cholfactv2.stan is the final gaussian process model

- model-amr.R fits and saves the models  
- model-amr-diagnostics?(-GP).qmd pulls and fits model diagnostics for the two
models  
- building-correlated-random-eff-mod.html outlines the gaussian process model and
the results from it 
- resistome.qmd is a description of the dataset and outputs of the random effect
models, focussing on the Gaussian process model (i.e. the final model).
