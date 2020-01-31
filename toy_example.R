observed_locs = matrix(rnorm(1000), ncol = 2)*10
predicted_locs = matrix(rnorm(1000), ncol = 2)*10
observed_field = t(chol(GpGp::matern_isotropic(locs =  observed_locs, c(1, 1, 1, 0))))%*%rnorm(500) + rnorm(500)

fit_mcmc_Vecchia(observed_locs, predicted_locs, observed_field)
