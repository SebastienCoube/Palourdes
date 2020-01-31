observed_locs = matrix(rnorm(2000), ncol = 2)*100
predicted_locs = matrix(rnorm(2000), ncol = 2)*100
observed_field = t(chol(GpGp::matern_isotropic(locs =  observed_locs, c(1, 1, 1, 0))))%*%rnorm(1000) + 0.1*rnorm(1000)

fit_mcmc_Vecchia(observed_locs = observed_locs, observed_field = observed_field, predicted_locs = predicted_locs, pc_prior_range = c(1, .5))
