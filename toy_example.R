observed_locs = matrix(rnorm(800), ncol = 2)*10
predicted_locs = matrix(rnorm(800), ncol = 2)*10

# low variance
observed_field = t(chol(GpGp::matern_isotropic(locs =  observed_locs, c(1, 1, 1, 0))))%*%rnorm(400) + sqrt(0.1)*rnorm(400)
test = fit_mcmc_Vecchia(observed_locs = observed_locs, observed_field = observed_field, predicted_locs = predicted_locs, pc_prior_range = c(1, .5))

# medium variance
observed_field = t(chol(GpGp::matern_isotropic(locs =  observed_locs, c(1, 1, 1, 0))))%*%rnorm(400) + rnorm(400)
test = fit_mcmc_Vecchia(observed_locs = observed_locs, observed_field = observed_field, predicted_locs = predicted_locs, pc_prior_range = c(1, .5))

# high variance
observed_field = t(chol(GpGp::matern_isotropic(locs =  observed_locs, c(1, 1, 1, 0))))%*%rnorm(400) + sqrt(5)*rnorm(400)
fit_mcmc_Vecchia(observed_locs = observed_locs, observed_field = observed_field, predicted_locs = predicted_locs, pc_prior_range = c(1, .5))


