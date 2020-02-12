remove(list = ls())
# number of obs
n = 200
# covariance parameters
log_range = 0
log_scale = 10
log_noise_variance = 8
beta_0 = 1
# observed locations
observed_locs = cbind(30*runif(n), 1)
# locatons to be predicted
predicted_locs = cbind(seq(0, 30, 0.01), 1)
# field simulation
true_field = beta_0 + exp(log_scale/2)*t(chol(GpGp::matern_isotropic(c(1, exp(log_range), 1, 0), observed_locs)))%*%rnorm(n) 
observed_field = true_field + sqrt(exp(log_noise_variance))*rnorm(n)
# reordering for Vecchia approx
locs_sample_obs = GpGp::order_maxmin(observed_locs)
locs_sample_pred = GpGp::order_maxmin(predicted_locs)
locs = rbind(observed_locs[locs_sample_obs,], predicted_locs[locs_sample_pred,])
observed_field = c(observed_field)[locs_sample_obs]
# just creating a variable of the number of observations
n_obs = length(observed_field)
n_simulated = nrow(locs)

# What our field does look like ? 
plot(observed_locs[locs_sample_obs,1], observed_field)
points(observed_locs[locs_sample_obs,1], true_field[locs_sample_obs], col = rep(3, n_obs))
legend("topright", legend = c("observed", "true"), fill = c(1, 3))

# Vecchia approximation graph information
vecchia_approx = list()
#extracting NNarray =  nearest neighbours for Vecchia approximation
vecchia_approx$m =  12
vecchia_approx$NNarray =  GpGp::find_ordered_nn(locs, vecchia_approx$m)
#computations from vecchia_approx$NNarray in order to create sparse Cholesky using Matrix::sparseMatrix
#non_NA indices from vecchia_approx$NNarray
vecchia_approx$NNarray_non_NA = !is.na(vecchia_approx$NNarray)
#column idx of the uncompressed sparse Cholesky factor
vecchia_approx$sparse_chol_column_idx = vecchia_approx$NNarray[vecchia_approx$NNarray_non_NA]
#row idx of the uncompressed sparse Cholesky factor
vecchia_approx$sparse_chol_row_idx = row(vecchia_approx$NNarray)[vecchia_approx$NNarray_non_NA]

# Distrbution
# Vecchia sparse cholesky factor of precision matrix for GP prior for a posteriori distribution at sampled and observed sites
compressed_sparse_chol = GpGp::vecchia_Linv(c(1, exp(log_range), 1, 0), "matern_isotropic", locs, vecchia_approx$NNarray)
# Vecchia sparse precision matrix for a posteriori distribution at sampled and observed sites
Q_AA =  Matrix::sparseMatrix(i = vecchia_approx$sparse_chol_row_idx, j = vecchia_approx$sparse_chol_column_idx, x = compressed_sparse_chol[vecchia_approx$NNarray_non_NA])
Q_AA = exp(-log_scale)*Matrix::crossprod(Q_AA)
Matrix::diag(Q_AA) = Matrix::diag(Q_AA) + c(rep(exp(-log_noise_variance), n_obs), rep(0, n_simulated-n_obs))
# a posteriori mean at sampled and observed sites
cond_mean = beta_0 - Matrix::solve(Q_AA, 
                                   - exp(-log_noise_variance)*
                                     (c(observed_field-beta_0, rep(0, n_simulated-n_obs))))
# Getting Choleky factor of a posteriori precision matrix
C0 = Matrix::Cholesky(Q_AA, perm = F)
expanded = Matrix::expand(C0)

# making one simulation of the denoised field
denoised_field = as.vector( cond_mean + Matrix::solve(Matrix::crossprod(expanded$L, expanded$P),  rnorm(n_simulated)) ) 
# making one simulation of the response field
response_field = denoised_field + rnorm(n_simulated, 0, exp(log_noise_variance/2))

# plotting at observed locations 
plot(locs[1:n_obs,1], observed_field[1:n_obs])
points(locs[1:n_obs,1], denoised_field[1:n_obs], col = rep(2, n_obs))
points(locs[1:n_obs,1], true_field[locs_sample_obs], col = rep(3, n_obs))
legend("topright", legend = c("observed", "true"), fill = c(1, 3))


# plotting at predicted locations + observed locations
plot(locs[1:n_obs,1], observed_field[1:n_obs])
points(locs[,1], denoised_field, col = rep(2, n_obs))
points(locs[1:n_obs,1], true_field[locs_sample_obs], col = rep(3, n_obs))
legend("topright", legend = c("observed", "true"), fill = c(1, 3))


# plotting response variable
plot(locs[,1], response_field)
