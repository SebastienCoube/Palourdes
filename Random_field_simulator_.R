
####################
# fit_vecchia_mcmc #
####################

### Desciption # ###
# Function that fits the Nearest Neighbor Gaussian Process using 
# delayed acceptance
# chromatic sampling
# parallel chains in order to improve convergence diagnostic
# automatic kernel variance adjustment and log-parametrized covariance parameters such as range and scale
# estimating model : field ~ N (beta_0, covmat), observed_field ~ N(beta_0, covmat + noise_variance I),  observed_field/field ~ N(0, noise_variance I)
# covmat (i, j) = GpGp's parametrization of no-nugget Matérn covariance

### Arguments ###
# benchmark is a list that gathers 
# geographic locations : a 2-column matrix under the name "locs"
# the observed field at the geographic locations : a numeric vector under the name "noisy_field"
# n_chains is the number of parallel MCMC chains
# m is the number of nearest neighbors in Vecchia approximation
# convergence_break is a vector of two stopping criterion for Gelman-Rubin-Brooks test
# first criterion : multivariate, second : univariate
# if only one criterion is met, the chain breaks
# applies to high-level parameters (covariance, intercept)
# pc_prior_scale, pc_prior_range : PC prior parameters
# see for example Fuglstad and al. https://arxiv.org/pdf/1503.00256.pdf
# if NULL, default guesses are picked from data
# n_cores : number of cores used to compute parallel chains
# if NULL : set equal to n_chains
# n_join is the number of iterations after which Gelman-Rubin-Brooks diagnostic is computed
# field_thinning is the proportion of field samples that are kept, it can save a lot of RAM...
# please keep it a divisor of n_join
# burn_in : the proportion of observations discarded to compute Gelman-Rubin-Brooks diagnostic
# n_delayed_acceptance : the number of observations that are used to do delayed acceptance
# if NULL it is set to the fifth of the number of observations, according to our experiments
# if you really want to touch it, keep it between m+1 and the number of observations...

### Outputs ###
# A list containing 3 objects
# chains : the result of the MCMC run. It contains n_chains lists, each corresponding to one chain. Each sub-list has the following bjects :
# iteration is a 2-colums matrix that records the iteration at the end of each chains join and the associated CPU time
# transition_kernel_sd is a list that stocks the (current) automatically-tuned transition kernels standard deviations
# params is a list that stocks the (current) parameters of the model, including covariance parameters, the value of the sampled field, etc
# records is a list that stocks the recorded parameters of the model, including covariance parameters, the value of the sampled field, etc. In terms of RAM, those are the biggest bit !
# estimates is the model's answer. It contains : 
# Information about the marginal distribution of each observation of the Gaussian field
# Information about the marginal distribution of the high-level parameters, with 2 parametrization
# INLA's parametrization
# log-parametrization of covariance parameters. In GpGp's parametrization, log(\sigma), log(\alpha), log(\sigma \tau)

fit_mcmc_Vecchia = function(observed_locs, predicted_locs, observed_field, 
                            n_iterations = 40000, n_chains = 3, m = 10, convergence_break = c(1.1, 1.1), 
                            pc_prior_range = NULL, pc_prior_sd = NULL,
                            n_cores = NULL, n_join = 800, field_thinning = .05, 
                            burn_in = .5, n_delayed_acceptance = NULL
)
{
  
  #######################
  # SOME INITIALIZATION #
  #######################
  
  # time  
  t_begin = Sys.time()
  # seed
  set.seed(1)
  # cleaning RAM
  gc()
  # ordering locations using max-min ordering. according to https://arxiv.org/pdf/1609.05372.pdf this works better
  locs_sample_obs = GpGp::order_maxmin(observed_locs)
  locs_sample_pred = GpGp::order_maxmin(predicted_locs)
  locs = rbind(observed_locs[locs_sample_obs,], predicted_locs[locs_sample_pred,])
  observed_field = c(observed_field)[locs_sample_obs]
  # just creating a variable of the umber of observations
  n_obs = length(observed_field)
  n_simulated = nrow(locs)
  # burn-in proportion, when  i iterations are done the (1 - burn_in) * i  last iterations are used
  if(is.null(n_delayed_acceptance)) n_delayed_acceptance = round(n_simulated/5) #number of observations used to "taste" Vecchia approximation
  # Computation parameters
  if(is.null(n_cores)) n_cores = n_chains #number of cores for parallel computing of chains
  # PC prior parameters
  if(is.null(pc_prior_range))pc_prior_range = c(30, .5)
  if(is.null(pc_prior_sd))pc_prior_sd = c(sqrt(var(observed_field)/2), .1)
  
  ###############################
  # Vecchia approximation setup #
  ###############################
  
  # This object gathers the NNarray table used by GpGp package and related objects
  
  vecchia_approx = list()
  #extracting NNarray =  nearest neighbours for Vecchia approximation
  vecchia_approx$NNarray =  GpGp::find_ordered_nn(locs, m)
  #computations from vecchia_approx$NNarray in order to create sparse Cholesky using Matrix::sparseMatrix
  #non_NA indices from vecchia_approx$NNarray
  vecchia_approx$NNarray_non_NA = !is.na(vecchia_approx$NNarray)
  #column idx of the uncompressed sparse Cholesky factor
  vecchia_approx$sparse_chol_column_idx = vecchia_approx$NNarray[vecchia_approx$NNarray_non_NA]
  #row idx of the uncompressed sparse Cholesky factor
  vecchia_approx$sparse_chol_row_idx = row(vecchia_approx$NNarray)[vecchia_approx$NNarray_non_NA]
  #color the graph iduced by Vecchia approx
  vecchia_approx$coloring = naive_greedy_coloring(vecchia_approx$NNarray)$cols
  
  
  #########################
  # Parallel chains setup #
  #########################
  
  
  chains = list()
  
  # creating sub-list : each sub-list is a chain
  for(i in seq(n_chains))
  {
    chains[[paste("chain", i, sep = "_")]] = list()
  }
  
  # for each chain, creating sub-lists in order to stock all the stuff that is related to one chain, including : 
  # iteration is a 2-colums matrix that records the iteration at the end of each chains join and the associated CPU time
  # transition_kernel_sd is a list that stocks the (current) automatically-tuned transition kernels standard deviations
  # params is a list that stocks the (current) parameters of the model, including covariance parameters, the value of the sampled field, etc
  # records is a list that stocks the recorded parameters of the model, including covariance parameters, the value of the sampled field, etc. In terms of RAM, those are the biggest bit !
  for(i in seq(n_chains))
  {
    
    chains[[i]]$iterations =c(0, 0)
    names(chains[[i]]$iterations) = c("iteration", "time")
    # Starting points for transition kernels, will be adaptatively tuned
    chains[[i]]$idx = i
    chains[[i]]$transition_kernel_sd = list()
    chains[[i]]$transition_kernel_sd$covariance_params = .05
    chains[[i]]$transition_kernel_sd$log_noise_variance = .1
    chains[[i]]$transition_kernel_sd$beta_0 = var(observed_field)
    
    #starting points for covariance parameters
    chains[[i]]$params$beta_0 = mean(observed_field)+rnorm(1, 0, sqrt(var(observed_field)/sqrt(n_obs)))
    chains[[i]]$params$log_scale = sample(log(var(observed_field))-log(seq(1, 50, 1)), 1)
    chains[[i]]$params$log_noise_variance = sample(log(var(observed_field))-log(seq(1, 50, 1)), 1)
    chains[[i]]$params$log_range = sample(log(max(dist(locs[1:100,])))-log(seq(10, 100, 10)), 1)
    #the field will be smoothed in order to match the randomized covariance parameters
    chains[[i]]$params$field  = c(observed_field, rep(0, nrow(predicted_locs)))
    
    # storing chain results
    chains[[i]]$records$params$beta_0 = rep(0, n_iterations)
    chains[[i]]$records$params$log_scale = rep(0, n_iterations)
    chains[[i]]$records$params$log_noise_variance = rep(0, n_iterations)
    chains[[i]]$records$params$log_range = rep(0, n_iterations)
    #matrix instead of vector since there is n field parameters = better thin the field
    chains[[i]]$records$params$field = matrix(0, nrow = round(n_iterations*field_thinning), ncol = n_simulated)
  }
  

  # End of initializing
  for(i in seq(n_chains) ) chains[[i]]$iterations[2] = Sys.time()- t_begin
  print(paste("Setup done,", as.numeric(Sys.time()- t_begin, units = "secs"), "s elapsed" ))
  
  
  
  iter = 1
  
  #################
  # MCMC SAMPLING #
  #################
  
  while(iter<n_iterations)
  {
    print(iter)
    # determine how many iterations will go on
    n_iterations_update = min(iter + n_join, n_iterations+1)-iter
    # export informations on number of update to the chains
    #update chains
    mcmc_update = parallel::mclapply(mc.cores = n_chains, 
        X = lapply(chains, function(chain)list(idx = chain$idx, params = chain$params, transition_kernel_sd = chain$transition_kernel_sd)), function(chain)
    {
      ####################
      # Initialize stuff #
      ####################
      
      # set seed 
      set.seed(100*iter+chain$idx)
      # copy iteration variable
      iter_begin = iter
      
      ##################
      # Vecchia factor #
      ##################
      # Initializing compressed inverse Cholesky using GpGp package. This form is mostly used in parameter updating
      compressed_sparse_chol = GpGp::vecchia_Linv(c(1, exp(chain$params$log_range), 1, 0), "matern_isotropic", locs, vecchia_approx$NNarray)
      # From compressed inverse Cholesky, the diagonal of the approximated precision matrix is computed. It is used in field updating
      sparse_precision_diag = (compressed_sparse_chol[vecchia_approx$NNarray_non_NA]^2)%*%Matrix::sparseMatrix(i = seq(length(vecchia_approx$sparse_chol_column_idx)), j = vecchia_approx$sparse_chol_column_idx, x = rep(1, length(vecchia_approx$sparse_chol_row_idx)))
      # Uncompressed form of the inverse Cholesky. This form is mostly used in field updating
      sparse_chol = Matrix::sparseMatrix(i = vecchia_approx$sparse_chol_row_idx, j = vecchia_approx$sparse_chol_column_idx, x = compressed_sparse_chol[vecchia_approx$NNarray_non_NA])
      
      
      ###################
      # Storing results #
      ###################
      # this part re-creates a small portion of the $records objects of each chain. It fills it with chain states during the run, and then updates each chain with the new values
      records = list()
      # storing chain results
      records$params$beta_0 = rep(0, n_iterations_update)
      records$params$log_scale = rep(0, n_iterations_update)
      records$params$log_noise_variance = rep(0, n_iterations_update)
      records$params$log_range = rep(0, n_iterations_update)
      #matrix instead of vector since there is n field parameters = better thin the field
      records$params$field = matrix(0, nrow = sum(round(seq(iter_begin, iter_begin+n_iterations_update-1)*field_thinning) == (seq(iter_begin, iter_begin+n_iterations_update-1)*field_thinning)), ncol = n_simulated)
      
      # acceptance results, used to tune proposal kernels
      records$acceptance$covariance_parameters = matrix(0, n_iterations_update, 2)
      records$acceptance$log_noise_variance = rep(0, n_iterations_update)
      records$acceptance$beta_0 = rep(0, n_iterations_update)
      
      Q_AA =  Matrix::sparseMatrix(i = vecchia_approx$sparse_chol_row_idx, j = vecchia_approx$sparse_chol_column_idx, x = compressed_sparse_chol[vecchia_approx$NNarray_non_NA])
      Q_AA = exp(-chain$params$log_scale)*Matrix::crossprod(Q_AA)
      Matrix::diag(Q_AA) = Matrix::diag(Q_AA) + exp(-chain$params$log_noise_variance)
      
      cond_mean = chain$params$beta_0 + Matrix::solve(Q_AA, 
                                                          exp(-chain$params$log_noise_variance)*
                                                            (c(observed_field-chain$params$beta_0, rep(0, n_simulated-n_obs))))
      chain$params$field = as.vector( cond_mean + Matrix::solve(Matrix::chol(Q_AA), rnorm(n_simulated)))
      
      #################
      # Gibbs sampler #
      #################
      # The Gibbs sampler is run for a given number of iterations
      for(iter in seq(1, n_iterations_update))
      {

        ###########################################
        # covariance parameters : scale and range #
        ###########################################
        #scale and range are block updated
        # adaptative random walk kernel if iteration is smaller than 3000
        if((iter_begin + iter <5000) & ((iter_begin + iter) / 200 == ((iter_begin+iter) %/% 200)))
        {
          if(mean(records$acceptance$covariance_parameters[seq(iter-199, iter), 2]) < .2) chain$transition_kernel_sd$covariance_params = chain$transition_kernel_sd$covariance_params * .8
          if(mean(records$acceptance$covariance_parameters[seq(iter-199, iter), 2]) > .3) chain$transition_kernel_sd$covariance_params = chain$transition_kernel_sd$covariance_params / .8
        }
        # proposing new values
        innovation = rnorm(2, 0, chain$transition_kernel_sd$covariance_params)
        new_log_range= chain$params$log_range+innovation[1]
        new_log_scale= chain$params$log_scale+innovation[2]
        new_compressed_sparse_chol  = GpGp::vecchia_Linv(c(1, exp(new_log_range), 1, 0), "matern_isotropic", locs, vecchia_approx$NNarray)
        new_sparse_chol = Matrix::sparseMatrix(i = vecchia_approx$sparse_chol_row_idx, j = vecchia_approx$sparse_chol_column_idx, x = compressed_sparse_chol[vecchia_approx$NNarray_non_NA])
        
        Q_AA =  Matrix::sparseMatrix(i = vecchia_approx$sparse_chol_row_idx, j = vecchia_approx$sparse_chol_column_idx, x = compressed_sparse_chol[vecchia_approx$NNarray_non_NA])
        Q_AA = exp(-chain$params$log_scale)*Matrix::crossprod(Q_AA)
        Matrix::diag(Q_AA) = Matrix::diag(Q_AA) + exp(-chain$params$log_noise_variance)
        
        cond_mean = chain$params$beta_0 + Matrix::solve(Q_AA, 
                                                        exp(-chain$params$log_noise_variance)*
                                                          (c(observed_field-chain$params$beta_0, rep(0, n_simulated-n_obs))))
        
        new_Q_AA =  Matrix::sparseMatrix(i = vecchia_approx$sparse_chol_row_idx, j = vecchia_approx$sparse_chol_column_idx, x = new_compressed_sparse_chol[vecchia_approx$NNarray_non_NA])
        new_Q_AA = exp(-new_log_scale)*Matrix::crossprod(new_Q_AA)
        Matrix::diag(new_Q_AA) = Matrix::diag(new_Q_AA) + exp(-chain$params$log_noise_variance)
          
        new_cond_mean = chain$params$beta_0 + Matrix::solve(new_Q_AA, 
                                                        exp(-chain$params$log_noise_variance)*
                                                          (c(observed_field-chain$params$beta_0, rep(0, n_simulated-n_obs)))
                                                        )
        new_field = as.vector( new_cond_mean + Matrix::solve(Matrix::chol(new_Q_AA), rnorm(n_simulated)))
        
        transition_ratio = 
          +.5*Matrix::determinant(new_Q_AA)$modulus -.5*sum(as.vector(Matrix::crossprod(new_field-new_cond_mean,new_Q_AA)) *(new_field-new_cond_mean))
          -.5*Matrix::determinant(Q_AA)$modulus     +.5*sum(as.vector(Matrix::crossprod(chain$params$field-cond_mean,Q_AA))*(chain$params$field-cond_mean))
        GP_ratio = 
          ll_compressed_sparse_chol(Linv= new_compressed_sparse_chol, log_scale = new_log_scale,          field = new_field - chain$params$beta_0         , NNarray =  vecchia_approx$NNarray) - 
          ll_compressed_sparse_chol(Linv= compressed_sparse_chol,     log_scale = chain$params$log_scale, field = chain$params$field - chain$params$beta_0, NNarray =  vecchia_approx$NNarray)
        field_response_ratio = 
          sum(dnorm(x = observed_field, mean = chain$params$field[seq(n_obs)], sd = exp(0.5*chain$params$log_noise_variance), log = T)-
              dnorm(x = observed_field, mean =          new_field[seq(n_obs)], sd = exp(0.5*chain$params$log_noise_variance), log = T))
        pc_prior_ratio = (new_log_range-chain$params$log_range)*(-ncol(locs)/2-1) +
          log(pc_prior_range[2])*pc_prior_range[1]^(ncol(locs)/2)*sqrt(8)*(exp(new_log_range)^(-ncol(locs)/2)-exp(chain$params$log_range)^(-ncol(locs)/2)) +
          log(pc_prior_sd[2])/pc_prior_sd[1]*(exp(new_log_scale/2)-exp(chain$params$log_scale/2))
        print("GP")
        print(GP_ratio)
        print("pc")
        print(pc_prior_ratio)
        print("response")
        print(field_response_ratio)
        print("transition")
        print(transition_ratio)
        print("sum")
        print(field_response_ratio+GP_ratio+pc_prior_ratio+transition_ratio)
        if(field_response_ratio+GP_ratio+pc_prior_ratio+transition_ratio > log(runif(1)))        
        {
          records$acceptance$covariance_parameters[iter,2] = 1
          #parameter updating
          chain$params$log_range = new_log_range
          chain$params$log_scale = new_log_scale
          chain$params$field = new_field
          #updating Vecchia cholesky
          compressed_sparse_chol = new_compressed_sparse_chol
          Q_AA = new_Q_AA
          cond_mean = new_cond_mean
          #Used in field updating
          sparse_chol = Matrix::sparseMatrix(i = vecchia_approx$sparse_chol_row_idx, j = vecchia_approx$sparse_chol_column_idx, x = compressed_sparse_chol[vecchia_approx$NNarray_non_NA])
        }
        
        ##############
        # Field mean #
        ##############
        # adaptative random walk kernel
        if((iter_begin + iter <5000) & ((iter_begin + iter) / 200 == ((iter_begin+iter) %/% 200)))
        {
          if(mean(records$acceptance$beta_0[seq(iter-199, iter)]) < .2) chain$transition_kernel_sd$beta_0 = chain$transition_kernel_sd$beta_0 * .8
          if(mean(records$acceptance$beta_0[seq(iter-199, iter)]) > .3) chain$transition_kernel_sd$beta_0 = chain$transition_kernel_sd$beta_0 / .8
        }
        # proposing a new parameter
        innovation = rnorm(1, 0, chain$transition_kernel_sd$beta_0)
        ll_ratio = 
          # ll with new parameters
          ll_compressed_sparse_chol(Linv= compressed_sparse_chol, log_scale = chain$params$log_scale, field = chain$params$field - chain$params$beta_0 - innovation, NNarray = vecchia_approx$NNarray)-
          # ll with old parameters
          ll_compressed_sparse_chol(Linv= compressed_sparse_chol, log_scale = chain$params$log_scale, field = chain$params$field - chain$params$beta_0             , NNarray = vecchia_approx$NNarray)
        # acceptance step
        if(ll_ratio >log(runif(1)))
        {
          records$acceptance$beta_0[iter] = 1
          chain$params$beta_0 = chain$params$beta_0 + innovation
        }    
        
        ##################
        # Noise variance #
        ##################
        # adaptative random walk kernel
        if((iter_begin + iter <5000) & ((iter_begin + iter) / 200 == ((iter_begin+iter) %/% 200)))
        {
          if(mean(records$acceptance$log_noise_variance[seq(iter-199, iter)]) < .2) chain$transition_kernel_sd$log_noise_variance = chain$transition_kernel_sd$log_noise_variance * .8
          if(mean(records$acceptance$log_noise_variance[seq(iter-199, iter)]) > .3) chain$transition_kernel_sd$log_noise_variance = chain$transition_kernel_sd$log_noise_variance / .8
        }
        # proposing a new parameter
        innovation = rnorm(1, 0, chain$transition_kernel_sd$log_noise_variance)
        # acceptance step
        #ll ratio : observed signal ~  n_obs(mu = simulated signal, sigma = noise variance * In)
        if(sum(dnorm(x = observed_field, mean = chain$params$field[seq(n_obs)], sd = exp(0.5*(chain$params$log_noise_variance + innovation)), log = T)-
               dnorm(x = observed_field, mean = chain$params$field[seq(n_obs)], sd = exp(0.5*chain$params$log_noise_variance)               , log = T))
           >log(runif(1)))
        {
          records$acceptance$log_noise_variance[iter] = 1
          chain$params$log_noise_variance = chain$params$log_noise_variance+innovation
        }

        ######################
        # Saving chain state #
        ######################
        
        records$params$log_scale[iter] = chain$params$log_scale
        records$params$beta_0[iter] = chain$params$beta_0
        records$params$log_noise_variance[iter] = chain$params$log_noise_variance
        records$params$log_range[iter] = chain$params$log_range
        if(round((iter_begin+iter)*field_thinning) == ((iter_begin+iter)*field_thinning)) records$params$field[match(0, records$params$field[,1])[1],] = chain$params$field - chain$params$beta_0
      }
      chain$records = records
      return(chain)
    })
    
    #################################
    # End of 1 Gibbs sampler update #
    #################################
    
    
    ###################
    # Updating chains #
    ###################
    
    
    for(i in seq(n_chains))
    {
      chains[[i]]$transition_kernel_sd = mcmc_update[[i]]$transition_kernel_sd
      chains[[i]]$params = mcmc_update[[i]]$params
      for(a in ls(chains[[i]]$records))
      {
        for(b in ls(chains[[i]][[a]]))
        {
          if(is.vector(chains[[i]][["records"]][[a]][[b]]))
          {
            chains[[i]][["records"]][[a]][[b]][seq(iter, iter+n_iterations_update-1)] = mcmc_update[[i]][["records"]][[a]][[b]]
          }
          if(is.matrix(chains[[i]][["records"]][[a]][[b]])&(b!="field"))
          {
            chains[[i]][["records"]][[a]][[b]][seq(iter, iter+n_iterations_update-1),] = mcmc_update[[i]][["records"]][[a]][[b]]
          }
          if("field" == b)
          {
            field_idx = ((seq(n_iterations_update) + iter)*field_thinning)[sapply(seq(n_iterations_update) + iter, function(i)(round(i*field_thinning) ==( i*field_thinning)))]
            chains[[i]][["records"]][[a]][[b]][field_idx,] = mcmc_update[[i]][["records"]][[a]][[b]]
          }
        }
      }
    }
    iter = iter+n_iterations_update
    for(i in seq(n_chains) ) chains[[i]][["iterations"]]  = rbind(chains[[i]][["iterations"]], c(iter, as.numeric(Sys.time()-t_begin, units = "secs")))
    
    ##################################
    # Gelman-Brooks-Rubin diagnostic #
    ##################################
    
    # computing within and between variance matrices for higher level parameters  
    within_variance = lapply(chains, function(chain)
      var(sapply(c("log_range", "log_scale", "log_noise_variance", "beta_0"), function(name)chain$records$params[[name]])[seq(burn_in*(iter-1), iter-1),])
    )
    within_variance = Reduce("+", within_variance)/n_chains
    means = sapply(chains, function(chain)
      apply(sapply(c("log_range", "log_scale", "log_noise_variance", "beta_0"), function(name)chain$records$params[[name]])[seq(burn_in*(iter-1), iter-1),], 2, mean)
    )
    between_variance = var(t(means))
    # multivariate diagnostic
    MPSRF = (iter - 2) / (iter - 1) + (n_chains + 1) / n_chains * svd(solve(within_variance)%*%between_variance)$d[1]
    names(MPSRF) = "Multivariate"
    # univariate diagnostics
    Individual_PSRF = (iter - 2) / (iter - 1) + (n_chains + 1) / n_chains * 1/diag(within_variance)*diag(between_variance)
    # show diagnostics
    print ("Gelan-Rubin-Brooks diag : ")
    print(c(MPSRF, Individual_PSRF))
    
    #########################
    # EFFECTIVE SAMPLE SIZE #
    #########################
    
    # prints Effective Sample Size for each high-level parameter
    ESS = sapply(c("log_range", "log_scale", "log_noise_variance", "beta_0"), function(name)sapply(chains, function(chain)coda::effectiveSize(
      chain$records$params[[name]][seq(burn_in*(iter-1), iter-1)])))
    ESS = rbind(ESS, apply(ESS, 2, sum))
    row.names(ESS) = c(sapply(seq(n_chains), function(i) paste("chain", i, sep = "_")), "total")
    print("ESS")
    print(ESS)
    
    #########
    # PLOTS #
    #########
    
    # plots chains of each high-level parameter. There will be as many plots as high-level parameters, and as many curves in each plot as chains. 
    # plot starts at burn-in
    par(mfrow = c(4, 1))
    # loop over parameters
    for(name in c("log_range", "log_scale", "log_noise_variance", "beta_0"))
    {
      to_be_plotted = lapply(chains, function(chain)chain$records$params[[name]][seq(burn_in*(iter-1), iter-1)])
      plot(seq(burn_in*(iter-1), iter-1), seq(burn_in*(iter-1), iter-1), type = "n", xlab = "iteration", ylab = name, main = name, 
           ylim = c(min(unlist(to_be_plotted)), max(unlist(to_be_plotted))))
      col = 1
      # loop over chains
      for(x in to_be_plotted)
      {
        lines(seq(burn_in*(iter-1), iter-1), x, col = col)
        col = col+1
      }
    }
    
    
    #############
    # ESTIMATES #
    #############
    
    # creating an empty list to store the estimates
    estimates = list()
    # field
    estimates$field = do.call(cbind, list(
      "mean" = apply(FUN = mean, MARGIN = 2, 
                     do.call(rbind, 
                             lapply(chains, function(chain) chain$records$params$field[seq(field_thinning*burn_in*(iter-1), field_thinning*(iter-1)),])
                     )), 
      "q025" = apply(FUN = function(x)quantile(x, c(0.025)), MARGIN = 2, 
                     do.call(rbind, 
                             lapply(chains, function(chain) chain$records$params$field[seq(field_thinning*burn_in*(iter-1), field_thinning*(iter-1)),])
                     )),
      "q975" = apply(FUN = function(x)quantile(x, c(0.975)), MARGIN = 2, 
                     do.call(rbind, 
                             lapply(chains, function(chain) chain$records$params$field[seq(field_thinning*burn_in*(iter-1), field_thinning*(iter-1)),])
                     )),
      "sd" = apply(FUN = sd, MARGIN = 2, 
                   do.call(rbind, 
                           lapply(chains, function(chain) chain$records$params$field[seq(field_thinning*burn_in*(iter-1), field_thinning*(iter-1)),])
                   ))))
    #parameters
    #log-params
    range_sample = do.call(c, lapply(chains, function(chain) chain$records$params[["log_range"]][seq(burn_in*(iter-1), (iter-1))]))
    scale_sample = do.call(c, lapply(chains, function(chain) chain$records$params[["log_scale"]][seq(burn_in*(iter-1), (iter-1))]))
    noise_variance_sample = do.call(c, lapply(chains, function(chain) chain$records$params[["log_noise_variance"]][seq(burn_in*(iter-1), (iter-1))]))
    beta_0_sample = do.call(c, lapply(chains, function(chain) chain$records$params[["beta_0"]][seq(burn_in*(iter-1), (iter-1))]))
    estimates$log_params = t(sapply(list(range_sample, scale_sample, noise_variance_sample, beta_0_sample), function(x)
      c(mean(x), quantile(x, c(.025, .975)), sd(x))))
    colnames(estimates$log_params) =c("mean", "q025", "q975", "sd")
    rownames(estimates$log_params) =c("log_range", "log_scale", "log_noise_variance", "beta_0")
    #INLA parametrization
    estimates$inla_params = t(sapply(list(exp(range_sample)*sqrt(8), exp(scale_sample/2), exp(-noise_variance_sample), beta_0_sample), function(x)
      c(mean(x), quantile(x, c(.025, .975)), sd(x))))
    colnames(estimates$inla_params) =c("mean", "q025", "q975", "sd")
    rownames(estimates$inla_params) =c("Range_for_spat", "Stdev_for_spat", "Precision_for_he_Gaussian_observations", "beta_0")
    
    #print the estimates
    print("Estimates")
    print(estimates$log_params)
    cat("\n")
    cat("\n")
    
    #################
    # Saving chains #
    #################
    # just saves the chains at each 10 joins
    if(iter / n_join == 10) saveRDS(chains, "Veccchia_run_chains.RDS")
    if((MPSRF<convergence_break[1])|all(Individual_PSRF<convergence_break[2]))break
  }
  
  return(list("chains" = chains, "estimates" = estimates, "locs_sample" = list("obs" = locs_sample_obs, "pred" = locs_sample_pred)))
  
}



###############################
# Vecchia likelihood function #
###############################

#Gaussian likelihood from Vecchia linv, for a 0-mean GP
ll_compressed_sparse_chol = function(Linv, field, NNarray, log_scale)
{
  chol_field = GpGp::Linv_mult(Linv = Linv, z = field, NNarray)
  sum(log(Linv[NNarray[,1],1]))-nrow(NNarray)*0.5*log_scale-0.5*sum(chol_field^2)/exp(log_scale)
}

