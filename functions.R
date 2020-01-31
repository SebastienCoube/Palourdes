####################################
# Moral graph coloring functions #
####################################

# "greedy" coloring algorithm for bigger n_obs 
naive_greedy_coloring = function(NNarray)
{
  #number of nodes
  n_obs = nrow(NNarray)
  #creating adjacency matrix
  M =  M =  Matrix::sparseMatrix(i = NNarray[!is.na(NNarray)], j = row(NNarray)[!is.na(NNarray)], x = 1)%*%Matrix::sparseMatrix(i = row(NNarray)[!is.na(NNarray)], j = NNarray[!is.na(NNarray)], x = 1)
  M = M!=0
  #deducting degrees
  degrees = as.vector(rep(1, n_obs)%*%M)
  #getting adjacent nodes of a given node
  idx = split(M@i+1, rep(seq_along(diff(M@p)),diff(M@p)))
  #creating a color * node matrix of incompatibilities
  incompatibilities = matrix(0, n_obs+1, max(degrees))
  cols = rep(0, n_obs)
  
  for(i in seq(n_obs))
  {
    cols[i] = match(0, incompatibilities[i,])
    incompatibilities[idx[[i]],cols[i]] = 1
  }
  return(list("cols" = cols, "M" = M))
}

naive_greedy_coloring = compiler::cmpfun(naive_greedy_coloring)

##########################
# Mixed parents function #
##########################

get_mixed_NNarray = function(m, m_knots, locs, n_knots)
{
  
  if(m_knots==0)return(GpGp::find_ordered_nn(locs, m-m_knots))
  n_obs = nrow(locs)
  if(m_knots <m)NNarray = cbind(GpGp::find_ordered_nn(locs, m-m_knots), matrix(NA, n_obs, m_knots))
  if(m_knots == m)NNarray = cbind(seq(n_obs), matrix(NA, n_obs, m_knots))
  
  dists = fields::rdist(locs[seq(n_knots),], locs)
  dists[row(dists)>=col(dists)]=Inf
  dists[cbind(NNarray[NNarray<n_knots+1],row(NNarray)[NNarray<n_knots+1])]=Inf
  
  for(i in seq(m-m_knots+2, n_obs))
  {
    for(j in seq(m-m_knots+2, min(i, m+1)))
    {
      k = which.min(dists[,i])
      dists[k,i] = Inf
      NNarray[i, j] = k
    }
  }
  colnames(NNarray) = c("ind_index", sapply(seq(m-m_knots), function(i)paste("NN", i, sep = "_")), sapply(seq(m_knots), function(i)paste("knot", i, sep = "_")))
  NNarray
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

#####