# Load needed packages

library(Rcpp)
library(foreach)
library(doSNOW)
library(smacof)
library(tidyverse)

# Simulating & Testing ----------------------------------------------------------


# Function: generate matrix with ndim cues and n items
#'   Input: @ndims number of latent dimensions/cues
#'          @n     number of itmes/rows
#'          @min   minimum possible cue value for the uniform distr.
#'          @max   maximum possible cue value for the uniform distr.
#   Output: n x ndim matrix

gen_lat_dim <- function(ndims=2, n=16, min=-1, max= 1){
  
  # Create matrix with cue values for ndims cues
  latent_dims <- matrix(runif(ndims*n,min,max),nrow=n,ncol=ndims)
  
  # Return results
  return(latent_dims)
  
}


# Function: compute pairwise distances between two vectors
#'   Input: @mat   matrix of latent dimensions
#'          @r     r parameter of minkowski space
#'          @diag  return diagonal matrix values
#'          @upper return upper matrix half
#   Output: matrix of pairwise distances


# Main function in Rcpp:
cppFunction('NumericMatrix pair_dist_cpp(NumericMatrix x, double r) {

  int nrow = x.nrow();
  int ndim = x.ncol();
  NumericMatrix out(nrow,nrow);


  for(int i = 0; i < nrow; i ++){

    for( int j = 0; j < nrow; j++){

      double t_p = 0;

      for(int k = 0; k < ndim; k++){

       t_p += pow(x(i,k)-x(j,k),r);

      }

      out(i,j) = pow(t_p,1/r);

    }

  }

  return out;
}')


# R Wrapper:
compute_dists <- function(mat,r=2,diag=F,upper=F){
  
  result <- pair_dist_cpp(mat,r)
  
  if(!diag){
    diag(result) <- NA
  }
  
  if(!upper){
    result[upper.tri(result)] <- NA
  }
  
  
  return(result)
  
}


# Function: generate mds data
#'   Input: @ndims number of latent dimensions/cues
#'          @n     number of itmes/rows
#'          @min   minimum possible cue value for the uniform distr.
#'          @max   maximum possible cue value for the uniform distr
#'          @r     r parameter of minkowski space
#'          @diag  return diagonal matrix values
#'          @upper return upper matrix half
#'          @return_latent should the latent matrix of underlying cues and cue values be returned?
#   Output: pairwise distance matrix (if return_latent=T list with both)


gen_data_mds <- function(ndims=2, n=16, min=-1, max= 1, r=2, diag=F, upper=F, return_latent=T){
  
  # Create matrix with cue values for ndims cues
  latent_dims <- gen_lat_dim(ndims, n, min, max)
  
  # Make pair-wise distance matrix out of this
  dist_mat    <- compute_dists(latent_dims,r,diag,upper)
  
  
  if(return_latent){
    res_list <- list("latent" = latent_dims,"dist_mat" = dist_mat)
    return(res_list)
  } else {
    return(dist_mat)
  }
  
}

# Function: run MDS cross-validation for one specific number of dimensions reps times
#'   Input: @x       pairwise distance matrix
#'          @ndim    number of dimensions
#'          @reps    number of repetitions
#'          @NA_prob probability of turning an entry of x to NA
#   Output: matrix of pairwise distances


# Cross-Validation     ----------------------------------------------------------


cross_mds <- function(x,reps=500,ndim=2,NA_prob=0.2,ncores = 4,criterion = "r"){
  

  
  # create empty results vector
  cor_reps <- vector("numeric",reps)
  
  # get dimensions of x
  dims <- dim(x)
  
  # include this for foreach
  RMSE <- function(x,y){
    sqrt(mean((x-y)^2))
  }
  
  # loop over repetitions

  # initialize clusters/workers
  cl       <- makeSOCKcluster(ncores)
  registerDoSNOW(cl)
  
  # loop over reps per dimension
  cor_reps <-  foreach(rep            = 1:reps,
                       .packages      = c("smacof","tidyverse"),
                       .combine       = c,
                       .errorhandling = "pass",
                       .verbose=F) %dopar% {
                         
                         # copy x
                         temp_mat <- x
                         
                         true <- c()
                         rows <- c()
                         cols <- c()
                         
                         for(col in 1:(dims[1]-1)){
                           for(row in (col+1):(dims[2])){
                             
                             # check indices
                             # print(paste0("[",row,",",col,"]"))
                             
                             if(runif(1)<NA_prob){
                               
                               true <- c(true,temp_mat[row,col])
                               rows <- c(rows,row)
                               cols <- c(cols,col)
                               temp_mat[row,col] <- NA
                               
                             }
                             
                           }
                         }
                         
                         mds_res  <- mds(temp_mat,ndim=ndim,"ordinal")
                         pred_mat <- mds_res$conf %>%  dist(.,method = "euclidean") |>  as.matrix()
                         
                         pred <- vector(mode="numeric",length=length(true))
                         
                         for(i in 1:length(true)){
                           
                           pred[i] <- pred_mat[rows[i],cols[i]]
                           
                         }
                         
                         #cor_reps[rep] =
                         if(criterion == "r"){
                           cor(true,pred)
                         } else {
                           RMSE(true,pred)
                         }
                         
                         #
                         
                       }
  
  stopCluster(cl)
  
  return(cor_reps)
}


# Quick performance check
# bench::mark("dist" =  mds_res$conf %>%  dist(.,method = "euclidean"),
#             "cpp"  =  mds_res$conf %>%  pair_dist_cpp(.,r = 2),
#             check  = F)

# Function: run MDS cross-validation over a number of dimensions
#'   Input: @x       pairwise distance matrix
#'          @max_dim maximum number of dimensions to test
#'          @reps    number of repetitions
#'          @NA_prob probability of turning an entry of x to NA
#'          @plot    return plot
#'          @verbose return descriptives
#   Output: matrix of pairwise distances


compute_cross_mds <- function(x,max_dim=5,reps=500,NA_prob=0.2,plot=T,verbose=T,ncores = 4,criterion = "r"){
  
  
  res <- data.frame("Dim"  = rep(0,max_dim),
                    "mean" = rep(0,max_dim),
                    "sd"   = rep(0,max_dim))
  
  for(i in 1:max_dim){
    
    temp <- cross_mds(x = x,reps=reps,ndim=i,NA_prob=NA_prob,ncores=ncores,criterion=criterion)
    
    res[i,] <- c(i,mean(temp),sd(temp))
    
    if(verbose){
      print(paste0("# Dim: ", i))
    }
    
  }
  
  
  if(plot){
    p <- ggplot(res,aes(x = Dim,y = mean)) +
          geom_line(lwd=1) +
          geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=0.1) +
          geom_point(size=3,shape=21,fill="grey") +
          labs(x = "Number of Dimensions",
               y = paste0("Average ",criterion," between training- and test-set"),
               caption = "Mean and +-1 standard deviation") +
          theme_bw() +
          theme(text = element_text(size=14))
    
    return(list("CV"=res,"plot"=p))
    
  } else {
    
    return(res)
    
  }
  
}



# BIC as in Lee (2001) ----------------------------------------------------------

# Function: compute BIC according to Lee 2001
#'   Input: @s       sample estimate of standard deviation
#'          @rss     residual sum-of-squares of current solution
#'          @n_dims  number of dimentions number
#'          @n_items number of items
#   Output: BIC value


BIC_mds <- function(s,rss,n_dims,n_items){
  
  BIC = 1/s^2 * rss + n_dims*n_items*log(n_items*(n_items-1)/2)
  
  return(BIC)
}


# Function: compute sigma^2 according to Lee (2001, p. 155)
#'   Input: @sds     vector of standard deviations for each matrix cell (pair)
#'          @nstim   number of stimuli
#   Output: sigma^2 estimate

s_lee <- function(sds,nstim){
  
  s <- sum(sds)/(nstim*(nstim-1)/2)
  
  return(s)
}




# Function: compute BIC-values according to Lee 2001
#'   Input: @d_list  list of pair-wise distance matrices (upper + diag should be NA)
#'          @max_dim maximum number of dimensions to test
#'          @plot    return plot
#   Output: BIC values


compute_BIC_mds <- function(d_list,max_dim=5,plot=T){
  
  
  check_NA <- sapply(d_list,function(x){
    
    d_T <- all(is.na(diag(x)))
    u_T <- all(is.na(x[upper.tri(x)]))
    
    all(d_T,u_T)
    
  }) %>% all()
  
  if(!check_NA){
    stop("Pair-wise distance matrix has non-NAs in upper triangle or diagonal")
  }
  
  # get number of items and persons
  n_items   <- sqrt(length(d_list[[1]]))
  n_persons <- length(d_list)
  
  # Calculate s
  d_by_ID <- sapply(X = d_list, as.vector, simplify=T) %>%
    na.omit() %>%
    as.data.frame()
  
  # Sum over row-wise SDs
  sum_sd_d <- sum(apply(d_by_ID,1,sd))
  
  # Compute s
  s <- 1/(n_items*(n_items-1)/2)*sum_sd_d
  #print(s)
  
  # compute aggregate pair-wise distance matrix is given
  d_aggr <- apply(simplify2array(d_list), 1:2, mean)
  
  
  # init empty results data.frame
  res <- data.frame("Dim"    = rep(0,max_dim),
                    "Stress" = rep(0,max_dim),
                    "BIC"    = rep(0,max_dim),
                    "P"      = rep(0,max_dim),
                    "s"      = rep(s,max_dim))
  
  # for i to max_dims
  for(i in 1:max_dim){
    # save stress as well, since it does not cost much computation
    temp <- mds(d_aggr,ndim=i,type="ordinal")
    
    # get stress
    stress <- temp$stress
    
    # get residual sum of squares of MDS soluation with i dims
    rss    <-temp$rss
    
    # calculate BIC
    BIC    <- BIC_mds(s,rss,i,n_items)
    
    # safe
    res[i,1:4] <- c(i,stress,BIC,rss,i*n_items)
    
  }
  
  
  if(plot){
    p <- ggplot(res,aes(x = Dim,y = log(BIC))) +
      geom_line(lwd=1) +
      geom_point(size=3,shape=21,fill="grey") +
      labs(x = "Number of Dimensions",
           y = "log(Bayesian Information Criterion)") +
      theme_bw() +
      theme(text = element_text(size=14))
    return(list("BIC"=res,"plot"=p))
  } else {
    return(res)
  }
  
}


# Helper functions     ----------------------------------------------------------


# Fill matrix
# Transform matrix with only lower triangle entries to full matrix with 0 on the diagonals

fill_mat <- function(lower_mat){
  
  diag(lower_mat) = 0
  lower_mat[upper.tri(lower_mat)] <- t(lower_mat)[upper.tri(lower_mat)]
  
  return(lower_mat)
  
}


