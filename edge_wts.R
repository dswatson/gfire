# Load libraries
library(ranger)
library(doMC)
registerDoMC(20L)

# Define edge_wts function
edge_wts <- function(dat, 
                     candidates = NULL,
                           mtry = round(p / 3L),
                          ntree = 2000L,
                         lambda = 0L,
                       directed = TRUE) {
  
  # Prelimz
  p <- ncol(dat)
  if (is.null(colnames(dat))) {
    colnames(dat) <- paste0('x', seq_len(p))
  }
  if (!is.null(candidates)) {
    if (is.character(candidates)) {
      candidates <- seq_len(p)[colnames(dat) %in% candidates]
    }
    others <- setdiff(seq_len(p), candidates)
    dat <- dat[, c(candidates, others)]
    x <- dat[, seq_len(candidates)]
    p <- ncol(x)
    p_names <- colnames(x)
  } else {
    p_names <- colnames(dat)
  }
  
  # Define rf_fit function
  rf_fit <- function(j) {
    out <- double(length = p)
    names(out) <- p_names
    y_name <- colnames(dat)[j]
    if (!is.null(candidates)) {
      if (j <= p) {
        dat <- x
      } else {
        dat <- cbind(x, dat[, j])
      }
    }
    fit <- ranger(data = dat, dependent.variable.name = y_name,
                  num.trees = ntree, mtry = mtry, 
                  importance = 'permutation', write.forest = FALSE,
                  scale.permutation.importance = TRUE,
                  num.threads = 1L, verbose = FALSE)
    # Use split.select.weights arg for prior info
    imp <- importance(fit)
    keep <- as.logical(imp > lambda)
    out[keep] <- imp[keep]
    return(out) 
  }
  
  # Execute in parallel
  adj_mat <- foreach(j = seq_len(p), .combine = cbind) %dopar% rf_fit(j)
  dimnames(adj_mat) <- list(p_names, colnames(dat))
  if (!directed) {
    adj_mat <- (adj_mat + t(adj_mat)) / 2L
  }
  
  # Export
  return(adj_mat)
  
}


