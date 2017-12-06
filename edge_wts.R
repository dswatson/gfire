# Define variable importance function
v_imp <- function(grf) {
  require(grf)
  split.freq <- split_frequencies(grf, 4L)
  split.freq <- split.freq / pmax(1L, rowSums(split.freq))
  weight <- seq_len(nrow(split.freq)) ^ -2L
  var.importance <- t(split.freq) %*% weight / sum(weight)
  var.importance <- as.numeric(var.importance)
  return(var.importance)
}

# Define edge_wts function
edge_wts <- function(dat, 
                     candidates = NULL,
                     mtry = NULL,
                     ntree = NULL,
                     directed = TRUE,
                     seed = NULL,
                     parallel = TRUE) {
  
  # Prelimz
  require(grf)
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
    x <- dat[, candidates]
    p <- ncol(x)
  } 
  if (is.null(mtry)) {
    mtry <- ceiling(2L * p / 3L)
  }
  if (is.null(ntree)) {
    ntree <- 2000L
  }
  
  # Define rf_fit function
  rf_fit <- function(j) {
    out <- double(length = p)
    names(out) <- colnames(dat)[seq_len(p)]
    if (!is.null(candidates) && j <= p) {
      x <- x[, -j]
    } else if (is.null(candidates)) {
      x <- dat[, -j]
    }
    y <- dat[, j]
    fit <- regression_forest(x, y, mtry = mtry, num.trees = ntree, 
                             honesty = FALSE, num.threads = 1L, seed = seed)
    out[colnames(x)] <- v_imp(fit) 
    return(out) 
  }
  
  # Execute in parallel
  if (parallel) {
    adj_mat <- foreach(j = seq_len(ncol(dat)), .combine = cbind) %dopar% rf_fit(j)
  } else {
    adj_mat <- sapply(seq_len(ncol(dat)), rf_fit)
  }
  dimnames(adj_mat) <- list(colnames(dat)[seq_len(p)], colnames(dat))
  if (!directed) {
    adj_mat <- pmax(adj_mat, t(adj_mat))
  }
  
  # Export
  return(adj_mat)
  
}


