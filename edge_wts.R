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
                     seed = NULL) {
  
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
    x <- dat[, seq_len(candidates)]
    p <- ncol(x)
    p_names <- colnames(x)
  } else {
    p_names <- colnames(dat)
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
    names(out) <- p_names
    if (!is.null(candidates) && j <= p) {
      x <- x[, -j]
    } else {
      x <- dat[, -j]
    }
    y <- dat[, j]
    fit <- regression_forest(x, y, mtry = mtry, num.trees = ntree, 
                             honesty = FALSE, num.threads = 1L, seed = seed) 
    imp <- v_imp(fit)
    out[colnames(x)] <- imp 
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


