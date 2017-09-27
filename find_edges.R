# Load libraries
library(randomForest)
library(doMC)
registerDoMC(20L)

# Define find_edges function
find_edges <- function(dat, 
                       mtry = NULL,
                      ntree = NULL,
                      focus = NULL, 
                   directed = FALSE) {
  
  # Prelimz
  if (!is.null(focus)) {
    if (length(focus) != 1L) {
      if (is.character(focus)) {
        focus <- seq_len(p)[colnames(dat) %in% focus]
      }
      not_focus <- setdiff(seq_len(ncol(dat)), focus)
      dat <- dat[, c(focus, not_focus)]
    }
    x <- dat[, seq_len(focus)]
    p <- ncol(x)
  } else {
    p <- ncol(dat)
  }
  if (is.null(mtry)) {
    mtry <- floor(p / 3L)
  }
  if (is.null(ntree)) {
    ntree <- 2000L
  }
  
  # Define edge_wts function
  edge_wts <- function(j) {
    out <- double(length = p)
    names(out) <- colnames(dat)[seq_len(p)]
    if (!is.null(focus) && j <= p) {
      x <- x[, -j]
    } else {
      x <- dat[, -j]
    }
    y <- dat[, j]
    fit <- randomForest(x, y, mtry = mtry, ntree = ntree, importance = TRUE)
    imp <- importance(fit, type = 1L, scale = TRUE)
    out[rownames(imp)] <- as.numeric(imp)
    return(out) 
  }
  
  # Execute in parallel
  adj_mat <- foreach(j = seq_len(p), .combine = cbind) %dopar% edge_wts(j)
  dimnames(adj_mat) <- list(colnames(dat)[seq_len(p)], colnames(dat))
  if (is.null(focus) && !directed) {
    adj_mat <- pmax(adj_mat, t(adj_mat))   # Or average?
  }
  
  # Export
  return(adj_mat)
  
}


