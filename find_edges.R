# Load libraries
library(data.table)
library(randomForest)
library(stabs)
library(doMC)
registerDoMC(20L)

# Define find_edges function
find_edges <- function(dat, 
                       focus = NULL, 
                      sparse = FALSE, 
                       alpha = 0.05, 
                    directed = FALSE, ...) {
  
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
  if (sparse) {
    EV <- floor(alpha * (p - 1L))
    rf <- function(x, y, q, ...) {  # ntree? mtry?
      fit <- randomForest(x, y, importance = TRUE, ...)
      imp <- importance(fit, type = 1L, scale = TRUE) 
      coef_ids <- order(imp[, 1L], decreasing = TRUE)[seq_len(q)]
      ret <- ifelse(seq_len(ncol(x)) %in% coef_ids, TRUE, FALSE)
      names(ret) <- colnames(x)
      return(list(selected = ret))
    }
  }
  
  # Define edge_wts function
  edge_wts <- function(j) {
    out <- double(length = ncol(dat))
    names(out) <- colnames(dat)
    if (!is.null(focus) && j <= p) {
      x <- x[, -j]
    } else {
      x <- dat[, -j]
    }
    y <- dat[, j]
    if (sparse) {
      ss <- stabsel(x, y, fitfun = rf, cutoff = 0.75, PFER = EV, 
                    sampling.type = 'SS', papply = lapply)
      edges <- ss$selected
      if (length(edges) > 0L) {
        fit <- randomForest(x, y, importance = TRUE)
        imp <- importance(fit, type = 1L, scale = TRUE)[names(edges), ]
        out[names(edges)] <- imp
      }
    } else {
      fit <- randomForest(x, y, importance = TRUE, ...)
      imp <- importance(fit, type = 1L, scale = TRUE)
      out[rownames(imp)] <- as.numeric(imp)
    }
    return(out) 
  }

  # Execute in parallel
  adj_mat <- foreach(j = seq_len(p), .combine = cbind) %dopar% edge_wts(j)
  dimnames(adj_mat) <- list(colnames(dat)[seq_len(p)], colnames(dat))
  if (is.null(focus) && !directed) {
    adj_mat <- pmax(adj_mat, t(adj_mat))
  }
  
  # Export
  return(adj_mat)
  
}
