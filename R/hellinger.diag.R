##################### Hellinger Distance Diagnostics ###########################
### UPDATE: 11/28/2020; Le Bao
### The code is partiall adapted from bmk package

#'  Hellinger Distance Diagnostics for MCMC
#'
#'  This computes the Hellinger distance for \code{mcmc} or \code{mcmc.list} objects.
#'
#'  @param x A \code{mcmc} or \code{mcmc.list} object. For mcmc object contains
#'    one single MCMC chain, the within-chain distance is computed. For \code{mcmc.list}
#'    object containing multiple MCMC chains, the between distribution distance(s)
#'    is computed.
#'  @param bins Number of bins for testing within-chain distance. Defaults to 10.
#'  @param binwidth Alternative specification for \code{bins}. The size of bin for each batch
#'   of the chain(s) for computing the distance between batches.
#'
#' @return A matrix of the estimates.
#'
#' @references Boone, Edward L., Jason RW Merrick, and Matthew J. Krachey (2014). "A
#'   Hellinger distance approach to MCMC diagnostics." \emph{Journal of Statistical Computation
#'   and Simulation}, 84 (4), 833-849, \code{doi:10.1080/00949655.2012.729588}.
#'
#' @examples
#' \dontrun{
#'
#' data(tobit.list)
#' summary(tobit.list[1])
#'
#' # FOR mcmc.list OBJECT
#' hellinger.diag.mcmc.list(tobit.list)
#'
#' # FOR mcmc OBJECT
#' hellinger.diag.mcmc(tobit.list[[1]])
#'
#' }
#'
#' @keywords internal
#'
hellinger.diag <- function(x) UseMethod("hellinger.diag")

#' @rdname hellinger.diag
#' @method hellinger.diag list
hellinger.diag.list <- function(x) {
  xa <- as.array.mcmc.list(x)
  n.iter <- dim(xa)[1]
  n.var <- dim(xa)[2]
  n.chain <- dim(xa)[3]
  commat <- combn(1:n.chain,2)
  mat <- matrix(NA, nrow = n.var, ncol = 2)
  for (i in 1:n.var){
    out <- NULL
    for (j in 1:ncol(commat)){
      if (j == 1) {out <- hdist(xa[,i,commat[1,j]], xa[,i,commat[2,j]])
      } else {
        out <- c(out, hdist(xa[,i,commat[1,j]], xa[,i,commat[2,j]]))
      }
    }
    mat[i,] <- c(min(out), max(out))
  }
  rownames(mat) <- colnames(xa[,,1])
  colnames(mat) <- c("Min", "Max")
  class(mat) <- c("hellinger.diag")
  mat
}

#' @rdname hellinger.diag
#' @method hellinger.diag mcmc.list
hellinger.diag.mcmc.list <- hellinger.diag.list

#' @rdname hellinger.diag
#' @method hellinger.diag matrix
hellinger.diag.matrix <- function(x, bins, binwidth) {
  nvar <- ncol(x)
  niter <- nrow(x)
  if (missing(bins) & missing(binwidth)) bins <- 10
  if (missing(bins) & !missing(binwidth)) bins <- round(niter/binwidth)
  size <- niter%/%bins
  if (niter%%bins != 0) {
    ndel <- niter - size * bins; x <- tail(x, size * bins); niter <- nrow(x)
    warning("The number of iterations is not a multiple of bins. The first ", paste(ndel), " iterations will be discarded.")
  }
  interval <- seq(size, niter, by = size)
  mat <- matrix(NA, nrow = nvar, ncol = length(interval)-1)
  mat
  for(i in 1:nvar){
    for (j in 1:(bins-1)){
      mat[i, j] <- hdist(x[1:interval[j],i], x[(interval[j]+1):interval[j+1], i])
    }
  }
  rownames(mat) <- colnames(x)
  colnames(mat) <- interval[1:ncol(mat)]
  class(mat) <- c("hellinger.diag")
  mat
}

#' @rdname hellinger.diag
#' @method hellinger.diag mcmc
hellinger.diag.mcmc <- hellinger.diag.matrix


#' @method print hellinger.diag
print.hellinger.diag <- function(x, digits = min(4, .Options$digits), ...) {
  class(x) <- "double"
  print.default(x, digits = digits, ...)
  invisible(x)
}

# HELPER: Hellinger distance
hdist <- function(v1, v2, bins=NULL){
  if (is.null(bins)) bins <- length(v1)
  vec <- c(v1, v2); min <- min(vec); max <- max(vec)
  P <- density(v1, n=bins, from=min, to=max)
  Q <- density(v2, n=bins, from=min, to=max)
  interval <- P$x[2] - P$x[1]
  hdist <- sqrt(sum((sqrt(P$y)-sqrt(Q$y))^2*interval)/2)
  hdist
}
