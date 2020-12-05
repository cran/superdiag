########################### SUPERDIAG FUNCTION #################################
### UPDATE: 11/28/2020; Le Bao and Jeff Gill
### UPDATE: 10/31/2011; Tsung-han Tsai and Jeff Gill

#' Test for Markov Chain Nonconvergence
#'
#' The \code{superdiag} function takes MCMC samples as input. It provides a comprehensive
#'   test suite for Markov chain nonconvergence, which integrates five standard
#'   empirical MCMC convergence diagnostics: Gelman-Rubin, Geweke, Heidelberger-Welch,
#'   Raftery-Lewis, Hellinger distance. It can also present trace plots and density
#'   histograms along with the diagnostics as options.
#'
#' @param mcmcoutput A \code{mcmc} or \code{mcmc.list} object.
#' @param burnin The number of burn-in iterations. Defaults to half of the chain(s).
#' @param terms The convergence diagnostic methods. Defaults to all five methods.
#'   Users can also specify one or more particular methods chosen from "\code{geweke}",
#'   "\code{heidel}", "\code{raftery}", "\code{gelman}", "\code{hellinger}".
#' @param plot Logical values indicates whether a graph of \code{superdiagPlot}
#'   should be presented with a run of the function. Defaults to \code{FALSE}.
#' @param confidence.gr (1-\eqn{\alpha})\% for the Gelman and Rubin test. The upper
#'  95\% credible interval is the default.
#' @param frac1.gw frac1 for the Geweke test. The proportion of the early era of
#'   the chain, defaulted to 0.1.
#' @param frac2.gw frac2 for the Geweke test. The proportion of the late era of
#'   the chain, defaulted to 0.5.
#' @param eps.hw \code{epsilon} for the Heidelberger and Welch test. The accuracy
#'   parameter determines whether the halfwidth is passed or not, defaulted to 0.1.
#' @param pvalue.hw p-value for the Heidelberger and Welch test. The halfwidth of
#'   the test calculates a (1-\eqn{\alpha})\% credible interval around the sample
#'   mean for each parameter dimension. The default is 0.05.
#' @param q.rl q-parameter for the Raftery and Lewis test. The posterior tail threshold
#'   of interest, defaulted to 0.025.
#' @param r.rl r-parameter for the Raftery and Lewis test. The tolerance for the
#'   tail threshold, defaulted to 0.0005.
#' @param s.rl s-parameter for the Raftery and Lewis test. The desired probability
#'   of being within the tolerance, defaulted to 0.95.
#' @param eps.rl convergence epsilon for the Raftery and Lewis test. The convergence
#'   tolerance value, which is used to determine a stopping point based on a parallel chain process, defaulted to 0.001.
#' @param bins Number of bins for within-chain Hellinger distance test. Defaults to 5.
#' @param binwidth Alternative specification for \code{bins}. The size of bin for each batch
#'   of the chain(s) for computing the distance between batches.
#'
#' @return A list object including the results for all the diagnostics. A \code{superdiagPlot}
#'   including both traceplot(s) and density histogram will also be returned if \code{plot=TRUE}.
#'
#' @details If only one chain is analyzed, the default settings in \code{boa} and
#'   \code{coda} for Geweke, Raftery-Lewis, and Heidelberger-Welch diagnostics are
#'   used. If multiple chains are provided, only the first chain uses the defaults
#'   and all other chain analyses get random values as follows. For Geweke test,
#'   random non-overlapping proportions up from the start of the chain and down
#'   from the end of the chain are generated. For Heidelberger-Welch test, the
#'   value of \code{pvalue.hw} is sampled with replacement from common \eqn{\alpha}
#'   values; the value of \code{eps.hw} is sampled uniformly in the interval [0.01:0.2].
#'   For Raftery-Lewis test, each of these four parameters are sample from a vector
#'   (changeable by users) of values around the defaults (larger and smaller) to
#'   provide a reasonable range of alternatives. For Hellinger distance, if only
#'   one chain is analyzed, only within-chain distance will be reported.
#'
#' @references
#'   Tsai, Tsung-han and Gill, Jeff (2012). \dQuote{superdiag: A Comprehensive Test
#'     Suite for Markov Chain Non-Convergence.} \emph{The Political Methodologist},
#'     19 (2), 12-18.
#'
#'   Plummer, Martyn, Nicky Best, Kate Cowles, and Karen Vines (2006). "CODA: convergence
#'      diagnosis and output analysis for MCMC." \emph{R news}, 6 (1), 7-11.
#'
#' @seealso \code{\link{superdiagPlot}}, \code{\link[coda]{gelman.diag}}, \code{\link[coda]{heidel.diag}},
#'   \code{\link[coda]{raftery.diag}}, \code{\link[coda]{gelman.diag}}
#'
#' @examples
#' \dontrun{
#'
#' data(tobit.list)
#' summary(tobit.list[1])
#'
#' # FOR mcmc.list OBJECT
#' superdiag(tobit.list, burnin=0)
#'
#' # FOR mcmc OBJECT
#' tobit.diag <- superdiag(tobit.list[[1]], burnin=0, plot=TRUE)
#' tobit.diag
#' }
#'
#' @importFrom coda as.mcmc as.mcmc.list as.array.mcmc.list geweke.diag heidel.diag raftery.diag gelman.diag
#' @importFrom stats runif density
#' @importFrom utils combn tail
#'
#' @export
superdiag <- function(mcmcoutput, burnin, terms = "all", plot = FALSE,
                      confidence.gr=0.95, frac1.gw=0.1, frac2.gw=0.5, eps.hw=0.1,
                      pvalue.hw=0.05, q.rl=0.025, r.rl=0.005, s.rl=0.95, eps.rl=0.001,
                      bins=5, binwidth=NULL) {

	# mcmcoutput: input chains from jags, bugs, etc.
	# burnin: the number of burn-in iterations for the sampler
	# confidence.gr: 1-alpha for the Gelman and Rubin test
	# frac1.gw: frac1 for the Geweke test
	# frac2.gw: frac2 for the Geweke test
	# eps.hw: epsilon for the Heidelberger and Welch test
	# pvalue.hw: p-value for the Heidelberger and Welch test
	# q.rl: q-parameter for the Raftery and Lewis test
	# r.rl: r-parameter for the Raftery and Lewis test
	# s.rl: s-parameter for the Raftery and Lewis test
	# eps.rl: convergence epsilon for the Raftery and Lewis test
  # bins: bins for within-chain Hellinger distance
  # binwidth: alternative for bins, defaults to NULL

	# THE INPUTS SHOULD BE SAVED AS "mcmc" OR "mcmc.list" CLASS
  if (!inherits(mcmcoutput, "mcmc") & !inherits(mcmcoutput, "mcmc.list") & !inherits(mcmcoutput, "list"))
    stop("The inputs have to be mcmc, mcmc.list, or list objects.")

	# CONVERT "mcmc" INTO "mcmc.list"
	if (inherits(mcmcoutput, "mcmc")) {
		mcmcoutput <- as.mcmc.list(mcmcoutput)
	}

  if (terms == "all") terms <- c("geweke", "heidel", "raftery", "gelman", "hellinger")

	para.names <- dimnames(mcmcoutput[[1]])[[2]]  # PARAMETERS NAMES
	n.chains <- length(mcmcoutput)  # THE NUMBER OF CHAINS
	dim.chain <- sapply(mcmcoutput, dim)  # THE DIMENSION OF EACH CHAIN
	t.iter <- dim.chain[1]	  # THE NUMBER OF ITERATIONS BEFORE DELETING THE BURN-IN PERIOD
	diff.dim <- dim.chain - dim.chain
	if (sum(diff.dim != 0) != 0) stop("The number of iterations or variables is not equal for all chains.")

	# DISCARD THE BURN-IN PERIOD
	if (missing(burnin)) burnin <- round(t.iter/2)

	if (burnin <= 0) {
	  if (burnin < 0) warning("The burn-in period is negative. 'burnin = 0' is used.")
	  mcmcburnin <- mcmcoutput
	} else {
	  if (burnin >= t.iter) {
	    stop("The number of iterations is less than the burn-in period.")
	  } else {
	    mcmcburnin <- lapply(mcmcoutput, burn, burnin=burnin)
	  }
	}

	# SAVE THE SAMPLES AS A MCMC LIST AFTER DISCARDING THE BURN-IN PERIOD
	mcmcburnin.list <- vector("list", n.chains)
	for (i in 1:n.chains) {
	  mcmcburnin.list[[i]] <- as.mcmc(mcmcburnin[[i]])
	}
	mcmcburnin.mcmclist <- as.mcmc.list(mcmcburnin.list)

	# THE TOTAL SAMPLES FOR ALL CHAINS
	t.samples <- dim(as.matrix(mcmcburnin.mcmclist))[1]

	# OUTPUT LIST
	out <- list(n.chains = n.chains, n.iter = t.iter, burnin = burnin,
	            total.sample = t.samples, terms = terms,
	            mcmc.samples = mcmcburnin.mcmclist)

	# ROW NAMES FOR CHAINS
	chain.name <- paste("chain ", 1:n.chains, sep="")

	# GEWEKE
	if ("geweke" %in% terms) {
	geweke.chains <- matrix(NA, nrow=n.chains, ncol=dim.chain[2])
	### SETUP DIFFERENT WINDOW SPECIFICATIONS FOR GEWEKE
	geweke.windows <- matrix(c(frac1.gw,frac2.gw),ncol=2)
	if (n.chains > 1){
	for (i in 2:n.chains) {
	  win1 <- runif(1,0,0.99); win2 <- 1-runif(1,win1,1)
	  geweke.windows <- rbind(geweke.windows,c(win1,win2))
	}
	}
	### RUN DIAGNOSTICS BY CHAIN IN ORDER TO AVOID ONE CHAIN RUINS ALL CHAINS
	for (i in 1:n.chains) {
	  geweke <- suppressWarnings(try(geweke.diag(mcmcburnin.mcmclist[[i]], geweke.windows[i,1], geweke.windows[i,2]), silent=TRUE));
	  if (class(geweke) == "geweke.diag")  geweke.chains[i,] <- t(geweke[1]$z);
	}
	colnames(geweke.chains) <- para.names
	colnames(geweke.windows) <- c("From start", "From stop")
	rownames(geweke.windows) <- chain.name
	rownames(geweke.chains) <- chain.name
	geweke.list <- list(windows = geweke.windows, chains = geweke.chains)
	out[["geweke"]] <- geweke.list
	}

	# HEIDELBERGER AND WELCH
	if ("heidel" %in% terms) {
	heidel.list <- vector("list", n.chains)
	### SETUP DIFFERENT PARAMETER SPECIFICATIONS FOR HEIDELBERGER AND WELCH
	heidel.params <- matrix(c(eps.hw,pvalue.hw),ncol=2)
	if (n.chains > 1){
	pvals <- c(0.1,0.05,0.025,0.01,0.005)
	for (i in 2:n.chains) {
	  param1 <- runif(1,0.01,0.2); param2 <- sample(x=pvals,size=1)
	  heidel.params <- rbind(heidel.params,c(param1,param2))
	}
	}
	rownames(heidel.params) <- 	chain.name
	colnames(heidel.params) <-  c("epsilon", "alpha")
	### RUN DIAGNOSTICS BY CHAIN IN ORDER TO AVOID ONE CHAIN RUINS ALL CHAINS
	for (i in 1:n.chains) {
	  heidel.list[[i]] <- suppressWarnings(try(heidel.diag(mcmcburnin.mcmclist[[i]],
	                                                       heidel.params[i,1],
	                                                       heidel.params[i,2]),
	                                           silent=TRUE))
	}
	out[["heidel"]] <- list(heidel.params = heidel.params, heidel.list = heidel.list)
	}

	# RAFTERY AND LEWIS
	if ("raftery" %in% terms) {
	raftery.list <- vector("list", n.chains)
	# SETUP DIFFERENT PARAMETER SPECIFICATIONS FOR RAFTERY AND LEWIS
	raft.params <- matrix(c(q.rl, r.rl, s.rl, eps.rl),ncol=4)
	if (n.chains > 1){
	qvals <-c(0.25,0.1,0.05,0.01,0.001)
	rvals <- c(0.001,0.0025,0.0005,0.001,0.005)
	svals <- c(0.9,0.95,0.975,0.99,0.999)
	evals <- c(0.005,0.0025,0.001,0.0005,0.0002)

	for (i in 2:n.chains) {
	  param1 <- sample(x=qvals,size=1); param2 <- sample(x=rvals,size=1)
	  param3 <- sample(x=svals,size=1); param4 <- sample(x=evals,size=1)
	  raft.params <- rbind(raft.params,c(param1,param2,param3,param4))
	}
	}
	rownames(raft.params) <- 	chain.name
	colnames(raft.params) <-  c("Quantile (q)", "Accuracy (r)", "Probability (s)", "Convergence eps")
	### RUN DIAGNOSTICS BY CHAIN IN ORDER TO AVOID ONE CHAIN RUINS ALL CHAINS
	for (i in 1:n.chains) {
	  raftery.list[[i]] <- suppressWarnings(try(raftery.diag(mcmcburnin.mcmclist[[i]], q=raft.params[i,1],r=raft.params[i,2],
	                                                         s=raft.params[i,3],converge.eps=raft.params[i,4]), silent=TRUE));
	}
	out[["raftery"]] <- list(raftery.params = raft.params, raftery.list = raftery.list)
	}

	if ("gelman" %in% terms & n.chains > 1){
	  gelman.list <- gelman.diag(mcmcburnin.mcmclist, confidence.gr)
	  out[["gelman"]] <- gelman.list
	}

	if ("hellinger" %in% terms){
	  if (!missing(bins) & !missing(binwidth)) warning("Both 'bins' and 'binwidth' are specified. 'binwidth' will be suppressed.")
	  if (missing(bins) & !missing(binwidth)) bins <- floor((t.iter-burnin)/binwidth)

	  if (n.chains < 2) {
	    hwithin.one <- suppressWarnings(try(hellinger.diag.matrix(as.matrix(mcmcburnin.mcmclist), bins=bins)))
	    hellinger.list <- list(within = list(chain1 = hwithin.one))
	  } else {
	    hbetween <- suppressWarnings(try(hellinger.diag.list(mcmcburnin.mcmclist)))
	    hwithin <- list()
	    for (i in 1:n.chains){
	      hwithin[[i]] <- suppressWarnings(try(hellinger.diag.matrix(mcmcburnin.mcmclist[[i]], bins=bins)))
	    }
	    names(hwithin) <- paste0("chain", 1:n.chains)
	    hellinger.list <- list(between = hbetween, within = hwithin)
	  }
	  out[["hellinger"]] <- hellinger.list
	}

	# PLOT
	if (plot == TRUE) do.call("superdiagPlot", list(mcmcoutput = out$mcmc.samples, burnin = 0))

	# OUTPUT
	class(out) <- c("superdiag", class(out))
	out
}

#' @method print superdiag
#' @export
print.superdiag <- function(x, digits = min(4, .Options$digits), ...) {
  cat(paste("Number of chains =", x$n.chains, "\n"))
  cat(paste("Number of iterations =", x$n.iter, "per chain before discarding the burn-in period\n"))
  cat(paste("Burn-in period =", x$burnin, "per chain\n"))
  cat(paste("Sample size in total =", x$total.sample, "\n"))
  cat("\n");

  if (x$n.chains < 2) {
    cat("****************** The Geweke diagnostic: ******************\n")
    cat(paste("Fraction from start =", x$geweke$windows[1,1], "\n"))
    cat(paste("Fraction from stop =", x$geweke$windows[1,2], "\n"))
    cat("\n")
    cat(paste("Z-scores:\n"))
    print(t(x$geweke$chains), digits = digits, ...)
    cat("\n")

    message("The Gelman-Rubin diagnostic is not reported since the number of chains is less than 2.\n")

    cat("************* The Heidelberger-Welch diagnostic ************\n")
    cat(paste0("epsilon = ", x$heidel$heidel.params[1,1], ", alpha = ", x$heidel$heidel.params[1,2], "\n"))
    print(x$heidel$heidel.list[[1]], digits = digits, ...)
    cat("\n")

    cat("*************** The Raftery-Lewis diagnostic ***************")
    print(x$raftery$raftery.list[[1]], digits = digits, ...)
    cat("\n")

    cat("************* The Hellinger distance diagnostic ************\n")
    cat(paste("Within chain: \n"))
    print(x$hellinger[[1]][[1]], digits = digits, ...)
    cat("\n")
    message("The between-chain Hellinger distance is not reported since the number of chains is less than 2.\n")

  } else {
    cat("****************** The Geweke diagnostic: ******************\n")
    cat(paste("Windows:\n"))
    print(round(t(x$geweke$windows), digits = digits), ...)
    cat(paste("\nZ-scores:\n"))
    print(t(x$geweke$chains), digits = digits, ...)
    cat("\n")

    cat("*************** The Gelman-Rubin diagnostic: ***************\n")
    cat(paste0("Potential scale reduction factors:\n"))
    print(x$gelman$psrf, digits = digits, ...)
    cat(paste0("\nMultivariate psrf: ", round(x$gelman[[2]], digits = digits), "\n"))
    cat("\n")

    cat("************* The Heidelberger-Welch diagnostic ************\n");
    #cat("\n")
    for (i in 1:x$n.chains)  {
      cat(paste("Chain ",i,":\nepsilon=",round(x$heidel$heidel.params[i,1],3),", alpha=",round(x$heidel$heidel.params[i,2],3),sep=""))
      print(x$heidel$heidel.list[[i]], digits = digits, ...)
      cat("\n")
    }

    cat("*************** The Raftery-Lewis diagnostic ***************\n");
    #cat("\n");
    for (i in 1:x$n.chains)  {
      cat(paste0("Chain ",i, ":\n"))
      cat(paste0("Convergence eps = ",x$raftery$raftery.params[i,4]))
      print(x$raftery$raftery.list[[i]], digits = digits, ...)
      ### Get rid of an empty line from print.raftery.diag()
      if (x$raftery$raftery.list[[i]]$resmatrix[1] == "Error") cat("\n")
    }

    cat("************* The Hellinger distance diagnostic ************\n")
    cat(paste("Between chains: \n"))
    print.hellinger.diag(x$hellinger$between)
    cat("\n")

    for (i in 1:x$n.chains)  {
    cat(paste0("Within chain ", i, ":\n"))
    print.hellinger.diag(x$hellinger$within[[i]], digits = digits, ...)
    if (i < x$n.chains) cat("\n")
    }
    cat("\n")
  }
}
