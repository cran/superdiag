######################## SUPERDIAG PLOT FUNCTION ###############################
### CREATED: 11/28/2020; Le Bao

#' Summary Plots for \code{mcmc} Objects
#'
#' The \code{superdiagPlot} function takes \code{mcmc} or \code{mcmc.list} as input.
#'   It provides summary plots including trace plot(s) and density histogram(s) for
#'   each variable in the MCMC chain(s).
#'
#' @param mcmcoutput A \code{mcmc} or \code{mcmc.list} object.
#' @param burnin The number of burn-in iterations. Defaults to half of the chain(s).
#' @param params The parameters to be summarized. Defaults to all. Users can specify
#'   a subset of variables by name or index to reduce the size of plots.
#' @param title A character vector specifies the title of the plot.
#' @param col A vector or list specifies the color schemes of the plot. For \code{mcmc}
#'   object, the default is a two-color grey scheme; Users can specify any numbers of color
#'   no more than the number of variables. For \code{mcmc.list} object, it requires a list
#'   or vector with multiple colors to examine the mixing. The default is a list of two sets
#'   of colors.
#' @param trace.options Additional options for trace plots. It can take arguments as
#'   in \code{plot} and \code{lines} functions.
#' @param density.options Additional options for density histograms. It can take arguments
#'   as in \code{hist} (e.g \code{breaks}, \code{freq}, \code{probability}) and \code{barplot}
#'   functions.
#' @param title.options Additional options for title. It can take arguments as in
#'   \code{mtext} function.
#'
#' @return A summary plot for each variable in the MCMC chain(s).
#'
#'
#' @seealso \code{\link{superdiag}}
#'
#' @examples
#' \dontrun{
#'
#' data(tobit.list)
#' summary(tobit.list[1])
#'
#' # FOR mcmc.list OBJECT
#' superdiagPlot(tobit.list, burnin=0)
#'
#' # FOR mcmc OBJECT
#' superdiagPlot(tobit.list[[1]], burnin=0, col=c("grey25", "dodgerblue"),
#'   title = "Tobit Model (Chain 1)", title.options=list(cex=1.2))
#' }
#'
#' @importFrom coda as.mcmc as.mcmc.list
#' @importFrom graphics plot lines hist barplot axis mtext layout matplot par box
#' @export
superdiagPlot <- function(mcmcoutput, burnin, params, title, col, trace.options,
                          density.options, title.options) {
  UseMethod("superdiagPlot")
}

#' @rdname superdiagPlot
#' @export
superdiagPlot.mcmc <- function(mcmcoutput, burnin, params = "all", title, col=c("grey25", "grey75"),
                               trace.options=list(trace.axis = TRUE, box = TRUE), density.options=list(),
                                title.options=list()) {
  # CHECK IMPUT
  ### COLLAPSE THE MCMC LIST WITH ONLY ONE CHAIN TO A MCMC
  if (inherits(mcmcoutput, "mcmc.list") & (length(mcmcoutput) == 1))  mcmcoutput <- mcmcoutput[[1]]
  if (!inherits(mcmcoutput, "mcmc")) stop("The input is not a 'mcmc' object")

  # DISCARD THE BURN-IN PERIOD
  t.iter <- nrow(mcmcoutput)
  if (missing(burnin)) burnin <- round(t.iter/2)
  if (burnin <= 0) {
    if (burnin < 0) warning("The burn-in period is negative. 'burnin = 0' is used.")
    mcmcoutput <- mcmcoutput
  } else {
    if (burnin >= t.iter) {
      stop("The number of iterations is less than the burn-in period.")
    } else {
      mcmcoutput <- burn(mcmcoutput, burnin=burnin)
    }
  }

  # ON EXIT SETTING
  suppressWarnings(on.exit(layout(1)))
  opar<-par(no.readonly = TRUE) #"mar","mgp", "oma",
  on.exit(par(opar),add = TRUE)

  # SELECTED PARAMETERS
  if (params != "all") {
    if (is.character(params))
      params <- c(1:dim(mcmcoutput[[1]])[2])[unlist(dimnames(mcmcoutput[[1]])[2]) %in% params]
    for (i in 1:length(mcmcoutput)) {
      mcmcoutput[[i]] <- mcmcoutput[[i]][,params]
    }
  }

  # STANDING PARAMETERS AFTER BURN IN AND SELECTION
  param.names <- colnames(mcmcoutput)
  n.var <- dim(mcmcoutput)[2]
  n.iter <- dim(mcmcoutput)[1]

  # OTHER SPECIFICATIONS FOR PLOT
  ### COLOR
  if (length(col) != 2 & length(col) != n.var){
    col <- c("grey25", "grey75")
    warning("The length of 'col' is not equal to the number of parameters. Use default values instead.")}
  if (length(col) != n.var) col <- rep(col, n.var)[1:n.var]

  ### AXES
  trace.axis <- ifelse("trace.options" %in% names(trace.options), trace.options$trace.axis, TRUE)
  trace.box <- ifelse("box" %in% names(trace.options), trace.options$box, TRUE)
  if (is.numeric(trace.axis)) {breaks <- trace.axis
  } else {size <- n.iter%/%5; breaks <- seq(0, n.iter, size)}

  ### OPTIONS FOR LINES
  trace.options <- trace.options[! names(trace.options) %in% c("trace.axis", "box")]
  lines.def <- list(xlab = '', xaxt='n', axes = FALSE, cex.axis = .8) #yaxt='n'
  trace.options <- c(trace.options, lines.def[!names(lines.def) %in% names(trace.options)])

  ### OPTIONS FOR HISTOGRAM
  hdata.options <- density.options[c("breaks", "freq", "probability") %in% names(density.options)]
  bar.user.options <- density.options[!names(density.options) %in% c("breaks", "freq", "probability")]
  hist.def <- list(axes = TRUE, xaxt='n', space = 0, horiz=TRUE)
  bar.options <- c(bar.user.options, hist.def[!names(hist.def) %in% names(bar.user.options)])

  ### TITLE
  title.def <- list(side = 3, line = 0, font = 2, outer = TRUE, padj=-1)
  title.options <- c(title.options, title.def[!names(title.def) %in% names(title.options)])

  # FOR CHAINS WITH MANY VARIABLES
  max.var <- 10 ### MAX NUMBER OF VARIABLES WITHOUT SEPARATING INTO MULTIPLE PAGES
  one.page <- 6 ### NUMBER OF VARS ON ONE PAGE

  ### CALCULATING THE NUMBER OF VARS
  if (n.var > max.var) {
    message("The number of parameters is too large. The plots will be shown on seperate pages.")
    multi <- ifelse(n.var%%one.page==0, n.var%/%one.page, n.var%/%one.page + 1)
    lay.ls <- list()
    for (i in 1:multi) {
    if(i < multi) {lay.ls[[i]] <- seq.int((i-1)*one.page+1,one.page*i,1)} else {lay.ls[[i]] <- seq.int((i-1)*one.page+1,n.var,1)}
    }
  } else {
    lay.ls <- list(1:n.var)
  }

  ### SUBSET
  for (j in 1:length(lay.ls)){
  new.nvar <- length(lay.ls[[j]])
  new.mat <- as.matrix(mcmcoutput[,lay.ls[[j]]])
  new.params <- param.names[lay.ls[[j]]]

  # PLOT
  ### LAYOUT
  layout.matrix <- cbind(matrix(rep(seq(1,2*new.nvar,2), each = 5), nrow = new.nvar, byrow = TRUE),
                         matrix(rep(seq(2,2*new.nvar,2), each = 1), nrow = new.nvar, byrow = TRUE))

  layout(layout.matrix)
  if (!missing(title)) {
    suppressWarnings(par(mgp=c(-1,1,0), oma=c(2,1,3,1) + 0.1))
  } else {
    suppressWarnings(par(mgp=c(-1,1,0), oma=c(2,1,1,1) + 0.1))
  }

  for (i in 1:new.nvar){
    var <- as.vector(new.mat[,i])
    if (trace.box == TRUE) {par(mar=c(0,3,0,0));ypara <- ""} else {
      par(mar=c(0,1,0,0));ypara<-new.params[i]}
    ### TRACE PLOT
    line.list <- c(list(x = var, type='l', col = col[i], ylab=ypara), trace.options)
    do.call(plot, line.list)
    if (trace.box == TRUE) mtext(new.params[i], side=2, line = 0, cex=trace.options$cex.axis-.1,
                           padj = -(trace.options$cex.axis)/0.225)

    ### AEXES
    if (trace.box == TRUE) axis(2, cex.axis=trace.options$cex.axis)
    if (i == new.nvar & (trace.axis[1] == TRUE | is.numeric(trace.axis))) {
      axis(1, at = breaks, cex.axis= trace.options$cex.axis+0.1, las=1, outer = TRUE)
    }
    if (trace.box == TRUE) box()
    ### DENSITY PLOT
    par(mar=c(0,.5,0,0))
    hdata.list <- c(list(x=var, plot=FALSE), hdata.options)
    hdata <- do.call(hist, hdata.list)
    hist.list <- c(list(height=hdata$counts, col = col[i]), bar.options)
    do.call(barplot, hist.list)
  }
  if (!missing(title)) {
    title.list <- c(list(text = title), title.options)
    do.call(mtext, title.list)
  }
  }
}

#' @rdname superdiagPlot
#' @export
superdiagPlot.mcmc.list <- function(mcmcoutput, burnin, params = "all", title, col,
                                     trace.options=list(trace.axis = TRUE, box = TRUE),
                                    density.options=list(), title.options=list()) {
  # CHECK IMPUT
  if (!inherits(mcmcoutput, "mcmc.list")) stop("The input is not a 'mcmc.list' object")
  ### IF LIST ONLY CONTAINS ONE CHAIN, PASS TO superdiagPlot.mcmc()
  if (length(mcmcoutput) == 1) {
    .args <- as.list(match.call()[-1])
    do.call(superdiagPlot.mcmc, .args)
  } else {

  # ON EXIT SETTING
    suppressWarnings(on.exit(layout(1)))
    opar<-par(no.readonly = TRUE) #"mar","mgp", "oma",
    on.exit(par(opar),add = TRUE)

  # ACTUAL FUNCTION FOR MCMC LIST STARTS HERE

  # PARAMETERS
  n.chains <- length(mcmcoutput)  # THE NUMBER OF CHAINS
  dim.chain <- sapply(mcmcoutput, dim)  # THE DIMENSION OF EACH CHAIN
  diff.dim <- dim.chain - dim.chain
  if (sum(diff.dim != 0) != 0) stop("The number of iterations or variables is not equal for all chains.")
  t.iter <- dim.chain[1]	  # THE NUMBER OF ITERATIONS BEFORE DELETING THE BURN-IN PERIOD

  # DISCARD THE BURN-IN PERIOD
  if (missing(burnin)) burnin <- round(nrow(mcmcoutput[[1]])/2)
  if (burnin <= 0) {
    if (burnin < 0) warning("The burn-in period is negative. 'burnin = 0' is used.")
    mcmcoutput <- mcmcoutput
  } else {
    if (burnin >= t.iter) {
      stop("The number of iterations is less than the burn-in period.")
    } else {
      mcmcoutput <- lapply(mcmcoutput, burn, burnin=burnin)
    }
  }

  # SELECTED PARAMETERS
  if (params != "all") {
    if (is.character(params))
      params <- c(1:dim(mcmcoutput[[1]])[2])[unlist(dimnames(mcmcoutput[[1]])[2]) %in% params]
    for (i in 1:length(mcmcoutput)) {
      mcmcoutput[[i]] <- mcmcoutput[[i]][,params]
    }
  }

  # STANDING PARAMETERS AFTER BURN IN AND SELECTION
  param.names <- dimnames(mcmcoutput[[1]])[[2]]
  n.var <- dim(mcmcoutput[[1]])[2]
  n.iter <- dim(mcmcoutput[[1]])[1]

  # COVERT TO NEW LISTS
  trace.mat.list <- list()
  for (j in 1:length(param.names)){
    tmp.mat <- matrix(NA, nrow=n.iter, ncol=n.chains)
    for (i in 1:n.chains){
      tmp.mat[,i] <- mcmcoutput[[i]][,j]
    }
    trace.mat.list[[j]] <- tmp.mat
  }
  names(trace.mat.list) <- param.names

  # OTHER SPECIFICATIONS FOR PLOT
  ### AXES
  trace.axis <- ifelse("trace.options" %in% names(trace.options), trace.options$trace.axis, TRUE)
  trace.box <- ifelse("box" %in% names(trace.options), trace.options$box, TRUE)
  if (is.numeric(trace.axis)) {breaks <- trace.axis
  } else {size <- n.iter%/%5; breaks <- seq(0, n.iter, size)}

  ### OPTIONS FOR LINES
  trace.options <- trace.options[! names(trace.options) %in% c("trace.axis", "box")]
  lines.def <- list(yaxt='n', axes = FALSE, xlab = '', cex.axis = .8)
  trace.options <- c(trace.options, lines.def[!names(lines.def) %in% names(trace.options)])

  ### OPTIONS FOR HISTOGRAM
  hdata.options <- density.options[c("breaks", "freq", "probability") %in% names(density.options)]
  bar.user.options <- density.options[!names(density.options) %in% c("breaks", "freq", "probability")]
  hist.def <- list(axes = TRUE, xaxt='n', space = 0, horiz=TRUE)
  bar.options <- c(bar.user.options, hist.def[!names(hist.def) %in% names(bar.user.options)])

  ### TITLE
  title.def <- list(side = 3, line = 0, font = 2, outer = TRUE, padj=-1)
  title.options <- c(title.options, title.def[!names(title.def) %in% names(title.options)])

  ### COLOR
  if (missing(col)) {col.list <- list(c(1:n.chains), c((n.chains+1):(2*n.chains)))
  } else if (!missing(col) & !is.list(col)) {
      if (length(col) < n.chains) col <- rep(col, n.chains*2)
      col.list[[1]] <- col[1:n.chains]; col.list[[2]] <- tail(col, length(col)-n.chains)
    } else {
      col.list <- col
    }
  if (length(col.list) < n.var) col.list <- rep(col.list, n.var)[1:n.var]

  max.var <- 10 ### MAX NUMBER OF VARIABLES WITHOUT SEPARATING INTO MULTIPLE PAGES
  one.page <- 6 ### NUMBER OF VARS ON ONE PAGE

  ### CALCULATING THE NUMBER OF VARS
  if (n.var > max.var) {
    message("The number of parameters is too large. The plots will be shown on seperate pages.")
    multi <- ifelse(n.var%%one.page==0, n.var%/%one.page, n.var%/%one.page + 1)
    lay.ls <- list()
    for (i in 1:multi) {
      if(i < multi) {lay.ls[[i]] <- seq.int((i-1)*one.page+1,one.page*i,1)} else {lay.ls[[i]] <- seq.int((i-1)*one.page+1,n.var,1)}
    }
  } else {
    lay.ls <- list(1:n.var)
  }

  # SUBSET
  for (j in 1:length(lay.ls)){
    new.nvar <- length(lay.ls[[j]])
    new.mcmc.list <- trace.mat.list[lay.ls[[j]]]
    new.params <- param.names[lay.ls[[j]]]

  # PLOT
    ### LAYOUT
    layout.matrix <- cbind(matrix(rep(seq(1,2*new.nvar,2), each = 5), nrow = new.nvar, byrow = TRUE),
                           matrix(rep(seq(2,2*new.nvar,2), each = 1), nrow = new.nvar, byrow = TRUE))

    layout(layout.matrix)
    if (!missing(title)) {
      suppressWarnings(par(mgp=c(-1,1,0), oma=c(2,1,3,1) + 0.1))
    } else {
      suppressWarnings(par(mgp=c(-1,1,0), oma=c(2,1,1,1) + 0.1))
    }

    for (i in 1:new.nvar){
      vmat <- new.mcmc.list[[i]]
      if (trace.box == TRUE) {par(mar=c(0,3,0,0));ypara <- ""} else {par(mar=c(0,.25,0,0));ypara<-new.params[i]}
      ### TRACE PLOT
      mat.list <- c(list(x = vmat, type='l', col = col.list[[i]], ylab=ypara), trace.options)
      do.call(matplot, mat.list)
      if (trace.box == TRUE) mtext(new.params[i], side=2, line = 0, cex=trace.options$cex.axis-.1,
                                   padj = -(trace.options$cex.axis)/0.225)

      ### AEXES
      if (trace.box == TRUE) axis(2, cex.axis=trace.options$cex.axis)
      if (i == new.nvar & (trace.axis[1] == TRUE | is.numeric(trace.axis))) {
        axis(1, at = breaks, las=1, outer = TRUE)
      }
      if (trace.box == TRUE) box()
      ### DENSITY PLOT
      par(mar=c(0,.5,0,0))
      hdata.list <- c(list(x=as.vector(vmat), plot=FALSE), hdata.options)
      hdata <- do.call(hist, hdata.list)
      hist.list <- c(list(height=hdata$counts, col = tail(col.list[[i]],1)), bar.options)
      do.call(barplot, hist.list)
    }
    if (!missing(title)) {
      title.list <- c(list(text = title), title.options)
      do.call(mtext, title.list)
    }
  }
  }
}

