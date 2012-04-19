\name{superdiag}
\alias{superdiag}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Test for Markov Chain Nonconvergence
%%  ~~function to do ... ~~
}
\description{
	The \code{superdiag} function takes MCMC samples as input. It provides a comprehensive test suite for Markov chain nonconvergence, which integrates four standard empirical MCMC convergence diagnostics: Gelman-Rubin, Geweke, Heidelberger-Welch, and Raftery-Lewis.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
superdiag(mcmcoutput, burnin = 10000, confidence.gr = 0.95, frac1.gw = 0.1, 
frac2.gw = 0.5, eps.hw = 0.1, pvalue.hw = 0.05, q.rl = 0.025, r.rl = 0.005, 
s.rl = 0.95, eps.rl = 0.001)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mcmcoutput}{input chains from \code{JAGS}, \code{BUGS}, etc.
%%     ~~Describe \code{mcmcoutput} here~~
}
  \item{burnin}{The number of burn-in iterations for the sampler.
%%     ~~Describe \code{burnin} here~~
}
  \item{confidence.gr}{(1-\eqn{\alpha})\% for the Gelman and Rubin test. The upper 95\% credible interval is the default.
%%     ~~Describe \code{confidence.gr} here~~
}
  \item{frac1.gw}{frac1 for the Geweke test. The proportion of the early era of the chain, defaulted to 0.1.
%%     ~~Describe \code{frac1.gw} here~~
}
  \item{frac2.gw}{frac2 for the Geweke test. The proportion of the late era of the chain, defaulted to 0.5.
%%     ~~Describe \code{frac2.gw} here~~
}
  \item{eps.hw}{\code{epsilon} for the Heidelberger and Welch test. The accuracy parameter determines whether the halfwidth is passed or not, defaulted to 0.1.
%%     ~~Describe \code{eps.hw} here~~
}
  \item{pvalue.hw}{p-value for the Heidelberger and Welch test. The halfwidth of the test calculates a (1-\eqn{\alpha})\% credible interval around the sample mean for each parameter dimension. The default is 0.05.
%%     ~~Describe \code{pvalue.hw} here~~
}
  \item{q.rl}{q-parameter for the Raftery and Lewis test. The posterior tail threshold of interest, defaulted to 0.025.
%%     ~~Describe \code{q.rl} here~~
}
  \item{r.rl}{r-parameter for the Raftery and Lewis test. The tolerance for the tail threshold, defaulted to 0.0005.
%%     ~~Describe \code{r.rl} here~~
}
  \item{s.rl}{s-parameter for the Raftery and Lewis test. The desired probability of being within the tolerance, defaulted to 0.95.
%%     ~~Describe \code{s.rl} here~~
}
  \item{eps.rl}{convergence epsilon for the Raftery and Lewis test. The convergence tolerance value, which is used to determine a stopping point based on a parallel chain process, defaulted to 0.001.
%%     ~~Describe \code{eps.rl} here~~
}
}
\details{
	If only one chain is analyzed, the default settings in \code{boa} and \code{coda} for Geweke, Raftery-Lewis, and Heidelberger-Welch diagnostics are used. If multiple chains are provided, only the first chain uses the defaults and all other chain analyses get random values as follows. For Geweke test, random non-overlapping proportions up from the start of the chain and down from the end of the chain are generated. For Heidelberger-Welch test, the value of \code{pvalue.hw} is sampled with replacement from common \eqn{\alpha} values; the value of \code{eps.hw} is sampled uniformly in the interval [0.01:0.2]. For Raftery-Lewis test, each of these four parameters are sample from a vector (changeable by users) of values around the defaults (larger and smaller) to provide a reasonable range of alternatives.
%%  ~~ If necessary, more details than the description above ~~
}

\references{
	\href{http://polmeth.wustl.edu/methodologist/tpm_v19_n2.pdf}{Tsai, Tsung-han and Gill, Jeff (2012). \dQuote{superdiag: A Comprehensive Test Suite for Markov Chain Non-Convergence.} \emph{The Political Methodologist}, 19 (2), 12-18.}
%% ~put references to the literature/web site here ~
}

\examples{
\dontrun{
{
data(tobit.list)
summary(tobit.list[1])
superdiag(tobit.list, burnin=0)
}
}	
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ MCMC }

