% File man/difconet.plot.histograms.heatmap2.Rd
\name{difconet.plot.histograms.heatmap2}
\alias{difconet.plot.histograms.heatmap2}
\title{PLOT A HEATMAP REPRESENTATION OF THE DISTRIBUTION OF CORRELATIONS OF MANY GENES}
\description{
  Draw a heatmap whose rows are genes and columns are segments of the histogram of the distribution of correlations per gene. The height/density of the histogram is shown in colors.
}
\usage{
difconet.plot.histograms.heatmap2(dObj, 
  genes=1:10, 
  stages=1:length(dObj$stages.data), 
  qprobs=c(0,.50,.975,.995), ...)
}
\arguments{
  \item{dObj}{The difconet object.}
  \item{genes}{Numeric or character. The gene indexes/rownames included.}
  \item{stages}{Numeric or character. The stages to be included.}
  \item{qprobs}{The quantiles used to draw the heatmap. Should be 4 points. Each has specific color codes. }
  \item{...}{Further parameters passed to plot/pairs.}
}
\details{
A heatmap is draw representing the distribution of correlations of several genes across stages.
}
\value{
Nothing.
}
\references{
	Gonzalez-Valbuena and Trevino 2017 Metrics to Estimate Differential Co-Expression Networks \emph{Journal Pending} \bold{volume} 00--10
	}
\author{Elpidio Gonzalez and Victor Trevino \email{vtrevino@itesm.mx}}
\seealso{
  \code{\link{difconet.run}}.
}
\examples{

xdata <- matrix(rnorm(1000), ncol=100)
xpredictor <- sample(c("A","B","C","D"),100,replace=TRUE)
dObj <- difconet.run(xdata, xpredictor, metric = 4, num_perms = 10,              
  comparisons = list(c("A","D"), c("A","B"), c("B","D")),
  perm_mode = "columns")

  #Top highest metric in first comparison but showing correlations in only 3 stages
  difconet.plot.gene.correlations(dObj, order(dObj$combstats[[1]][,"M4.dist"], 
    decreasing=TRUE)[1], type="s", stages=1:3)
  #Bottom lowest metric in second comparison showing all stages
  difconet.plot.gene.correlations(dObj, order(dObj$combstats[[2]][,"M4.dist"], 
    decreasing=TRUE)[1], type="d")
  #Another specific gene (1), showing densities of correlations
  difconet.plot.gene.correlations(dObj, 1, type="d")
}
