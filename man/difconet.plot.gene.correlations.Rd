% File man/difconet.plot.gene.correlations.Rd
\name{difconet.plot.gene.correlations}
\alias{difconet.plot.gene.correlations}
\title{PLOTS THE CORRELATIONS OF A SPECIFIC GENE}
\description{
  Draw scatter plots of the correlations of a specific gene.
}
\usage{
difconet.plot.gene.correlations(dObj, gene, 
  stages=1:length(dObj$stages.data), type=c("density","scatter")[1], 
  main=rownames(dObj$stages.data[[1]])[gene], 
  legends=TRUE, plot=TRUE, ... )
}
\arguments{
  \item{dObj}{The difconet object.}
  \item{gene}{Numeric or character. The gene index/rowname whose correlations will be drawn.}
  \item{stages}{Numeric or character. The stages to be included. If \strong{type="scatter"} and more than two stages, a call to pairs is used instead of plot.}
  \item{type}{Character. The type of plot \strong{density} or \strong{scatter}.}
  \item{main}{Character. The main title passed to plot.}
  \item{legends}{Logical. Specifies whether the legends are drawn when \strong{type="density"}.}
  \item{plot}{Logical. Specifies whether the plots are actually drawn (to get the correlations).}
  \item{...}{Further parameters passed to plot/pairs.}
}
\details{
Run the whole process of estimation differences in correlations for a given dataset. The estimations are done for all \strong{metric} values, all \strong{cutoff} values across all \strong{comparisons}.
}
\value{
The correlations of the gene across stages (invisible).
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
#Another specific gene (3), showing densities of correlations
difconet.plot.gene.correlations(dObj, 3, type="d")

}
