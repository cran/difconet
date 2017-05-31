% File man/difconet.run.Rd
\name{difconet.run}
\alias{difconet.run}
\alias{difconet}
\title{RUNS A DIFCONET ANALYSIS}
\description{
  Estimates the DIFferential COrrelation NETworks analysis from a given dataset.
}
\usage{
difconet.run(data, predictor, metric=c(1,2,3,4,5,6), cutoff=0.3, blocs=5000, 
  num_perms=10, comparisons="all", perm_mode="columns", use_all_perm = TRUE,
  save_perm=FALSE, speedup=0, verbose=TRUE, metricfunc=NULL, 
  corfunc=function(a,b) cor(a,b,method="spearman") )
}
\arguments{
  \item{data}{data.frame or matrix represent the dataset. Genes in rows, samples in columns. }
  \item{predictor}{Factor or numeric vector representing the classes of each column in data. The correlations will be estimated for each class separately.}
  \item{metric}{The metrics needed to be calculated. Valid values are 1 to 6 and 8. 1 to 6 are already implemented and shown in details. 8 specifies a user-defined metric specified in metricfunc. }
  \item{cutoff}{Cut off values used for metric 1 and/or 3.}
  \item{blocs}{Number of rows per block. Because of memory issues, the correlations are estimated by blocks of genes. This value represent the size of the block. Larger values requires more memory if needed. Lower values requiere more cycles and therefore it is slower but makes it computable depending on database size and memory.}
  \item{num_perms}{Number of permutations.}
  \item{comparisons}{Character or list. If character, it could be "all" to specify all possible combinations of classes. If set to "seq", classes are taken in order and comparisons are done by first versus second, second versus third, and so on. If this is a list containing vectors of two elements, the estimations are done for the specific comparisons included (numeric or character).}
  \item{perm_mode}{Character. It determines the how the permutated data is generated. It can be permutated by "columns",  permutated by "rows" (all classes/stages), or permutated by rows within each class separately using "rows.class", or "all" in which all data is shuffled. }
  \item{use_all_perm}{Logical. If TRUE, it uses all permutated data to estimate the p-value, otherwise it uses only the same row permutations to estimate the p-value (it requires a lot more permutations).}
  \item{save_perm}{Logical. If TRUE, it save all permutated data. It may require more memory.}
  \item{speedup}{Numeric. Determines whether the calculation will be sped up. This is experimental. The value specify which metric will be used to speed up. This is done by modeling the dependency of the metric and p-value using 1 percent of the rows.}
  \item{verbose}{Logical. Determines if printing progress information.}
  \item{metricfunc}{Function. Specify the function to be used if a metric==8 is included. The function should receive dObj, a, and b which correspond to the difconet object and the a and b vectors of correlations needed to estimate the value of the metric. It is assumed a distance-like measure (non-negative) and values close to 0 means no difference whereas larger values represent more dissimilar correlations.}
  \item{corfunc}{Function. Specify the function that estimates the correlations, similar to the cor function. The default uses cor and spearman coefficients.}
}
\details{
Run the whole process of estimation differences in correlations for a given dataset. The estimations are done for all \strong{metric} values, all \strong{cutoff} values across all \strong{comparisons}.
}
\value{
A difconet object represented as a list. The items are the followings:
\item{stage}{Vector. A copy of \strong{predictor} (classes).}
\item{labels}{Vector. The levels or values of the different classes. }
\item{comparisons}{The specified \strong{comparisons} parameter.}
\item{num_perms}{The specified number of permutations \strong{num_perms} parameter.}
\item{perm_mode}{The specified number of permutations \strong{perm_mode} parameter.}
\item{use_all_perm}{The specified number of permutations \strong{use_all_perm} parameter.}
\item{speedup}{The specified \strong{speedup} parameter.}
\item{verbose}{The specified \strong{verbose} parameter.}
\item{metricfunc}{The specified \strong{metricfunc} parameter.}
\item{combinations}{A data.frame of the combinations that were compared.}
\item{stages.data}{A list of datasets. This is only the original \strong{data} split by classes.}
\item{combstats}{A list of all comparisons made. Each element contains a matrix whose rows represent the genes and columns represent the results of all metrics (metric.dist : metric value, metric.p : p-value, metric.q : q-value, metric.expr.p : p-value of differential expression for comparison purposes, metric.expr.q : q-value of differential expression.)}
\item{combdens}{A list of the densities of the metric for observed data and permutations. This can be used to compare the estimated metric statistics.}
\item{permutations}{List. If \strong{save_perm==TRUE}, it saves all permutated data.}
}
\references{
	Gonzalez-Valbuena and Trevino 2017 Metrics to Estimate Differential Co-Expression Networks \emph{Journal Pending} \bold{volume} 00--10
	}
\author{Elpidio Gonzalez and Victor Trevino \email{vtrevino@itesm.mx}}
\seealso{
  \code{\link{difconet.build.controlled.dataset}}.
}
\examples{

xdata <- matrix(rnorm(1000), ncol=100)
xpredictor <- sample(c("A","B","C","D"),100,replace=TRUE)
dObj <- difconet.run(xdata, xpredictor, metric = 4, num_perms = 10,              
  comparisons = list(c("A","D"), c("A","B"), c("B","D")),
  perm_mode = "columns")

\dontrun{
  #xpredictor contains A, B, C, and D.
  #xdata contains the data matrix
  dObj <- difconet.run(xdata, xpredictor,
  metric = c(1,2,4),
  cutoff = 0.6,
  blocs = 7000,
  num_perms = 10,              
  comparisons = list(c("A","D"), c("A","B"), c("B","D")),          
  perm_mode = "columns")
}

}
