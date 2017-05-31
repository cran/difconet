% File man/difconet.noise.inspection.Rd
\name{difconet.noise.inspection}
\alias{difconet.noise.inspection}
\title{PLOT ESTIMATED CORRELATION DISTRIBUTION AFTER ADDING NOISE}
\description{
  Plots the estimated correlation distribution of a normal dataset after adding different levels of gaussian noise. It is used to estimate the level of noise needed to be added to a normal dataset to match the correlation distribution of a tumor dataset. This assumes that the correlation distribution of the tumor dataset is sharper around zero.
}
\usage{
difconet.noise.inspection(ndata, tdata, sigma=c(0.5, 0.75, 1.25), maxgenes=5000, 
  corfunc=function(a,b) cor(a,b,method="spearman"))
}
\arguments{
  \item{ndata}{The normal dataset. Rows are genes and columns are samples.}
  \item{tdata}{The tumor dataset. Rows are genes and columns are samples. Rows of tumor and normal datasets should be the same.}
  \item{sigma}{Levels of gaussian noise to be added (at zero mean).}
  \item{maxgenes}{Number of genes used to estimate the correlation distribution. If the number of rows in normal/tumor datasets are larger than maxgenes, maxgenes random genes are used for the estimation.}
  \item{corfunc}{Correlation method used.}

}
\details{
Plots the estimated density of correlation distributions of normal, tumor, and normal after adding sigma levels of noise.
}
\value{
Nothing.
}
\references{
	Gonzalez-Valbuena and Trevino 2017 Metrics to Estimate Differential Co-Expression Networks \emph{Journal Pending} \bold{volume} 00--10
	}
\author{Elpidio Gonzalez and Victor Trevino \email{vtrevino@itesm.mx}}
\seealso{
  \code{\link{difconet.build.controlled.dataset}}.
  \code{\link{difconet.run}}.
}
\examples{

\dontrun{difconet.noise.inspection(normaldata, tumordata, sigma=0:15/10)}

}
