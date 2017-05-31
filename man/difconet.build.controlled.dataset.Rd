% File man/difconet.build.controlled.dataset.Rd
\name{difconet.build.controlled.dataset}
\alias{difconet.build.controlled.dataset}
\title{GENERATES A DATASET CONTROLLING FOR NOISE AND GENES CONNECTED IN NETWORKS}
\description{
  This function takes a normal dataset and generate simulated tumor stages by adding progressive levels of noise. It may add artificial networks of genes connected at given correlations that can progressively increase or decrease their level of correlation.
}
\usage{
difconet.build.controlled.dataset(data,
    noise.genes = round(nrow(data)*0.1),
    noise.sigma = c(0.0, 0.1, 0.2), 
    nonoise.sigma = c(0.0, 0.01, 0.01), 
    netcov = matrix(c(
      0.90, 0.90, 0.75, 0.75, 0.60, 0.60, 0.45, 0.45, 0.30, 0.30, 
      0.15, 0.15, 0.30, 0.30, 0.45, 0.45, 0.60, 0.60, 0.75, 0.75,
      0.95, 0.95, 0.80, 0.80, 0.65, 0.65, 0.50, 0.50, 0.35, 0.35, 
      0.10, 0.10, 0.25, 0.25, 0.40, 0.40, 0.55, 0.55, 0.70, 0.70,
      1.00, 1.00, 0.85, 0.85, 0.70, 0.70, 0.55, 0.55, 0.40, 0.40, 
      0.05, 0.05, 0.20, 0.20, 0.35, 0.35, 0.50, 0.50, 0.65, 0.65
      ), ncol=3),
    genes.nets = 10,
    corfunc=function(a,b) cor(a,b,method="spearman"),
    verbose = TRUE)
}
\arguments{
  \item{data}{data.frame or matrix representing the normal dataset. Rows are genes and columns are samples.}
  \item{noise.genes}{the number of genes from data that will noised.}
  \item{noise.sigma}{Levels of gaussian noise to be added (at zero mean) expressed in a cumulative manner. }
  \item{nonoise.sigma}{Levels of gaussian noise to be added (at zero mean) for the rest of the genes.}
  \item{netcov}{numeric matrix of correlation levels for networks, rows represent networks and columns represent stages.}
  \item{genes.nets}{The number of genes in each generated network.}
  \item{corfunc}{Correlation method used.}
  \item{verbose}{Print progress.}
}
\details{
This function generates a simulated tumor progression dataset based on normal data. The progression is done by stages. The number of stages is given by the length of noise.sigma. Each stage will have the same dimensions than data (plus the networks). The stages will be N, T1, T2, and so on. The N is meant to be the data itself with no noise but for generality, the first element of noise.sigma specifies the level of noise for N (default to 0). The next values of noise.sigma will be used to generate T1, T2, and so on. Thus the returned data will be estimated by N=data+noise.sigma[1], T1=N+noise.sigma[2], T2=T1+noise.sigma[3], and so on. Note that noise.sigma will be added only to a specific number of rows given by noise.genes. The value returned is a list of the generated matrices. In top of that, the nonoise.sigma specify the level of noise added to those genes not selected to be noised. This is meant to be lower levels of noise than noise.sigma to avoid that data in stages is just a copy of previous data. This function also adds full connected networks of genes connected at netcov levels. The data added has mean=0 and sd=1. The number of rows represent the networks added. The columns represent the stages. 
}
\value{
List of stages.
}
\references{
	Gonzalez-Valbuena and Trevino 2017 Metrics to Estimate Differential Co-Expression Networks \emph{Journal Pending} \bold{volume} 00--10
	}
\author{Elpidio Gonzalez and Victor Trevino \email{vtrevino@itesm.mx}}
\seealso{
  \code{\link{difconet.noise.inspection}}.
  \code{\link{difconet.run}}.
}
\examples{

\dontrun{difconet.noise.inspection(normaldata, tumordata, sigma=0:15/10)}

}
