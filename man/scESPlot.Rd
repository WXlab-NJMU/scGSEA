\name{scESPlot}
\alias{scESPlot}
\title{
Plot showing (normalized) enrichment scores of a specified gene set for each cell in 
low-dimention, in a single-cell experiment 
}
\description{
The function is to generate a scatter plot of enrichment scores of a specified gene set 
for each cell in low-dimention, in a single-cell experiment. 
}
\usage{
scESPlot(sc.obj, 
  gs.obj, 
  gene.set=gs.obj@GSNames, 
  cols.use = c("lightgrey", "blue"), 
  reduction.use = "tsne", ...) 
}
\arguments{
  \item{sc.obj}{}
  \item{gs.obj}{}
}
\author{
Wei Sheng, wsheng@stu.njmu.edu.cn
}
\seealso{
\code{\link{GSEnrichAnalyze}},
\code{\link{scESPlotGeneSet}}
}
\examples{
}
