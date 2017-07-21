\name{limma_voom_DE}
\alias{limma_voom_DE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Differential Expression of count data using limma-voom
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
limma_voom_DE(infile, a, b, method)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{infile}{
Input file should be in order of condition(first set will be one condition and next condition will be next set)
}
  \item{a}{
Number of replicates for first set of condition
}
  \item{b}{
Number of replicate for second set of condition
}
  \item{method}{
FDR
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (infile, a, b, method) 
{
    file = basename(infile)
    data = read.table(infile, header = T, row.names = 1)
    group = as.factor(c(rep("x", a), rep("y", b)))
    head(data)
    nf <- calcNormFactors(data, method = "TMM")
    design <- model.matrix(~group)
    design
    y <- voom(data, design, plot = TRUE, lib.size = colSums(data) * 
        nf)
    fit <- lmFit(y, design)
    fit <- eBayes(fit)
    toptable = topTable(fit, coef = ncol(design), number = 20000, 
        adjust.method = "BH")
    p.adjusted <- p.adjust(fit$p.value[, 2], method = method)
    toptable_out = paste0(Sys.Date(), file, "_toptable_rsults.txt")
    write.table(toptable, toptable_out, quote = FALSE, row.names = FALSE)
    results_limma <- cbind(fit$coeff, fit$p.value[, 2], p.adjusted)
    colnames(results_limma) <- c("av_expr", "2LogFC", "pvalue", 
        "adjusted_pvalue")
    results_limma <- results_limma[order(p.adjusted), ]
    results_limma_genes = cbind(rownames(results_limma), combat)
    dim(results_limma_genes)
    colnames(results_limma_genes) <- c("Genes", colnames(results_limma))
    output = paste0(Sys.Date(), file, "_", "DE_results.txt")
    write.table(results_limma, output, quote = FALSE, row.names = FALSE)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
