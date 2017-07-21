rm(list=ls())

requie("limma")
requie("edgeR")
requie("statmod")
requie("DESeq")


limma_voom_DE <- function(infile,a,b,method){
  if(missing(infile)){stop("input file data required.\n")}
  if(missing(a)){stop("Give number of replicate for first set of condition.\n ")}
  if(missing(b)){stop("Give number of replicate for second set of condition.\n")}
  if(missing(method)){stop("Give FDR as method.\n")}
  
  
file = basename(infile)

data = read.table(infile, header = T,row.names = 1) # Reading use input filename

group = as.factor(c(rep("x",a),rep("y",b))) # Giving number of replicate information for each condition.

head(data)

nf <- calcNormFactors(data,method = "TMM") # TMM Normalization
design <- model.matrix(~group) # Base on conditions and plicate information designing the info for your count matrix
design

y <- voom(data, design, plot=TRUE, lib.size=colSums(data)*nf) # voom function

fit <- lmFit(y,design) # Fitting to linear model for voom output
fit <- eBayes(fit)
toptable = topTable(fit,coef=ncol(design),number=20000,adjust.method="BH") # Getting toptable significant hits
#toptable = topTable(fit,coef=1,number=20000,genelist=fit$genes,adjust.method="fdr",sort.by="B",resort.by=NULL,p.value=0.05,lfc=0)
#toptable = topTableF(fit,number=10,genelist=fit$genes,adjust.method="BH",sort.by="F",p.value=1)
#toptable = topTable(fit,coef=NULL,adjust.method="fdr",p.value=0.05)
p.adjusted <- p.adjust(fit$p.value[,2], method=method) # Getting q value(p-adjusted value)

# Results Writing to file
toptable_out = paste0(Sys.Date(),file,"_toptable_rsults.txt")
write.table(toptable,toptable_out,quote = FALSE,row.names = FALSE)

results_limma <- cbind(fit$coeff, fit$p.value[,2], p.adjusted)
colnames(results_limma) <- c("av_expr", "2LogFC", "pvalue", "adjusted_pvalue")
results_limma <- results_limma[order(p.adjusted),]

results_limma_genes=cbind(rownames(results_limma),combat); dim(results_limma_genes)
colnames(results_limma_genes) <- c("Genes",colnames(results_limma))
output = paste0(Sys.Date(),file,"_","DE_results.txt")
write.table(results_limma,output,quote = FALSE,row.names = FALSE)

}
