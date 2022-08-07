



######################################################################
############ Normalization
######################################################################

# Normalize MotifHitsMatrix
.normalize.x = function(Xi){
  Xi - matrix(rep(colMeans(Xi), each=dim(Xi)[1]), nrow=dim(Xi)[1])
}


# Normalize SignalMatrix
.normalize.y = function(Y,
                        transform.meth=FALSE){
  if(transform.meth){
    if(max(Y)<=1){
      Y = Y*100
    }
    #if(model=="multitask"){Y$methylation = 100-Y$methylation}
    Y <- log2((Y+1)/(100-Y+1))
  }
  # normalize signal matrix
  prep.y = function(N){
    N = N - matrix(rep(colMeans(N), each=dim(N)[1]), nrow=dim(N)[1])
    N = N - matrix(rep(rowMeans(N), each=dim(N)[2]), nrow=dim(N)[1], byrow=TRUE)
    N
  }
  prep.y(Y)
}


######################################################################
############ Modelling
######################################################################


.model.rank.lm = function(X,Y){
  d=as.data.frame(cbind(Y=Y, X))
  colnames(d) = c("Y", colnames(X))
  fit = Rfit::rfit(Y ~ ., data=d)
  coef=summary(fit)$coefficients[-1,]
  df=data.frame(coef=coef[,1], pval=coef[,4])
  rownames(df) = colnames(X)
  return(df)
}




.model.rank.lm.per.sample = function(X,Y,cores=1){


  # remove motifs that dont have any variance
  wh.0=which(apply(X, 1, var) == 0)
  if(sum(wh.0)>0){
    message("Removing motifs with variance equal to 0..")
    X = X[,-wh.0]
  }

  activity.per.sample = mclapply(1:ncol(Y),function(i){

    .model.rank.lm(X, Y[,i])

  },mc.cores=cores)



  activity.per.sample.coef <- data.frame(matrix(unlist(lapply(activity.per.sample,
                                                              function(x) x$coef)),
                                      ncol=length(activity.per.sample),
                                      byrow=FALSE))
  activity.per.sample.pval <- data.frame(matrix(unlist(lapply(activity.per.sample,
                                                              function(x) x$pval)),
                                                ncol=length(activity.per.sample),
                                                byrow=FALSE))
  colnames(activity.per.sample.coef) = colnames(Y)
  rownames(activity.per.sample.coef) = colnames(X)

  colnames(activity.per.sample.pval) = colnames(Y)
  rownames(activity.per.sample.pval) = colnames(X)

  return(list(coef=activity.per.sample.coef,
              pval=activity.per.sample.pval))
}


#' # Get a model
#'
#' @param signal a SignalMatrix object or a numeric matrix
#' @param motifhits a MotifHitsMatrix object or a numeric  matrix
#' @param cores a number of cores
#' @param transform.meth a boolean whether log transform input \code{signal}
#'                       using a following transformation: log2((Y+1)/(100-Y+1))
#'
#' @details
#' https://journal.r-project.org/archive/2012-2/RJournal_2012-2_Kloke+McKean.pdf
#'
#' @examples
#' 
#' signalmatrix_example_path = system.file("extdata", 
#'                                    "signalmatrix_example.txt", 
#'                                    package = "motifActivity")
#' signalmatrix_example = read.table(signalmatrix_example_path, header=TRUE)
#' 
#' motifhitsmatrix_example_path = system.file("extdata", 
#'                                        "motifhitsmatrix_example.txt", 
#'                                        package = "motifActivity")
#' motifhitsmatrix_example = read.table(motifhitsmatrix_example_path, header=TRUE)
#'
#' model=getModel(signalmatrix_example,
#'               motifhitsmatrix_example,
#'               cores=1,
#'               transform.meth=TRUE)
#'
#' @return returns a \code{Model} object
#' @docType methods
#' @export
getModel = function(signal,
                    motifhits,
                    cores=1,
                    transform.meth=FALSE){

  # prepare input data
  X=.normalize.x(motifhits)
  Y=.normalize.y(signal,transform.meth=transform.meth)

  # modeling
  mod=.model.rank.lm.per.sample(X, Y, cores)

  # performance of the model
  # 

  # output
  list(norm.signal = Y,
      norm.motif.hits = X,
      coef = mod$coef,
      pval = mod$pval)

}


######################################################################
############ Estimate which DMRs / genomic regions are targets
############ of motifs:
######################################################################



######################################################################
## . get diffferentially expressed genes
######################################################################



.GetSignificantDEGenes = function(DEresults,
                                  condition,
                                  case,
                                  control,
                                  rowsums=1,
                                  padj=0.1,
                                  log2foldchange=1){

  des <- DEresults[ rowSums(DESeq2::counts(DEresults)) > rowsums, ]
  res <- results(des, contrast=c(condition, case, control))
  #remove genes with NA values
  DE <- res[!is.na(res$padj),]
  #select genes with adjusted p-values below 0.1
  DE <- DE[DE$padj < padj,]
  #select genes with absolute log2 fold change above 1 (two-fold change)
  DE <- DE[abs(DE$log2FoldChange) > log2foldchange,]
  DE
}


#' # Perform differential gene expression analysis using DESeq2 in order to obtain
#' # active/expressed genes.
#'
#' It requires minimum of 2 biological replicates for each condition.
#'
#' @param gene.expression a numeric matrix indicating gene expression.
#'                        In columns are samples and in rows genes.
#' @param conditions a data.frame with a one column `condition` and row names
#'                   indicate names of samples the same as in column names
#'                   of the `gene.expression` matrix
#' @param case
#' @param control
#' @param covariates a data.frame that contains covariates for differential
#'                   expression analysis. Each column indicates a covariate
#'                   and in rows are samples.
#' @param rowsums
#' @param padj
#' @param log2foldchange
#'
#' @return returns
#' @docType methods
#' @export
getActiveGenes = function(gene.expression,
                          conditions,
                          case,
                          control,
                          covariates=NULL,
                          rowsums=1,
                          padj=0.1,
                          log2foldchange=1){

  colData = conditions
  colData$condition = as.factor(colData$condition)

  if(!is.null(covariates)){
    if(any(rownames(conditions)!=rownames(covariates))){
      stop("Names of rows of `conditions` and `covariates` are not the same")
    }
    colData = cbind(colData, covariates)
    my.formula = paste0("~condition+",
                        paste0(colnames(colData)[2:ncol(colData)],
                               collapse="+"))
  }else{
    my.formula = "~condition"
  }

  rownames(colData)<-NULL
  dds_covs <- DESeqDataSetFromMatrix(countData = gene.expression,
                                     colData = colData,
                                     design = as.formula(my.formula))
  dds_covs <- estimateSizeFactors(dds_covs)


  if(length(levels(colData$condition)) == 1){
    stop("There are not enough conditions of samples to compare")
  }
  if(length(levels(colData$condition)) == 2){
    # perform DESeq
    dds <- DESeq(dds_covs)
  }else{
    # perform DESeq ANOVA
    # https://support.bioconductor.org/p/61563/#61564
    if(is.null(covariates)){
       de_formula="~ 1"
     }else{
       de_formula=paste0("~",
                         paste0(colnames(colData)[2:ncol(colData)],
                                collapse="+"))
          }
    #dds = DESeq(dds_covs, test = "LRT", reduced = de_formula)
    dds = DESeq(dds_covs,
                test = "LRT",
                reduced = ~immune.score+microenvironment.score)

  }

  sde=.GetSignificantDEGenes(dds,
                         "condition",
                         case,
                         control,
                         rowsums=rowsums,
                         padj=padj,
                         log2foldchange=log2foldchange)
  list(DESeqDataSet = dds,
       significant.genes=sde)
}


















