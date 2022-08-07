

######################################################################
############ Identification of regulatory sites - signal matrix
######################################################################


#' # Get a score in each window
#'
#' @param target a path to a tabix file
#' @param windows a GRanges to be overlapped with ranges in target
#' @param strand.aware if TRUE (default: FALSE), the strands of the windows
#'                     will be taken into account. If the strand of a window is -,
#'                     the values of the bins for that window will be reversed
#' @param save.db a logical to decide whether regional counts within given windows
#'                should be saved as flat file database or not
#' @param suffix a character string to append to the name of the output
#'               flat file database, only used if save.db is true
#' @param cov.bases number minimum bases covered per region (Default:0).
#'                  Only regions with base coverage above this threshold are returned.
#' @param chunk.size Number of rows to be taken as a chunk for processing the methylDB
#'                   objects (default: 1e6)
#'
#' @examples
#'
#' library(methylKit)
#' file.list=list( system.file("extdata", "test1.myCpG.txt", package = "methylKit"),
#'                system.file("extdata", "test2.myCpG.txt", package = "methylKit"),
#'                system.file("extdata", "control1.myCpG.txt", package = "methylKit"),
#'                system.file("extdata", "control2.myCpG.txt", package = "methylKit") )
#' myobj=methRead( file.list,
#'            sample.id=list("test1","test2","ctrl1","ctrl2"),assembly="hg18",
#'            pipeline="amp",treatment=c(1,1,0,0))
#' methylbase.obj=unite(myobj)
#' methylbase.objDB=makeMethylDB(methylbase.obj, dbdir = "methylDB")
#'
#' windows=GRanges("chr21",
#'	IRanges(start=c(9853296,9906604), end=c(9860126,9906616)))
#'
#'
#' sig.mtx=getSignalMatrix(
#'              Rsamtools::TabixFile(methylbase.objDB@dbpath),
#'	            windows,
#'	            methylbase.objDB@sample.ids,
#'	            methylbase.objDB@assembly,
#'	            methylbase.objDB@treatment)
#'
#' @return returns a \code{SignalMatrix} object that is
#' 		   a numeric matrix of scores such as DNA methylation values
#' 		   in each target region (rows) in each sample (columns) 
#' @docType methods
#' @export
setGeneric("getSignalMatrix",
           function(target,
                    windows,
                    sample.ids=NULL,
                    assembly="",
                    treatment=NULL,
                    strand.aware=FALSE,
                    save.db=FALSE,
                    suffix=NULL,
                    cov.bases=0,
                    chunk.size=10000,
                    context="CpG",
                    destranded=TRUE,
                    resolution="base")
           standardGeneric("getSignalMatrix") )


#' @aliases getSignalMatrix,TabixFile,GRanges,character,character-method
#' @rdname getSignalMatrix-methods
#' @usage  \S4method{getSignalMatrix}{TabixFile,GRanges,character,character,numeric}(target,windows,sample.ids,assembly,treatment)
setMethod("getSignalMatrix",signature("TabixFile","GRanges","character","character","numeric"),
          function(target,
                   windows,
                   sample.ids=NULL,
                   assembly="",
                   treatment=NULL,
                   strand.aware=FALSE,
                   save.db=FALSE,
                   suffix=NULL,
                   cov.bases=0,
                   chunk.size=10000,
                   context="CpG",
                   destranded=TRUE,
                   resolution="base"){

            methylBase.obj = methylKit:::readMethylBaseDB(
              dbpath = target$path,
              dbtype = "tabix",
              sample.ids=sample.ids,
              assembly=assembly,
              context=context,
              treatment=treatment,
              destranded=destranded,
              resolution=resolution)
            
            # regions.methylBase=selectByOverlap(methylBase.obj, windows)
            # regions.meth = percMethylation(regions.methylBase,
            #                                chunk.size=chunk.size)
            # Get percent of methylation within given windows
            # using a methylBase object compressed to a tabix file
            regions.methylBase=regionCounts(object = methylBase.obj,
                                            regions = windows,
                                            cov.bases = cov.bases,
                                            strand.aware=strand.aware,
                                            save.db = save.db,
                                            suffix=suffix,
                                            chunk.size=chunk.size)
            regions.meth = methylKit::percMethylation(regions.methylBase,
                                                      chunk.size=chunk.size)
            
            # Get average percent methylation within windows
            # 	meth.ave.fun = function(dt,numCs.index,numTs.index){
            # 	      dt[, numCs.index,with=FALSE]/
            # 	      (dt[,numCs.index,with=FALSE] +
            # 	         dt[,numTs.index,with=FALSE] )
            #   	}
            #   	windows_percmeth=methylKit:::applyTbxByChunk(
            # 						            tbxFile = regions.methylBase@dbpath,
            # 	                      chunk.size=chunk.siz,
            # 	                      return.type="data.table",
            # 	                      FUN = meth.ave.fun,
            #                         numCs.index=regions.methylBase@numCs.index,
            # 	                      numTs.index=regions.methylBase@numTs.index)
            #   	colnames(windows_percmeth) = methylBase.obj@sample.ids #######@TODO
            
            
            # rownames(regions.meth) = paste(regions.methylBase@.Data[[1]],
            #                                regions.methylBase@.Data[[2]],
            #                                regions.methylBase@.Data[[3]],
            #                                sep=".")
            
            return(regions.meth)

          })

######################################################################
## Filter input regulatory regions by their correlation
## with gene expression
######################################################################


#'
#' @param regions a GRanges object corresponding to
#' @param species a character indicating a genome version
#' @param rule a character indicating
#' @param plot a boolean value indicating whether to plot ...
#' @param ... a parameters passed to the GREAT::submitGreatJob function
#'
#'
.assoc.region2gene = function(regions,
                              species="hg38",
                              rule="basalPlusExt",
                              plot = FALSE,
                              ...){
  #require(rGREAT)
  mybed = data.frame(chr=seqnames(regions),
                     start=start(regions),
                     end=end(regions))
  job = submitGreatJob(mybed,
                       species = species,
                       rule=rule,
                       ...)
  #tb.rgreat = getEnrichmentTables(job)
  genes=plotRegionGeneAssociationGraphs(job,
                                        plot = plot)
  return(genes)
}


#' # Filter rows in SignalMatrix by the correlation with their nearby genes
#'
#' @param sig.mtx a SignalMatrix object or a numeric matrix
#' @param signal.granges a path to a gtf file
#' @param gene.matrix a numeric matrix of gene expression with the same dimensions as `signal`
#' @param corr.method a character string indicating which correlation coefficient
#'                   is to be used for the test. One of "pearson", "kendall", or
#'                   "spearman", can be abbreviated.
#' @param species species name supported by rGREAT such as "hg38", "hg19", "mm10", "mm9" etc.
#' @param rule a rule how to associate genomic regions to genes supported by RGREAT, such as
#'             "basalPlusExt", "twoClosest", "oneClosest"
#' @param plot a boolean value indicating whether rGREAT::plotRegionGeneAssociationGraphs()
#'             will generates plots
#' @param cores a number of cores
#'
#' @details
#'
#' @examples
#'
#' @return returns a \code{Model} object
#' @docType methods
#' @export
getCorrRegions = function(SignalMatrix,
                          signal.granges,
                          gene.matrix,
                          corr.method="spearman",
                          species="hg38",
                          rule="basalPlusExt",
                          plot = FALSE,
                          cores=1){
  
  assoc=.assoc.region2gene(signal.granges,
                           species=species,
                           rule=rule,
                           plot=plot)
  
  # make sure that columns in `region.cal` and `gene.expression`
  # refer to the same samples
  sample.ids=intersect(colnames(SignalMatrix),
                       colnames(gene.matrix))
  gene.expression = gene.matrix[,match(sample.ids,
                                       colnames(gene.matrix))]
  regions.val = SignalMatrix[,match(sample.ids,
                                    colnames(SignalMatrix))]
  
  
  # one DMR can be associated with many genes
  # from the previous analysis it looks like usually
  # with two genes
  if(is.null(assoc$gene)){
    stop("rGREAT::submitGreatJob didn't return column `gene`")
  }
  cor3=
    mclapply(1:length(signal.granges), function(i){
      region.i=signal.granges[i]
      assoc.i = subsetByOverlaps(assoc,region.i)
      regions.val.i = regions.val[i,]
      gene.expression.i = gene.expression[match(assoc.i$gene,
                                                rownames(gene.expression)),]
      if(is.data.frame(gene.expression.i)){
        my.cor=sapply(1:nrow(gene.expression.i), function(j){
          cor(regions.val.i,
              as.numeric(gene.expression.i[j,]),
              method=corr.method)
        })
        
      }else{
        my.cor=cor(regions.val.i,
                   gene.expression.i,
                   method=corr.method)
      }
      return(my.cor)
      
    }, mc.cores=cores)
  
  mcols(assoc)$cor = unlist(cor3)
  return(assoc)
  
}



#' # Get a score in each window
#'
#' @param target a path to a tabix file
#' @param windows a GRanges to be overlapped with ranges in target
#' @param strand.aware if TRUE (default: FALSE), the strands of the windows
#'                     will be taken into account. If the strand of a window is -,
#'                     the values of the bins for that window will be reversed
#' @param save.db a logical to decide whether regional counts within given windows
#'                should be saved as flat file database or not
#' @param suffix a character string to append to the name of the output
#'               flat file database, only used if save.db is true
#' @param cov.bases number minimum bases covered per region (Default:0).
#'                  Only regions with base coverage above this threshold are returned.
#' @param chunk.size Number of rows to be taken as a chunk for processing the methylDB
#'                   objects (default: 1e6)
#'
#' @examples
#'
#' library(methylKit)
#' file.list=list( system.file("extdata", "test1.myCpG.txt", package = "methylKit"),
#'                system.file("extdata", "test2.myCpG.txt", package = "methylKit"),
#'                system.file("extdata", "control1.myCpG.txt", package = "methylKit"),
#'                system.file("extdata", "control2.myCpG.txt", package = "methylKit") )
#' myobj=methRead( file.list,
#'            sample.id=list("test1","test2","ctrl1","ctrl2"),assembly="hg18",
#'            pipeline="amp",treatment=c(1,1,0,0))
#' methylbase.obj=unite(myobj)
#' methylbase.objDB=makeMethylDB(methylbase.obj, dbdir = "methylDB")
#'
#' windows=GRanges("chr21",
#'	IRanges(start=c(9853296,9906604), end=c(9860126,9906616)))
#'
#'
#' sig.mtx=getSignalMatrix(
#'              Rsamtools::TabixFile(methylbase.objDB@dbpath),
#'	            windows,
#'	            methylbase.objDB@sample.ids,
#'	            methylbase.objDB@assembly,
#'	            methylbase.objDB@treatment)
#'	            
#' gene.matrix =  matrix(10:17, 2,4)   
#' 
#' fil.gene.sign.mtx = 
#'    filterSignalMatrix(sig.mtx,
#'                       windows,
#'                       gene.matrix,
#'                       corr.method="spearman",
#'                       species="hg38",
#'                       rule="basalPlusExt",
#'                       plot = FALSE,
#'                       cores=1)        
#'
#' @return returns a list of \code{SignalMatrix} object that is
#' 		   a numeric matrix of scores such as DNA methylation values
#' 		   in each target region (rows) in each sample (columns) and
#' 		   a GRanges object
#' @docType methods
#' @export
filterSignalMatrix = function(SignalMatrix,
                              signal.granges,
                              gene.matrix,
                              abs.corr.threshold = 0.5,
                              corr.method="spearman",
                              species="hg38",
                              rule="basalPlusExt",
                              plot = FALSE,
                              cores=1){
  
  corr.obj=getCorrRegions(SignalMatrix,
                          signal.granges,
                          gene.matrix,
                          corr.method=corr.method,
                          species=species,
                          rule=rule,
                          plot = plot,
                          cores=cores)
  
  # filter signal regions by correlation with gene expression
  dmrs.assoc=unique(corr.obj[which(!is.na(corr.obj$cor)),])
  dmrs.assoc.corr=dmrs.assoc[which( dmrs.assoc$cor > abs.corr.threshold )]
  
  # filter SignalMatrix by rows according to correlation
  signal.granges$id=1:length(signal.granges)
  tmp=mergeByOverlaps(signal.granges,
                      dmrs.assoc.corr,
                      type="equal")
  org.dmrs.assoc.corr = tmp$`signal.granges`
  mcols(org.dmrs.assoc.corr) = cbind(mcols(org.dmrs.assoc.corr),
                                     mcols(tmp$`dmrs.assoc.corr`))
  SignalMatrix.corr=SignalMatrix[org.dmrs.assoc.corr$id,]
  
  return(list(SignalMatrix=SignalMatrix.corr,
              signal.granges=org.dmrs.assoc.corr))
  
}

