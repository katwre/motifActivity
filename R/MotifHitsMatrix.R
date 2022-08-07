

##################################################################################
############ Identification of regulatory sites - matrix that contains motif hits
##################################################################################


#' # Get nucleotide sequences from a given windows
#'
#' @param windows a GRanges object
#' @param genome a BSgenome object
#'
getSequences = function(windows, genome,
                        background=FALSE, order=1,
                        cores=1){

  # trim sequences if needed
  seqinfo(windows) = seqinfo(genome)[seqnames(seqinfo(windows))]
  windows= GenomicRanges::trim(windows)

  # get nucleotide sequences
  ShortRead::clean( Biostrings::getSeq(genome, windows) )
}



#' # Calculate number of occurences of DNA motifs in given sequences
#'
#' @param pfms a list of PWM matrices
#' @param signal genomic regions either as a GRanges object.
#'               If a GRanges if provided then genome argument is required.
#' @param bg a BSgenome object
#' @param order order of the Markov models that will be used as the background model.
#'              Default: order = 1
#' @param normalize normalizes a PFM and optionally adds pseudo-evidence to each entry
#'                  of the matrix. It internally uses the motifcounter::normalizeMotif
#'                  function.
#' @param singlestranded
#' @param cores a numeric indicating number of cores. Default: 1.
#' @param genome a BSgenome object
#' @param ... parameters to pass to the motifcounter::motifcounterOptions function
#'
#' @examples
#'
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#' genes=genes(txdb)
#' genes=keepSeqlevels(genes, "chr12", pruning.mode="coarse")[1:50]
#' library(RCAS)
#' bg=createControlRegions(genes)
#' 
#' DNAmotifs.pfms <- as.list(query(query(MotifDb, "Hsapiens"),
#'                                 "jaspar2018"))[1:10]
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' genome=Hsapiens
#' 
#' # run finding motif hits in given windows
#' motifhits = getMotifsHits(DNAmotifs.pfms,
#'                           genes,
#'                           bg,
#'                           normalize=TRUE,
#'                           singlestranded=FALSE,
#'                           cores=1,
#'                           genome=genome,
#'                           order=0,
#'                           alpha=0.0001, gran=0.05)
#'
#'
#' @return returns a \code{MotifsHitsMatrix} object that is
#' 		   a matrix of occurrences DNA motifs in given genomic regions
#' @docType methods
#' @export
#'
setGeneric("getMotifsHits",
           function(pfms,
                    signal,
                    bg,
                    order=1,
                    normalize=TRUE,
                    singlestranded=FALSE,
                    cores=1,
                    genome=NULL,
                    ...)
             standardGeneric("getMotifsHits") )

#' @aliases getMotifsHits,list,GRanges-method
#' @rdname getMotifsHits-methods
#' @usage  \\S4method{getMotifsHits}{list,GRanges,"GRanges"}(pfms,signal,bg)
setMethod("getMotifsHits",signature("list","GRanges","GRanges"),
          function(pfms,
                   signal,
                   bg,
                   order=1,
                   normalize=TRUE,
                   singlestranded=FALSE,
                   cores=1,
                   genome,
                   ...){

            #require(motifcounter)
            motifcounterOptions() # let's use the default
            if(length(list(...))>0)
              motifcounterOptions(...) # let's use given parameters if any


            # make sure that motifs are in a format that in
            # rows correspond to nucleotides
            nrow.motifs=sapply(pfms, nrow)
            if(sum(nrow.motifs)!=(length(nrow.motifs)*4)){
              pfms = lapply(pfms, t)
            }


            # normalize DNA motifs
            if(normalize)
              pfms=lapply(pfms,function(x) normalizeMotif(x)) #normalizeMotif(t(x)))


            # get sequences if not provided
            if(!is.null(genome)){

              message("Getting signal sequences..")
              if(is(signal, "GRanges"))
                signal=getSequences(signal, genome, cores=cores)

              message("Getting background sequences..")
              if(is(bg, "GRanges"))
                bg=getSequences(bg, genome, cores=cores)

            }

            # estimate a background model from a set of DNA sequences
            message("Estimate a background model..")
            bg.object=readBackground(bg, order)

            # calculate number of hits of motifs in given input sequences
            message("Calculate number of motifs hits in given input sequences..")
            numMotifHits = mclapply(1:length(pfms),
                                    function(i){
                                      motifcounter:::numMotifHits(signal,
                                                                  pfms[[i]],
                                                                  bg.object,
                                                                  singlestranded = singlestranded)
                                    }, mc.cores=cores)

            # convert motif hits into a matrix
            res.numMotifHits <- lapply(numMotifHits, function(x) x$numofhits)
            output <- matrix(unlist(res.numMotifHits),
                             ncol = length(pfms),
                             byrow = TRUE)
            colnames(output) = names(pfms)

            return(output)

          })


