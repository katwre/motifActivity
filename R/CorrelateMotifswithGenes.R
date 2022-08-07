

##################################################################################
############ Filter motifs activities by their correlation with expression
############ of their TFs
##################################################################################


#' # Filter motifs in activity matrix by their correlation with their TF expression
#'
#' @param coef a numeric matrix indicating an activity motif matrix
#' @param gene.expression a numeric matrix indicating a gene expression matrix
#'                        that is used to correlate with a `coef` entries
#' @param threshold a numeric value indicating correlation threshold
#' @param method a character indicating a method used for correlation.
#'               By default "pearson", for more details check out
#'               the `method` argument of the `stats:cor` function.
#' @param cores a numeric indicating number of cores. Default: 1.
#' @param abs a boolean value indicating whether calculate absolute value of
#'            correlation. Defalt: FALSE.
#' @param split.sign a character value
#'
#' @return returns a \code{gtable} object from a ggplot R package
#' @docType methods
#' @export
filterActivityByTFexpression =
  function(coef,
           gene.expression,
           threshold=0.6,
           method="spearman", # spearman for methylation data
           cores=1,
           #fun=function(x) max(abs(x)),
           abs=FALSE,
           split.sign=NULL
  ){
 
    
    #require(data.table)
    # Match names of samples of coef and gene.expression
    # matrices (columns)
    #if(colnames(gene.expression)==colnames(coef))
    #   warning(paste0("columns names `gene.expression` are not ",
    #               "the same as column names of `coef`"))
    inter.samples = intersect(colnames(coef),
                              colnames(gene.expression))
    ge.matched=gene.expression[,match(inter.samples,
                                      colnames(gene.expression))]
    coef.matched=coef[,match(inter.samples,
                             colnames(coef))]


    # Remove rows that variation is equal to 0 -> if I do it then I
    # wil change the order of rows more likely and mess up
    # correlation part....
    #ge.matched = ge.matched[apply(ge.matched, 1, var) != 0, ]
    #coef.matched = coef.matched[apply(coef.matched, 1, var) != 0, ]


    ## Match names of motifs of coef and gene.expression
    ## matrices (rows)
    # some motifs might have few TFs that bind to them
    # if at least such 1 TF is correlated with TF expression
    # above given threshold then report it

    corr.motifcoef.tfexpr =
      mclapply(1:nrow(coef.matched),
               function(i){
        # get activity of an i-th motif across samples
        activity.motif.i = coef.matched[i,]
        # get a name of a TF corresponding to an i-th motif
        if(is.null(split.sign)){
          motif.TF.id = rownames(coef.matched)[i]
        }else{
          motif.TF.id = strsplit(rownames(coef.matched)[i],
                                 split.sign)[[1]]
        }

        # get gene expression of a TF corresponding to
        # an i-th motif across samples
        gmtch=match(motif.TF.id, rownames(ge.matched))
        if(any(is.na(gmtch))){return(NULL)}
        TFexp.motif.i = ge.matched[gmtch,]

        # if a DNA motif has 1 corresponding TF
        if(length(motif.TF.id )==1){
          my.cor=cor(as.numeric(activity.motif.i),
                     as.numeric(TFexp.motif.i),
                     method=method)
        }

        # if a DNA motif has more than 1 corresponding TF
        if(length(motif.TF.id )>1){
          my.cor=sapply(1:nrow(TFexp.motif.i),
                        function(k){
                          cor(activity.motif.i,
                              as.numeric(TFexp.motif.i[k,]),
                              method=method)
                        })
        }


        if(!all(is.na(my.cor))){
          if(any(is.na(my.cor))){
            my.cor.nona=my.cor[-which(is.na(my.cor))]
          }else{
            my.cor.nona=my.cor
          }
        if(any(abs(my.cor.nona)>=threshold)){

            data.frame(index=i,
                       cor=my.cor.nona, ### if len(my.cor)>1 then `i`-th row
                       ### will be repeated in the output `coef` var.
                       TF=motif.TF.id,
                       TF_initial=rownames(coef.matched)[i])

        }}

      }, mc.cores=cores)

    # Remove NULL in a list
    corr.motifcoef.tfexpr =
      corr.motifcoef.tfexpr[lengths(corr.motifcoef.tfexpr) != 0]

    # Merge results
    df.cor=rbindlist(corr.motifcoef.tfexpr)

    return(
      list(coef=coef[df.cor$index,],
           cor=df.cor$cor,
           TF=df.cor$TF,
           TF_initial=df.cor$TF_initial,
           index_initial=df.cor$index)
    )

  }





