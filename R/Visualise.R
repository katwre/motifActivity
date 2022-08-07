
##################################################################################
############ Plot motifs activities
##################################################################################


#' # Plot motif activity for a given DNA motif
#'
#' @param motif.name a character indicating a name of a DNA motif
#' @param coef a matrix of coefficients. In columns are samples, in rows motifs.
#' @param title a string indicating a title of the plot. By default it's `motif.name`.
#' @param colour.line a character indicating a color of a line
#' @param colour.xaxis a character indicating a color of labels on X axis
#' @param legend.position a character indicating a position of legend, such as 
#'                        "left", "top", "right", or "bottom".
#'
#' @example 
#' coefficients_example_path = system.file("extdata", 
#'                                    "coefficients_example.txt", 
#'                                    package = "motifActivity")
#' coefficients_example = read.table(coefficients_example_path, 
#'                                   header=TRUE)
#' plotMotifActivity("MA0006.1_Ahr..Arnt",
#'                   coefficients_example)
#'
#' @return returns a \code{gg} object from a ggplot R package
#' @docType methods
#' @export
plotMotifActivity = function(motif.name,
                             coef,
                             title=motif.name,
                             colour.line="black",
                             colour.xaxis="black",
                             ylab="Motif activity",
                             match.clustering=TRUE,
                             legend.position="top"){
  #require(ggplot2)
  #require(reshape2)

  coef.m=melt(t(coef))
  # get coefficients for a given motif
  coef.m.motif=coef.m[which(as.character(coef.m$Var2)==motif.name),]

  #print(legend.position)
  g=ggplot(data=coef.m.motif, aes(Var1, round(value,5), group=1))+
    geom_point() +
    geom_line(#linetype = "dashed",
              colour=colour.line) +
    ylab(ylab)+xlab("Samples")+
    ggtitle(title)+
    theme(axis.text.x=element_text(angle=45, hjust=1,
                                   colour = colour.xaxis))#+
    #scale_colour_manual(values=unique(colour.line))+
    #scale_color_manual("Line.Color", values=c(red="red",green="green",blue="blue"),
    #                   labels=paste0("Int",1:3))
  return(g)

}



#' # Plot multiple motifs activities
#'
#' @param motif.name a character indicating a name of a DNA motif
#' @param coef a matrix of coefficients. In columns are samples, in rows motifs.
#' @param title a character vector indicating a title of the plot
#' @param colour.line a character indicating a colour of a line
#' @param colour.xaxis a character indicating a colour of labels on X axis
#' @param legend.position a character indicating a position of legend, such as 
#'                        "left", "top", "right", or "bottom".
#'
#' @example 
#' coefficients_example_path = system.file("extdata", 
#'                                    "coefficients_example.txt", 
#'                                    package = "motifActivity")
#' coefficients_example = read.table(coefficients_example_path, 
#'                                   header=TRUE)
#' plotMotifActivity(c("MA0006.1_Ahr..Arnt",
#'                      "MA0854.1_Alx1"),
#'                   coefficients_example)
#'                   
#' @return returns a \code{gtable} object from a ggplot R package
#' @docType methods
#' @export
#'
plotMotifsActivity = function(motif.names,
                              coef,
                              title=motif.names,
                              colour.line="black",
                              colour.xaxis="black",
                              ylab="Motif activity",
                              legend.position="top",
                              ...){
  #require(gridExtra)
  #require(ggplot2)
  plist = lapply(1:length(motif.names),
                 function(i){
                   x=motif.names[i]
                   plotMotifActivity(x,
                                     coef, #@TODO coef.m=melt(t(coef)) -> time consuming
                                     title[i],
                                     colour.line=colour.line,
                                     colour.xaxis=colour.xaxis,
                                     ylab=ylab,
                                     legend.position=legend.position)
                 })
  grid.arrange( grobs = plist,
                ...)
}







##################################################################################
############ Plot interactions
##################################################################################


plotSTRINGDBNetwork = function(genes,
                               version,
                               species=9606 # human
                               ){
  #require(STRINGdb)
  string_db <- STRINGdb$new( version=version,
                             species=species,
                             score_threshold=200,
                             input_directory="")
  string.mna<-string_db$map( data.frame(gene=genes),
                             "gene",
                             removeUnmappedRows = TRUE)
  hits <- string.mna$STRING_id
  string_db$plot_network( hits )

}


#' # Plot transcription factors in a PPI network build from the STRING database
#'
#' @param genes.hgnc a character vector indicating names of genes in the HUGO Gene Nomenclature
#' @param species a number corresponding to a species required for the STRINGdb$new function
#' @param version a number indicating a version of the STRINGdb database
#' @param score_threshold a number indicating a STRING combined score threshold 
#'                        (the network loaded contains only interactions having a combined 
#'                        score greater than this threshold)
#' @param nodes.color a character vector indicating color of nodes
#' @param mark.color a character vector indicating color of nodes
#' @param cluster.nodes a boolean value indicating whether to cluster nodes and color-code them
#'                      together
#' @param plot_dendrogram a boolean value indicating whether to plot dendogram of clustered
#'                        nodes
#' @param cut_at a numeric value indicating a cut off for number of nodes
#' @param graph_sub a graph object from the igraph R package indicating graph of TFs
#'
#' @return
#' @docType methods
#' @export
plotNetwork = function(genes.hgnc,
                       species=9606,
                       version="11",
                       score_threshold=200,
                       input_directory="./",

                       nodes.color = NULL,
                       mark.color=NULL,

                       cluster.nodes=TRUE,
                       plot_dendrogram=FALSE,
                       cut_at=NULL,
                       
                       add.missing=FALSE,

                       graph_sub=NULL,
                       ...){

  #require(STRINGdb)
  #require(igraph)
   
  if(!is.null(graph_sub)){
    human_graph_sub = graph_sub
    # get sub-network with only target genes
    string_db <- STRINGdb$new(species=species,
                              version=version,
                              score_threshold=score_threshold,
                              input_directory=input_directory)
    string_db_targetgenes<-string_db$map(
      data.frame(gene=genes.hgnc),
      "gene",
      removeUnmappedRows = TRUE)
    human_graph_sub=string_db$get_subnetwork(string_db_targetgenes$STRING_id)
    

  }else{
    string_db <- STRINGdb$new(species=species,
                              version=version,
                              score_threshold=score_threshold,
                              input_directory=input_directory)
    human_graph <- string_db$get_graph()

    # get sub-network with only target genes
    string_db_targetgenes<-string_db$map(
      data.frame(gene=genes.hgnc),
      "gene",
      removeUnmappedRows = TRUE)
    human_graph_sub=string_db$get_subnetwork(string_db_targetgenes$STRING_id)
  }

  # replace STRINGdb ids of genes to HGNC ids
  net <- set.vertex.attribute(human_graph_sub,
                              "name",
                              value=string_db_targetgenes[
                                match(V(human_graph_sub)$name,
                                      string_db_targetgenes$STRING_id),
                              ]$gene)

  # Remove edges that are redundant, the output graph is not directed
  net=simplify(net, remove.multiple=TRUE)

  if(add.missing){
    indx.miss=which(genes.hgnc %in% V(net)$name)
    if(length(indx.miss)>0){
      missing.nodes=genes.hgnc[-which(genes.hgnc %in% V(net)$name)]
      net=net+vertices(missing.nodes)
    }
  }
  
  # if clusters of net should be color-coded
  if(cluster.nodes){

    ceb <- cluster_edge_betweenness(net)
    if(plot_dendrogram)  dendPlot(ceb, mode="hclust")

    ceb <- cluster_edge_betweenness(net)
    if(!is.null(cut_at)){
      cut <- cutat(ceb,cut_at)
      colors <- rainbow(cut_at)

      # color-coded halo effect around nodes
      df=data.frame(index=1:length(V(net)$name),
                    name=V(net)$name,
                    group=cut)
      df.split <- split(df, df$group)
      df.split = lapply(df.split, function(x) x$index)
      return(
          plot(net,
               mark.groups=df.split,
               mark.col=rainbow(length(df.split), alpha=0.5),
               mark.border=NA,
               ...)
      )
    }
  }

    # color-code network by a given data.frame
  if(is.data.frame(nodes.color)){
    V(net)$color <- nodes.color$color[match(V(net)$name, nodes.color$gene)]
  }

  if(is.data.frame(mark.color)){

    mark.color2=data.frame(gene=genes.hgnc,
                           color=mark.color[match(genes.hgnc, mark.color$gene),]$color)
    mark.color2=mark.color2[-is.na(mark.color2$color),]
    mark.color2$id=1:nrow(mark.color2)
    my.mark.groups = lapply(split(mark.color2, mark.color2$color),
                            function(x) x$id)

    return(
      plot.igraph(net,
                  mark.groups=my.mark.groups,
                  mark.col=names(my.mark.groups),
                  ...)
    )
  }

  return(
    plot.igraph(net,
                ...)
  )

}



##################################################################################
############ Heatmaps
##################################################################################


#' # Plot a heatmap of motifs activities
#'
#' @param values a numeric matrix
#' @param top.annotation a list of vectors to plot as a top annotation
#' @param row.annotation a list of vectors to plot as a row annotation
#' @param col.top.annotation a list of colors that correspond to `top.annotation`
#' @param col.row.annotation a list of colors that correspond to `row.annotation`
#'
#' @return
#' @docType methods
#' @export
plotHeatmap = function(values,
                       col,
                       top.annotation=NULL,
                       row.annotation=NULL,
                       col.top.annotation=NULL,
                       col.row.annotation=NULL,
                       title="Heatmap",
                       cluster_columns =TRUE,
                       match.heat.dend=FALSE){
  #require(ComplexHeatmap)

  .get.heatmap = function(values,
                          col,
                         top.annotation,
                         row.annotation,
                         col.top.annotation,
                         col.row.annotation,
                         title="Heatmap",
                         cluster_columns,
                         column_dend_reorder=TRUE){

    if(!is.null(top.annotation)){
      top.annotation = HeatmapAnnotation(
        df=as.data.frame(top.annotation),
        col=col.top.annotation,
        show_annotation_name = TRUE,
        na_col = "grey"
      )
    }

    if(!is.null(row.annotation)){
      row.annotation = rowAnnotation(
        df=as.data.frame(row.annotation),
        col=col.row.annotation,
        show_annotation_name = TRUE,
        na_col = "grey")
    }


        Heatmap(values,
                col=col,
                name = title,
                show_column_names = TRUE,
                cluster_columns = cluster_columns,
                right_annotation = row.annotation,
                top_annotation = top.annotation,
                row_names_gp = gpar(fontsize = 8),
                column_dend_reorder=column_dend_reorder)
  }


  if(!is.list(values)){

    ht_list1 = .get.heatmap(values,
                          col=col,
                         top.annotation,
                         row.annotation,
                         col.top.annotation,
                         col.row.annotation,
                         title="Heatmap",
                         cluster_columns)


  }else{

    my.column.order=TRUE
    ht_list = lapply(1:length(values), function(i){

     if(match.heat.dend){
       if(i==1){
         heat=.get.heatmap(values[[i]],
                           col=col[[i]],
                           top.annotation[[i]],
                           row.annotation[[i]],
                           col.top.annotation[[i]],
                           col.row.annotation[[i]],
                           title=title[[i]],
                           cluster_columns=cluster_columns)
         my.column.order = column_order(heat)
       }else{
         heat=.get.heatmap(values[[i]],
                           col=col[[i]],
                           top.annotation[[i]],
                           row.annotation[[i]],
                           col.top.annotation[[i]],
                           col.row.annotation[[i]],
                           title=title[[i]],
                           cluster_columns=cluster_columns,
                           column_dend_reorder=my.column.order)
       }

      # no matching columns of input matrices
     }else{
       .get.heatmap(values[[i]],
                    col=col[[i]],
                    top.annotation[[i]],
                    row.annotation[[i]],
                    col.top.annotation[[i]],
                    col.row.annotation[[i]],
                    title=title[[i]],
                    cluster_columns=cluster_columns,
                    column_dend_reorder=my.column.order)
     }
     })

   ht_list1 = NULL  ## Heatmap(...) + NULL gives you a HeatmapList object
   for(s in ht_list) {
     ht_list1 = ht_list1 + s
   }

  }

  ht_list1
 #draw(ht_list1, ht_gap = unit(3, "mm"))

}












