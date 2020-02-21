# these fx need proper citation. They are collected from here and there and customized as needed.


genome_loadings <- function (component = NULL, max.items = 20, label.size = 3, label.repel = 1, 
          label_both = TRUE, label_X = FALSE, min_loading = 0.01, gene_locations = NULL, 
          chromosome_lengths = NULL, hide_unknown = FALSE, highlight_genes = NULL, 
          label_genes = NULL, GeneLoadings = GeneLoadings) 
{
  # component = GeneLoadings[1,]
  
  # label_both = T; 
  # max.items = 10; 
  # gene_locations =   gene_locations;
  # chromosome_lengths = chromosome.lengths; 
  # GeneLoadings = GeneLoadings
  
  
  temp <- data.table(loading = component, gene_symbol = names(component))
  
  if (is.null(gene_locations)) {
    gene_locations <- load_gene_locations(colnames(GeneLoadings))
  }
  if (is.null(chromosome_lengths)) {
    message("Using mmusculus_gene_ensembl chromosome lengths")
    chromosome_lengths <- load_chromosome_lengths(organism = "mmusculus_gene_ensembl")
  }
  setkey(temp, gene_symbol)
  setkey(gene_locations, gene_symbol)
  temp <- merge(temp, gene_locations, all.x = TRUE)
  temp$chromosome <- factor(temp$chromosome_name, levels = c(1:22, 
                                                             "X", "Y", "MT", "?"))
  temp[is.na(chromosome)]$chromosome <- "?"
  temp[chromosome == "?"]$start_position <- sample(1:chromosome_lengths[chromosome == 
                                                                          "?"]$length, nrow(temp[chromosome == "?"]))
  setkey(temp, chromosome)
  temp <- chromosome_lengths[temp]
  temp[, `:=`(genomic_position, genomic_offset + start_position)]
  if (label_both == TRUE) {
    label_data <- temp[abs(loading) > min_loading][order(-abs(loading))][1:max.items]
  }
  else {
    which_size_max <- as.logical(temp[abs(loading) == max(abs(loading)), 
                                      sign(loading)] - 1)
    label_data <- temp[abs(loading) > min_loading][order(-loading, 
                                                         decreasing = which_size_max)][1:max.items]
  }
  if (!is.null(label_genes)) {
    label_data <- unique(rbind(label_data, temp[gene_symbol %in% 
                                                  label_genes]))
  }
  if (label_X == TRUE) {
    label_data <- rbind(label_data, temp[abs(loading) > min_loading][chromosome_name == 
                                                                       "X"])
  }
  levels(temp$chromosome)[levels(temp$chromosome) == "MT"] <- "M"
  cl <- chromosome_lengths
  levels(cl$chromosome)[levels(cl$chromosome) == "MT"] <- "M"
  if (hide_unknown) {
    temp <- temp[chromosome != "?"]
    cl <- cl[chromosome != "?"]
  }
  P <- ggplot(temp, aes(genomic_position, loading, size = abs(loading)^2)) + 
    geom_point(stroke = 0, aes(alpha = (abs(loading))^0.7, 
                               color = chromosome)) + scale_colour_manual(values = c(rep_len(c("black", 
                                                                                               "cornflowerblue"), length(levels(temp$chromosome))), 
                                                                                     "grey")) + xlab("Genomic Coordinate") + ylab("Loading") + 
    theme_minimal() + theme(legend.position = "none") + scale_x_continuous(breaks = cl$center, 
                                                                           labels = cl$chromosome, minor_breaks = NULL) + geom_label_repel(data = label_data, 
                                                                                                                                           aes(label = gene_symbol), size = label.size, box.padding = unit(0.5, 
                                                                                                                                                                                                           "lines"), point.padding = unit(0.1, "lines"), force = label.repel, 
                                                                                                                                           segment.size = 0.2, segment.color = "blue")
  if (!is.null(highlight_genes)) {
    P <- P + geom_point(data = temp[gene_symbol %in% highlight_genes], 
                        color = "red")
  }
  return(P)
}



# print_loadings_scores <- function(SDAScores = CellScores, SDALoadings = GeneLoadings, ComponentN=NA, ColFac = NA, Prefix="SDAV", GeneLoc=NA,
#                                   chromosome.lengths = chromosome.lengths){
#   library(ggthemes)
#   library(scales)
#   if(!is.factor(ColFac)) ColFac <- factor(ColFac)
#   
#   # SDAScores    <- SDAResult$scores
#   # SDALoadings <- SDAResult$loadings[[1]]
#   
#   
#   #cowplot::plot_grid(
#   g1 <- ggplot(data.table(cell_index = 1:nrow(SDAScores), 
#                           score = SDAScores[, paste0(Prefix, ComponentN)], experiment = gsub("_.*", 
#                                                                                              "", gsub("[A-Z]+\\.", "", rownames(SDAScores))), ColFac = ColFac), 
#                aes(cell_index, score, colour = ColFac)) + 
#     geom_point(size = 0.5, stroke = 0) + 
#     xlab("Cell Index") + ylab("Score") + 
#     #scale_color_brewer(palette = "Paired") + 
#     
#     
#     theme_bw() + 
#     theme(legend.position = "bottom") + 
#     guides(colour = guide_legend(override.aes = list(size=2, alpha=1))) +
#     scale_colour_manual(values = colorRampPalette(solarized_pal()(8))(length(levels(ColFac))),
#                         guide = guide_legend(nrow=2)) +
#     #guides(color = guide_legend(ncol = 2, override.aes = list(size = 2))) + 
#     ggtitle(paste0(Prefix, ComponentN)) 
#   #,
#   
#   # g2 <- genome_loadings(SDALoadings[ComponentN,], 
#   #                       label_both = T, 
#   #                       max.items = 10, 
#   #                       gene_locations =   GeneLoc,
#   #                       chromosome_lengths = chromosome.lengths)
#   # #, ncol = 1)
#   print(g1)
#   # print(g2)
# }


simplify <- theme(legend.position = "none",
                  axis.title.x=element_blank(),
                  axis.title.y=element_blank(),
                  axis.text.x = element_blank(),
                  axis.text.y = element_blank(),
                  axis.ticks = element_blank())

simplify2 <- theme(legend.position = "right",
                  axis.title.x=element_blank(),
                  axis.title.y=element_blank(),
                  axis.text.x = element_blank(),
                  axis.text.y = element_blank(),
                  axis.ticks = element_blank())

print_gene_list <- function(i, PrintRes = F, PosOnly = F, NegOnly = F, AbsLoad = T, TopN = 150, GeneLoadings = GeneLoadings) {
  # colnames(GeneLoadings)
  # GeneLoadings <- t(GeneLoadings)
  
  if(AbsLoad)  tmp <- data.table(as.matrix(GeneLoadings[i,]), keep.rownames = TRUE)[order(-abs(V1))][1:TopN]
  
  if(PosOnly)  tmp <- data.table(as.matrix(GeneLoadings[i,]), keep.rownames = TRUE)[order(-(V1))][1:TopN]
  
  if(NegOnly)  tmp <- data.table(as.matrix(GeneLoadings[i,]), keep.rownames = TRUE)[order((V1))][1:TopN]
  
  
  setnames(tmp, c("Gene.Name","Loading"))
  setkey(tmp, Gene.Name)
  
  # Add is testis enriched annotations
  #tmp <- merge(tmp, fantom_summary_subset, all.x = TRUE)
  #tmp$Is_Testis_Enriched <- gsub("FALSE","",tmp$fold_difference > 10)
  
  # Add infertility gene annotations
  #tmp$Is_Infertility_Gene <- gsub("FALSE","",tmp$Gene.Name %in% infertility_genes)
  #setcolorder(tmp, c("Gene.Name", "Loading", "Is_Infertility_Gene", "Is_Testis_Enriched", #"fold_difference"))
  
  # Display Result
  if(PrintRes) print(tmp[order(-abs(Loading))]) else return(tmp[order(-abs(Loading))])
}


go_volcano_plot <- function(x=GO_data, component="V5N", extraTitle=""){
  if(extraTitle=="") extraTitle = paste("Component : ", component, sep="")
  
  #print(
  ggplot(data.table(x[[component]]), aes(GeneOdds/BgOdds, -log(pvalue), size=Count)) +
    geom_point(aes(colour=p.adjust<0.05)) +
    scale_size_area() +
    geom_label_repel(data = data.table(x[[component]])[order(p.adjust)][1:30][p.adjust<0.7], aes(label = Description, size=0.25), size = 3, force=2) + 
    ggtitle(paste("",extraTitle, sep="\n") ) +
    xlab("Odds Ratio") +
    scale_x_log10(limits=c(1,NA), breaks=c(1,2,3,4,5,6,7,8))
  #)
}



plotEnrich <- function(GeneSetsDF, GeneVec, plotTitle="", xLab="", N = NULL, k = NULL, ReturnPval=T, BiPlot = F){
  
  
  GeneVec_overlap <- apply(GeneSetsDF, 2, function(x){
    length(which(x %in% GeneVec))
  })
  table(GeneVec_overlap)
  
  
  
  m = length(GeneVec)
  n = N - m
  
  ## Random expectation
  marked.proportion <- m / N; marked.proportion
  exp.x <- k * marked.proportion; exp.x
  
  x = GeneVec_overlap
  
  ## Fold enrichment, as computed by David
  fold.enrichment <-  (x / k ) / (m / N)
  
  # barplot(fold.enrichment, las=2)
  
  
  
  # ggplot(data=data.frame(x=names(fold.enrichment),
  #            y=fold.enrichment), aes(x=x, y=y)) + theme_bw()  +
  #   geom_bar(stat="identity") +
  #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Fold-enrichment \n Sertoli cell-only syndrome vs control DE genes")
  
  p.value <-  phyper(q=x-1, m=m, n=n, k=k, lower.tail=FALSE)
  # p.value <-  p.adjust(p.value, "BH")
  p.value <-  p.adjust(p.value, "fdr")
  
  fold.enrichment <- fold.enrichment[gtools::mixedsort(names(fold.enrichment))]
  
  p.value <-  p.value[names(fold.enrichment)]
  
  tDF <- data.frame(x=(fold.enrichment),
             y=p.value)
  rownames(tDF) <- names(fold.enrichment)
  # tDF[tDF$y <0.01,]
  
  # print(head(p.value))
  # print(head(names(fold.enrichment)))
  # print(names(p.value)[p.value<0.01])
  # print(paste(names(p.value)[p.value<0.01]))
  # print(rownames(tDF[tDF$y <0.01,]))
  
  tempSigComps <- rownames(tDF[tDF$y <0.01,])
  
  tempLen <- length(tempSigComps)
  
  tempKeepLen <- tempLen - length(grep("Removed", tempSigComps))
  
  if(tempLen> 16){
    plotTitle <- paste0(plotTitle, "\nSig Comps ", tempKeepLen, "/", tempLen, ":", 
                        paste(rownames(tDF[tDF$y <0.01,])[1:7], collapse = ", "), "\n", 
                        paste(rownames(tDF[tDF$y <0.01,])[8:15], collapse = ", "), "\n", 
                        paste(rownames(tDF[tDF$y <0.01,])[16:tempLen], collapse = ", "))
    
  } else  if(tempLen> 8){
    plotTitle <- paste0(plotTitle, "\nSig Comps: ", tempKeepLen, "/", tempLen, ":",
                        paste(rownames(tDF[tDF$y <0.01,])[1:7], collapse = ", "), "\n", 
                        paste(rownames(tDF[tDF$y <0.01,])[8:tempLen], collapse = ", "))
    
  } else {
    plotTitle <- paste0(plotTitle, "\nSig Comps: ", tempKeepLen, "/", tempLen, ":",
                        paste(rownames(tDF[tDF$y <0.01,]), collapse = ", "))
    
  }
  
  if(BiPlot) print(ggplot(data=tDF, aes(x=x, y=y, label = ifelse(y < 0.01, "*", ""))) + theme_bw()  +
          geom_point() + geom_line() +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
          geom_text(vjust = 0) +
          ggtitle(plotTitle) +
          xlab("Fold-enrichment") +
          ylab("1 - P-value = 1 - P(X>=x) = P(X<x)"))
  
  tDF <- data.frame(x=factor(names(fold.enrichment), levels = names(fold.enrichment)),
                       y=fold.enrichment,
                       p=p.value)
  
  rownames(tDF) <- names(fold.enrichment)
  
  

  print(ggplot(data=tDF, aes(x=x, y=y, label = ifelse(p < 0.01, "*", ""))) + theme_bw()  +
          geom_bar(stat="identity") +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
          geom_text(vjust = 0) + 
          ggtitle(plotTitle) + 
          xlab(xLab) + 
          ylab("Fold-enrichment"))
  
  
  if(ReturnPval) return(1-p.value)
  
}
