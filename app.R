
# ShinySDAscSeq (ShinSASe): A Shiny template for scRNASeq analysis processed with SDA

library(shiny)
library(rclipboard)
library(shinydashboard)

library(shiny)
library(ggplot2)
library(data.table)
library(ggrepel)
library(viridis)
library(ggnewscale)
library(RColorBrewer)
library(grid)
library(gridExtra) 
library(dplyr)

source("./Fxs.R")


list2env(readRDS( "./data/ShinyServerDataLS.rds"), envir = globalenv())

#ShinyServerDataLS is a list of data:

#"datat" a dataframe with cell barcode + meta data including 1 dimreduc (x,y) like tsne. The other paramters differ per study so custom.
#"results" : an SDA object post processing and label correction holds the cell scores and gene loading
#"chromosome.lengths" : a dataframe with "chromosome", "length", "length_padded", "genomic_offset", "center" from mapping genes to biomart
#"StatFac": A dataframe of "SDA" comp and annotations found in processing, labels to be shown when a comp is selected
#"GO_data": a list of GO annotations, 1 per every component in the resuts
#"gene_locations": from biomart
#"col_vector" : vector of colors to be used



ui <- dashboardPage(skin="red",
  dashboardHeader(title = "EM_KO_BB_EXP1",
                  dropdownMenu(type = "messages",
                               messageItem(
                                 from = "A",
                                 message = "Task 1 due."
                               ),
                               messageItem(
                                 from = "B",
                                 message = "question?",
                                 icon = icon("question"),
                                 time = "13:45"
                               ),
                               messageItem(
                                 from = "C",
                                 message = "Woohoo!",
                                 icon = icon("life-ring"),
                                 time = "2014-12-01"
                               )
                  ),
                  dropdownMenu(type = "notifications",
                               notificationItem(
                                 text = "ABC",
                                 icon("users")
                               ),
                               notificationItem(
                                 text = "XYZ",
                                 icon("rep"),
                                 status = "success"
                               ),
                               notificationItem(
                                 text = "IJK",
                                 icon = icon("exclamation-triangle"),
                                 status = "warning"
                               )
                  ),
                  #red, yellow, aqua, blue, light-blue, green, navy, teal, olive, lime, orange, fuchsia, purple, maroon, black.
                  dropdownMenu(type = "tasks", badgeStatus = "success",
                               taskItem(value = 99, color = "green",
                                        "Preprocessing"
                               ),
                               taskItem(value = 99, color = "aqua",
                                        "SDA processing"
                               ),
                               taskItem(value = 99, color = "yellow",
                                        "DE analysis"
                               ),
                               taskItem(value = 70, color = "olive",
                                        "Documentation"
                               ),
                               taskItem(value = 15, color = "black",
                                        "Manuscript"
                               ),
                               taskItem(value = 70, color = "red",
                                        "Overall project"
                               )
                  )
                  
  ),
  
  
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Main Tab", tabName = "combodash", icon = icon("dashboard"),
               badgeLabel = "underconst.", badgeColor = "yellow"),
      menuItem("Tab 2. (dev)", tabName = "germdash", icon = icon("affiliatetheme"),
               badgeLabel = "soon", badgeColor = "red"),
      menuItem("Tab 3. (dev)", tabName = "somadash", icon = icon("allergies"),
               badgeLabel = "soon", badgeColor = "red"),
      menuItem("Tab 3. (dev)", tabName = "pseudotime", icon = icon("arrows-alt"),
               badgeLabel = "soon", badgeColor = "red"),
      menuItem("Enrichment Analysis", tabName = "Enrichment", icon = icon("dashboard"),
               badgeLabel = "underconst.", badgeColor = "yellow"),
      menuItem("Source code", icon = icon("file-code-o"), 
               href = "https://")
    )
  ),
  
  dashboardBody(
    tags$head(
      tags$style(HTML("
                      .content-wrapper {
                      background-color: black !important;
                      }
                      .main-sidebar {
                      background-color: black !important;
                      }
                      "))),
    tabItems(
      # First tab content
      tabItem(tabName = "combodash",
              fluidRow(
                
                box(
                  title = "Inputs", status = "warning", solidHeader = TRUE,
                  #"Box content here", br(), "More box content",
                  #sliderInput("ComponentN", "Slider input:", 1, 40, 1),
                  
                  textInput("ComponentNtext", "Numerical input:", "1"),
                  textInput("Genetext", "Text input:", "CD69"),
                  radioButtons("data", "Data origin:",  # ADD or REMOVE depending on what is in the data
                               c("tSNE - Batch-removed SDA-CellScores" = "tsneBrSDACS",
                                 "tSNE - Original SDA-CellScores" = "tsneDSNDGE",
                                 "tSNE - Batch-removed-Imputed-SER-DGE" = "tsneImpDGE"
                                 )),
                  radioButtons("metaselect", "Metadata Selection:",
                               c("SubjID" = "SubjectId",
                                 "SampleDate" = "SampleDate",
                                 "Cluster" = "CustClus",
                                 "Cluster2" = "CustClus2",
                                 "BarcodePrefix" = "BarcodePrefix",
                                 "SingleR_Labels" = "SingleR_Labels"
                                 
                               )),
                  textInput("NoOfGenes", "No. of Genes to output:", "20"),
                  actionButton("C2Cpos", "Copy2ClipPosGenes"),
                  actionButton("C2Cneg", "Copy2ClipNegGenes"),
                  actionButton("PrevSDA", "Prev SDA"),
                  actionButton("NextSDA", "Next SDA"),
                  width = 3
                ),
               
                
                valueBoxOutput("cellinfo1", width = 3),
                
                valueBoxOutput("GeneName", width = 3)
                
                
              ), 
              
              box(
                #actionButton("button", "Next"),
                title = "Cell score", status = "primary", solidHeader = TRUE,
                collapsible = TRUE,
                plotOutput("plot1"), #plotlyOutput
                width = 5, background = "black"
              ),
              box(
                title = "Gene expr", status = "primary", solidHeader = TRUE,
                collapsible = TRUE,
                plotOutput("plot4"), #plotlyOutput
                width = 5, background = "black"
              ),
              box(
                title = "Pheno - tSNE", status = "primary", solidHeader = TRUE,
                collapsible = TRUE,
                plotOutput("plot3"), #plotlyOutput
                width = 5, background = "black"
              ),
              box(
                title = "Pheno Legend", status = "primary", solidHeader = TRUE,
                collapsible = TRUE,
                plotOutput("plot2"),
                width = 5, background = "black"
              ),
              
              box(
                title = "Pos. Loadings GO", status = "primary", solidHeader = TRUE,
                collapsible = TRUE,
                plotOutput("GOpos"), #plotlyOutput
                width = 5, background = "black"
              ),
              box(
                title = "Neg. Loadings GO", status = "primary", solidHeader = TRUE,
                collapsible = TRUE,
                plotOutput("GOneg"),
                width = 5, background = "black"
              ),
              
              box(
                title = "Cell Scores Across", status = "primary", solidHeader = TRUE,
                collapsible = TRUE,
                plotOutput("SDAScoresAcross"),
                width = 10, background = "black"
              ),
              
              box(
                title = "Chrom. Location", status = "primary", solidHeader = TRUE,
                collapsible = TRUE,
                plotOutput("ChrLoc"),
                width = 10, background = "black"
              ),
              
              box(title = "Pos. Top Genes", status = "info", solidHeader = TRUE, width = 4, background = "black",
                  tableOutput("packageTablePos")
              ),
              box(title = "Neg. Top Genes", status = "info", solidHeader = TRUE, width = 4, background = "black",
                  tableOutput("packageTableNeg")
              )
              
              
              
      ),
      
      # Germ tab content
      tabItem(tabName = "germdash",
              h2("Germ only tab content")
      ),
      
      # Soma tab content
      tabItem(tabName = "somadash",
              h2("Soma only tab content")
      ),
      
      # pseudotime tab content
      tabItem(tabName = "pseudotime",
              h2("Pseudotime tab content")
      ),
      
      # Enrichment tab content
      tabItem(tabName = "Enrichment",
              fluidRow(
                box(
                  title = "Inputs", status = "warning", solidHeader = TRUE,
                  "Multiple formatting of gene sets accepted", 
                  br(), "List can be seperated by comma e.g. from ", 
                  br(), "   or spaces e.g. from Excel", 
                  br(), "Also, single or double quotes or not",
                  #sliderInput("ComponentN", "Slider input:", 1, 40, 1),
                  textInput("GeneSet", "A set of genes", "'GZMA', 'GZMK', 'MKI67'"),
                  width = 8
                ),
                
                box(
                  title = "Positive Loadings", status = "primary", solidHeader = TRUE,
                  collapsible = TRUE,
                  plotOutput("plot5"),
                  width = 10
                ),
                
                box(
                  title = "Negative Loadings", status = "primary", solidHeader = TRUE,
                  collapsible = TRUE,
                  plotOutput("plot6"),
                  width = 10
                )
              )
      )
      
    )
  )
  
  
)


server <- function(input, output, session) {
  
  
  datat <- as.data.frame(cbind(datat, results$scores[rownames(datat),])); rownames(datat) <- datat$barcode 
  
  datat <- data.table(datat)
  
  
  
  ScoreAcrossDT <- reactive({
    #these depend on metadata 
    
    if(input$metaselect == "SubjectId") {
      ColFac_DONR.ID <- sort(as.character(datat$SubjectId))
      names(ColFac_DONR.ID) <- rownames(datat)[order(as.character(datat$SubjectId))]
      
      ColFac_DONR.ID
    } else {
      if(input$metaselect == "SampleDate"){
        ColFac_DONR.ID <- sort(as.character(datat$SampleDate))
        names(ColFac_DONR.ID) <- rownames(datat)[order(datat$SampleDate)]
        
      } else {
        if(input$metaselect == "CustClus"){
          ColFac_DONR.ID <- sort(as.character(datat$CustClus))
          names(ColFac_DONR.ID) <- rownames(datat)[order(datat$CustClus)]
          
        } else {
          if(input$metaselect == "CustClus2"){
            ColFac_DONR.ID <- sort(as.character(datat$CustClus2))
            names(ColFac_DONR.ID) <- rownames(datat)[order(datat$CustClus2)]
            
          } else {
            if(input$metaselect == "BarcodePrefix"){
              ColFac_DONR.ID <- sort(as.character(datat$BarcodePrefix))
              names(ColFac_DONR.ID) <- rownames(datat)[order(datat$BarcodePrefix)]
              
            } else {
              if(input$metaselect == "SingleR_Labels"){
                ColFac_DONR.ID <- sort(as.character(datat$SingleR_Labels))
                names(ColFac_DONR.ID) <- rownames(datat)[order(datat$SingleR_Labels)]
                
              } else {
                
              }
            }
          }
        }
      }
    }
    

    
    SDAScores <- results$scores
    ComponentN <- as.numeric(input$ComponentNtext)
    
    DTt <- data.frame(cell_index = 1:nrow(SDAScores), 
                      score = SDAScores[, paste0("SDAV", ComponentN)], 
                      experiment = gsub("_.*", "", gsub("[A-Z]+\\.", "", rownames(SDAScores))),
                      ColFac = as.character(ColFac_DONR.ID),
                      row.names = names(ColFac_DONR.ID))
    # DTt$ColFac <- "unk"
    
    # DTt[names(ColFac_DONR.ID),]$ColFac <- as.character(ColFac_DONR.ID)
    DTt
  })
  
  
  

  tDF <- reactive({

    
    if(input$data == "tsneBrSDACS") {
      tempDF <- as.data.frame(datat)[, c("Tsne1_SDAQC2", "Tsne2_SDAQC2")]; colnames(tempDF) <- c("tSNE1", "tSNE2")
      tempDF
    } else {
        if(input$data == "tsneDSNDGE") {
          tempDF <- as.data.frame(datat)[, c("Tsne1_SDAQC", "Tsne2_SDAQC")]; colnames(tempDF) <- c("tSNE1", "tSNE2")
         } else {
            if(input$data == "tsneImpDGE"){
              tempDF <- as.data.frame(datat)[, c("tSNE_SerImp_1", "tSNE_SerImp_2")]; colnames(tempDF) <- c("tSNE1", "tSNE2")
              tempDF$tSNE1 <- as.numeric(tempDF$tSNE1)
              tempDF$tSNE2 <- as.numeric(tempDF$tSNE2)
              
         }
        
         }
      tempDF
    }
    
    rownames(tempDF) <- datat$barcode
    

    tempDF$GeneExpr <- rep(0, nrow(tempDF))

    return(tempDF)
  })
  

    observeEvent(input$NextSDA, {
    Val <- as.character(min(c(40, as.numeric(input$ComponentNtext)+1)))
    
    updateTextInput(session, "ComponentNtext", value = Val)
  })
  
  observeEvent(input$PrevSDA, {
    Val <- as.character(max(c(1, as.numeric(input$ComponentNtext)-1)))
    
    updateTextInput(session, "ComponentNtext", value = Val)
  })
  
  observeEvent(input$C2Cpos, {
 
    
    Out1 <- print_gene_list(as.numeric(input$ComponentNtext), PosOnly = T) %>%
      as.data.frame() %>%
      head(as.numeric(input$NoOfGenes)) 
    Out1 <- Out1$Gene.Name
    
    clipr::write_clip(Out1)
    
  })
  observeEvent(input$C2Cneg, {
  
    
    Out2 <- print_gene_list(as.numeric(input$ComponentNtext), NegOnly = T) %>%

      as.data.frame() %>%
      head(as.numeric(input$NoOfGenes)) 
    Out2 <- Out2$Gene.Name
    
    clipr::write_clip(Out2)
    
  })
  
  

  
  output$packageTablePos <- renderTable({
    print_gene_list(as.numeric(input$ComponentNtext), PosOnly = T) %>%
      as.data.frame() %>%
      head(as.numeric(input$NoOfGenes))
  }, digits = 1)
  
  output$packageTableNeg <- renderTable({
    print_gene_list(as.numeric(input$ComponentNtext), NegOnly = T) %>%
      as.data.frame() %>%
      head(as.numeric(input$NoOfGenes))
  }, digits = 1)
  
  
  
  output$cellinfo1 <- renderValueBox({
     valueBox(
      value = StatFac[paste0("SDAV", input$ComponentNtext, sep=""),2], #format(Sys.time(), "%a %b %d %X %Y %Z"),
      subtitle = StatFac[paste0("SDAV", input$ComponentNtext, sep=""),6],
      icon = icon("area-chart"),
      color = "yellow" #if (downloadRate >= input$rateThreshold) "yellow" else "aqua"
    )
  })
  
  output$GeneName <- renderValueBox({
    if(input$Genetext %in% colnames(results$loadings[[1]])){
      # results$loadings[[1]][,"PRM1"]
      GeneName <- input$Genetext
    } else {
      GeneName <- paste0(input$Genetext, " Not Found")
      
    }
    
    valueBox(
      value = GeneName, 
      subtitle = "Gene Name",
      icon = icon("area-chart"),
      color = "yellow" 
    )
  })
  
  
  #renderPlotly
  output$plot1 <- renderPlot({
    #ggplotly
    tempDF <- tDF()
    (ggplot(cbind(tempDF, SDAComp=datat[,get(paste0("SDAV", input$ComponentNtext, sep=""))]), 
            aes(tSNE1, tSNE2, color=cut(asinh(SDAComp), breaks = c(-Inf, -1, -.5, 0, .5, 1, Inf)))) +
       geom_point(size=0.1) + theme_bw() +
       scale_color_manual("CS", values = rev(c("red", "orange", "yellow", "lightblue", "dodgerblue", "blue")) ) + 
       guides(colour = guide_legend(override.aes = list(size=2, alpha=1))) +
       theme(legend.position = "bottom", aspect.ratio=1,
             panel.background = element_rect(fill = "black",
                                             colour = "black",
                                             size = 0.5, linetype = "solid")) + 
       ggtitle(paste0("SDAV", input$ComponentNtext, " \n", 
                      StatFac[paste0("SDAV", input$ComponentNtext, sep=""),2], sep="")) + 
       simplify2 + coord_cartesian(xlim = NULL, ylim = NULL, expand = FALSE))
    
  })
  
  output$plot2 <- renderPlot({
    tempDF <- tDF()
    
    if(input$metaselect == "SubjectId") {
      MetaFac <- (datat$SubjectId)
    } else {
        if(input$metaselect == "SampleDate"){
          MetaFac <- (datat$SampleDate)
        } else {
          if(input$metaselect == "CustClus"){
            MetaFac <- (datat$CustClus)
          } else {
            if(input$metaselect == "CustClus2"){
              MetaFac <- (datat$CustClus2)
            } else {
              if(input$metaselect == "BarcodePrefix"){
                MetaFac <- (datat$BarcodePrefix)
              } else {
                if(input$metaselect == "SingleR_Labels"){
                  MetaFac <- (datat$SingleR_Labels)
                } else {
                  
                }
              }
            }
          }
        }
    }


    
    ggFph <- ggplot(tempDF, aes(tSNE1, tSNE2, color=factor(MetaFac))) +
      geom_point(size=0.1)+ theme_bw() +
      theme(legend.position = "bottom", aspect.ratio=1,
            legend.title = element_blank()) +
      ggtitle("t-SNE - Final Pheno") +
      scale_color_manual(values=rev(col_vector)) + guides(colour = guide_legend(override.aes = list(size=2, alpha=1), ncol =2))  +
      coord_cartesian(xlim = NULL, ylim = NULL, expand = FALSE)
    legend <- cowplot::get_legend(ggFph)
    
    #grid.newpage()
    grid.draw(legend)
    
  })
  #renderPlotly
  output$plot3 <- renderPlot({
    tempDF <- tDF()
    
    if(input$metaselect == "SubjectId") {
      MetaFac <- (datat$SubjectId)
    } else {
      if(input$metaselect == "SampleDate"){
        MetaFac <- (datat$SampleDate)
      } else {
        if(input$metaselect == "CustClus"){
          MetaFac <- (datat$CustClus)
        } else {
          if(input$metaselect == "CustClus2"){
            MetaFac <- (datat$CustClus2)
          } else {
            if(input$metaselect == "BarcodePrefix"){
              MetaFac <- (datat$BarcodePrefix)
            } else {
              if(input$metaselect == "SingleR_Labels"){
                MetaFac <- (datat$SingleR_Labels)
              } else {
                
              }
            }
          }
        }
      }
    }
    
    #ggplotly
    (ggplot(tempDF, aes(tSNE1, tSNE2, color=factor(MetaFac))) +
       geom_point(size=0.1)+ theme_bw() +
       theme(legend.position = "none", aspect.ratio=1,
             panel.background = element_rect(fill = "black",
                                             colour = "black",
                                             size = 0.5, linetype = "solid")) +
       ggtitle("t-SNE - Final Pheno") +
       scale_color_manual(values=rev(col_vector))   +
       coord_cartesian(xlim = NULL, ylim = NULL, expand = FALSE))
    
    
  })
  
  output$plot4 <- renderPlot({
    
    tempDF <- as.data.frame(tDF())
    
    if(input$Genetext %in% colnames(results$loadings[[1]])){
      # results$loadings[[1]][,"PRM1"]
      GeneExpr <- results$scores %*% results$loadings[[1]][,as.character(input$Genetext)]
    } else {
      GeneExpr <- results$scores %*% rep(0, nrow(results$loadings[[1]]))

    }
    #ggplotly

    LoadOrdVal <- round(results$loadings[[1]][,as.character(input$Genetext)][order(abs(results$loadings[[1]][,as.character(input$Genetext)]), decreasing = T)], 3)
    

    tempDF[rownames(GeneExpr), ]$GeneExpr <- GeneExpr[,1]
    
    (ggplot(tempDF, 
            aes(tSNE1, tSNE2, color=cut(asinh(GeneExpr), breaks = c(-Inf, -1, -.5, 0, .5, 1, Inf)))) +
       geom_point(size=0.1) + theme_bw() +
       scale_color_manual("EX", values = rev(c("red", "orange", "yellow", "lightblue", "dodgerblue", "blue")) ) + 
       guides(colour = guide_legend(override.aes = list(size=2, alpha=1))) +
       theme(legend.position = "bottom", aspect.ratio=1,
             panel.background = element_rect(fill = "black",
                                             colour = "black",
                                             size = 0.5, linetype = "solid")) + 
       
       simplify2 + 
        coord_cartesian(xlim = NULL, ylim = NULL, expand = FALSE)) + 
      labs(title = paste("Gene: ", input$Genetext, sep=""), 
           subtitle = paste("Found in comps: \n",
                            paste(names(LoadOrdVal)[1:5], collapse = ", "), 
                            "\n",
                            paste(LoadOrdVal[1:5], collapse = ", "), 
                            "\n",
                            paste(names(LoadOrdVal)[6:10], collapse = ", "), 
                            "\n",
                            paste(LoadOrdVal[6:10], collapse = ", "), 
                            "\n"), 
           caption = "Caption here")
    
    

    
    
  })
  
  output$GOpos <- renderPlot({
    
    
    if(! (as.numeric(input$ComponentNtext) %in% 1:40)){
      print("No GO")
    } else {
      go_volcano_plot(component = paste("V", input$ComponentNtext, "P", sep=""))+ theme_bw()+ theme(aspect.ratio = 1)
      
    }
    
  })
  
  
  output$GOneg <- renderPlot({
    
    if(! (as.numeric(input$ComponentNtext) %in% 1:40)){
      print("No GO")
    } else {
      go_volcano_plot(component = paste("V", input$ComponentNtext, "N", sep=""))+ theme_bw()+ theme(aspect.ratio = 1)
      
    }
    
  })
  
  output$ChrLoc <- renderPlot({
    
    if(! (as.numeric(input$ComponentNtext) %in% 1:40)){
      print("No Comp")
    } else {
      pgl <- genome_loadings(results$loadings[[1]][as.numeric(input$ComponentNtext),], 
                             label_both = T, 
                             max.items = 10, 
                             gene_locations =   gene_locations,
                             chromosome_lengths = chromosome.lengths)+ theme(aspect.ratio = .5)
      print(pgl)
      
    }
    
  })
  
  
  
  output$SDAScoresAcross <- renderPlot({
    
    DTt <- ScoreAcrossDT()
    ComponentN <- as.numeric(input$ComponentNtext)
    SDAScores <- results$scores
    
    if(! (as.numeric(input$ComponentNtext) %in% 1:40)){
      print("No Comp")
    } else {
      
      
      pgl <- ggplot(data.table(DTt), 
                    aes(cell_index, score, colour = factor(ColFac))) + 
        geom_point(size = 0.7, stroke = 0) + 
        xlab("Cell Index") + ylab("Score") + 
        #scale_color_brewer(palette = "Paired") + 
        theme_bw() + 
        theme(legend.position = "bottom",
              panel.background = element_rect(fill = "black",
                                              colour = "black",
                                              size = 0.5, linetype = "solid")) + 
        guides(colour = guide_legend(override.aes = list(size=2, alpha=1))) +
        scale_colour_manual(values =rev(col_vector),
                            guide = guide_legend(nrow=2)) +
        #guides(color = guide_legend(ncol = 2, override.aes = list(size = 2))) + 
        ggtitle(paste0("SDAV", ComponentN))+
        ylim(quantile(SDAScores[, paste0("SDAV", ComponentN)],.01), quantile(SDAScores[, paste0("SDAV", ComponentN)],.99))
      
      print(pgl)
      
    }
    
  })
  
  output$messageMenu <- renderMenu({
    # Code to generate each of the messageItems here, in a list. This assumes
    # that messageData is a data frame with two columns, 'from' and 'message'.
    msgs <- apply(messageData, 1, function(row) {
      messageItem(from = row[["from"]], message = row[["message"]])
    })
    
    # This is equivalent to calling:
    #   dropdownMenu(type="messages", msgs[[1]], msgs[[2]], ...)
    dropdownMenu(type = "messages", .list = msgs)
  })
  
  output$plot5 <- renderPlot({
    # N = total number of genes (usually not entire genome, since many have unk func)
    N=8025
    # k = number of genes submitted, top 100
    k = 40 #100 #needs stats thinking and adjustment to improve enrichment calls # maybe change to dynamically genes sumbited and in universe/data
    
    GeneSet <- input$GeneSet
    if(length(grep(",", GeneSet)) == 0){
      
      if(length(grep('"', GeneSet)) + length(grep("'", GeneSet))>0) {
        GeneSet <- unlist(strsplit(gsub("'", '', gsub('"', '', GeneSet)), " "))
      } else {
        GeneSet <- unlist(strsplit(GeneSet, " "))
      }
      
      #print(GeneSet)
    }else {
      GeneSet <- (unlist(strsplit(gsub(" ", "", gsub("'", '', gsub('"', '', GeneSet))), ",")))
      #print(GeneSet)
    }
    
    print("length of your genes:")
    print(length(GeneSet))
    GeneSet <- GeneSet[GeneSet %in% colnames(results$loadings[[1]][,])]
    print("length of your genes in this dataset:")
    print(length(GeneSet))
    
    plotEnrich(GeneSetsDF=SDA_Top100pos, 
               GeneVec = GeneSet, 
               plotTitle="Gene-set enrichment\n SDA top 40 pos loadings\n Cust. Input. genes \n Hypergeometric test: * adj.p < 0.01",
               xLab = "SDA Comps",
               N=N,
               k=k)
    
    
  })
  
  output$plot6 <- renderPlot({
    # N = total number of genes (usually not entire genome, since many have unk func)
    N=8025
    # k = number of genes submitted, top 100
    k = 40 #100, correction needed
    GeneSet <- input$GeneSet

    
    if(length(grep(",", GeneSet)) == 0){
      
      if(length(grep('"', GeneSet)) + length(grep("'", GeneSet))>0) {
        GeneSet <- unlist(strsplit(gsub("'", '', gsub('"', '', GeneSet)), " "))
      } else {
        GeneSet <- unlist(strsplit(GeneSet, " "))
      }

      #print(GeneSet)
    }else {
      GeneSet <- (unlist(strsplit(gsub(" ", "", gsub("'", '', gsub('"', '', GeneSet))), ",")))
      #print(GeneSet)
    }

    GeneSet <- GeneSet[GeneSet %in% colnames(results$loadings[[1]][,])]

    
    plotEnrich(GeneSetsDF=SDA_Top100neg, 
               GeneVec = GeneSet, 
               plotTitle="Gene-set enrichment\n SDA top 40 neg loadings\n Cust. Input. genes \n Hypergeometric test: * adj.p < 0.01",
               xLab = "SDA Comps",
               N=N,
               k=k)
    
    
  })
  
  
}


shinyApp(ui, server)

