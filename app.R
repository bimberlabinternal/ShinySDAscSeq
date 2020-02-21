
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

lookup.path <- "../../../Expts/TetCombo2/data"
SerObj.paths <- list.files(lookup.path, full.names = T)[grepl('FinalSerObj', list.files(lookup.path))]
SerObj.2load = 1

SerObj.ShinySDAOut <- readRDS(SerObj.paths[SerObj.2load])

CellScores   <- SerObj.ShinySDAOut@reductions$SDA@cell.embeddings
GeneLoadings <- t(SerObj.ShinySDAOut@reductions$SDA@feature.loadings)

MetaDF      <- SerObj.ShinySDAOut@misc$SDA_processing_results$MetaDF

# names(SerObj.ShinySDAOut@misc$SDA_processing_results)
# names(SerObj.ShinySDAOut@reductions$SDA@misc)



GO_data             <- SerObj.ShinySDAOut@misc$SDA_processing_results$GO_data
chromosome.lengths  <- SerObj.ShinySDAOut@misc$SDA_processing_results$chromosome.lengths
gene_locations      <- SerObj.ShinySDAOut@misc$SDA_processing_results$gene_locations
Remove_comps        <- as.numeric(SerObj.ShinySDAOut@misc$SDA_processing_results$Remove_comps)
Keep_comps          <- setdiff(1:ncol(CellScores), Remove_comps)
SDA_TopNpos             <- SerObj.ShinySDAOut@misc$SDA_processing_results$SDA_TopNpos
SDA_TopNneg             <- SerObj.ShinySDAOut@misc$SDA_processing_results$SDA_TopNneg
TopN                    <- SerObj.ShinySDAOut@misc$SDA_processing_results$TopN



MetaDF$tsne1_CS_raw <- SerObj.ShinySDAOut@misc$SDA_processing_results$tsne_CS_raw$Y[,1]
MetaDF$tsne2_CS_raw <- SerObj.ShinySDAOut@misc$SDA_processing_results$tsne_CS_raw$Y[,2]

MetaDF$tsne1_CS_br <- SerObj.ShinySDAOut@reductions$tSNECSBR@cell.embeddings[,1]
MetaDF$tsne2_CS_br <- SerObj.ShinySDAOut@reductions$tSNECSBR@cell.embeddings[,2]
StatFac <- data.frame(ord=1:ncol(CellScores), Name=colnames(CellScores), row.names = colnames(CellScores))
StatFac$Meta1 <- rep(0, nrow(StatFac))
StatFac$Meta2 <- rep(0, nrow(StatFac))
StatFac$Meta3 <- rep(0, nrow(StatFac))
StatFac$Meta4 <- rep(0, nrow(StatFac))
StatFac$Meta5 <- rep(0, nrow(StatFac))
StatFac$Meta6 <- rep(0, nrow(StatFac))

col_vector <- OOSAP::ColorTheme()$col_vector

# list2env(readRDS( "./data/ShinyServerDataLS.rds"), envir = globalenv())

#ShinyServerDataLS is a list of data:

#"MetaDF" a dataframe with cell barcode + meta data including 1 dimreduc (x,y) like tsne. The other paramters differ per study so custom.
#"results" : an SDA object post processing and label correction holds the cell scores and gene loading
#"chromosome.lengths" : a dataframe with "chromosome", "length", "length_padded", "genomic_offset", "center" from mapping genes to biomart
#"StatFac": A dataframe of "SDA" comp and annotations found in processing, labels to be shown when a comp is selected
#"GO_data": a list of GO annotations, 1 per every component in the resuts
#"gene_locations": from biomart
#"col_vector" : vector of colors to be used



ui <- dashboardPage(skin="red",
  dashboardHeader(title = "ShinSaSe"
                  # ,
                  # dropdownMenu(type = "messages",
                  #              messageItem(
                  #                from = "A",
                  #                message = "Task 1 due."
                  #              ),
                  #              messageItem(
                  #                from = "B",
                  #                message = "question?",
                  #                icon = icon("question"),
                  #                time = "13:45"
                  #              ),
                  #              messageItem(
                  #                from = "C",
                  #                message = "Woohoo!",
                  #                icon = icon("life-ring"),
                  #                time = "2014-12-01"
                  #              )
                  # ),
                  # dropdownMenu(type = "notifications",
                  #              notificationItem(
                  #                text = "ABC",
                  #                icon("users")
                  #              ),
                  #              notificationItem(
                  #                text = "XYZ",
                  #                icon("rep"),
                  #                status = "success"
                  #              ),
                  #              notificationItem(
                  #                text = "IJK",
                  #                icon = icon("exclamation-triangle"),
                  #                status = "warning"
                  #              )
                  # ),
                  # #red, yellow, aqua, blue, light-blue, green, navy, teal, olive, lime, orange, fuchsia, purple, maroon, black.
                  # dropdownMenu(type = "tasks", badgeStatus = "success",
                  #              taskItem(value = 99, color = "green",
                  #                       "Preprocessing"
                  #              ),
                  #              taskItem(value = 99, color = "aqua",
                  #                       "SDA processing"
                  #              ),
                  #              taskItem(value = 99, color = "yellow",
                  #                       "DE analysis"
                  #              ),
                  #              taskItem(value = 70, color = "olive",
                  #                       "Documentation"
                  #              ),
                  #              taskItem(value = 15, color = "black",
                  #                       "Manuscript"
                  #              ),
                  #              taskItem(value = 70, color = "red",
                  #                       "Overall project"
                  #              )
                  # )
                  
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
                               c("tSNE - raw SDA-CellScores" = "tsneCSraw",
                                 "tSNE - batch-removed SDA-CellScores" = "tsneCSbr"
                                 )),
                  radioButtons("metaselect", "Metadata Selection:",
                               c("SubjID" = "SubjectId",
                                 "SampleDate" = "SampleDate",
                                 #"Cluster" = "CustClus",
                                 #"Cluster2" = "CustClus2",
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
                plotOutput("CellScoresAcross"),
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
  
  
  # MetaDF <- as.data.frame(cbind(MetaDF, CellScores[rownames(MetaDF),])); rownames(MetaDF) <- MetaDF$barcode 
  # 
  # MetaDF <- data.table(MetaDF)
  
  
  
  ScoreAcrossDT <- reactive({
    #these depend on metadata 
    
    if(input$metaselect == "SubjectId") {
      ColFac_DONR.ID <- sort(as.character(MetaDF$SubjectId))
      names(ColFac_DONR.ID) <- rownames(MetaDF)[order(as.character(MetaDF$SubjectId))]
      
      # ColFac_DONR.ID
    } else {
      if(input$metaselect == "SampleDate"){
        ColFac_DONR.ID <- sort(as.character(MetaDF$SampleDate))
        names(ColFac_DONR.ID) <- rownames(MetaDF)[order(MetaDF$SampleDate)]
        
      } else {
        if(input$metaselect == "BarcodePrefix"){
          ColFac_DONR.ID <- sort(as.character(MetaDF$BarcodePrefix))
          names(ColFac_DONR.ID) <- rownames(MetaDF)[order(MetaDF$BarcodePrefix)]
          
        } else {
          if(input$metaselect == "SingleR_Labels"){
            ColFac_DONR.ID <- sort(as.character(MetaDF$SingleR_Labels))
            names(ColFac_DONR.ID) <- rownames(MetaDF)[order(MetaDF$SingleR_Labels)]
            
          }  else {
                
              }
            
          
        }
      }
    }
    

    # ComponentN = 1
    ComponentN <- as.numeric(input$ComponentNtext)
    
    DTt <- data.frame(cell_index = 1:nrow(CellScores), 
                      score = CellScores[, paste0("SDA_", ComponentN)], 
                      experiment = gsub("_.*", "", gsub("[A-Z]+\\.", "", rownames(CellScores))),
                      ColFac = as.character(ColFac_DONR.ID),
                      row.names = names(ColFac_DONR.ID))
    # DTt$ColFac <- "unk"
    
    # DTt[names(ColFac_DONR.ID),]$ColFac <- as.character(ColFac_DONR.ID)
    DTt
  })
  
  
  
  
  ### tSNE DF selection

  tDF <- reactive({

    
    if(input$data == "tsneCSraw") {
      tempDF <- as.data.frame(MetaDF)[, c("tsne1_CS_raw", "tsne2_CS_raw")]; colnames(tempDF) <- c("tSNE1", "tSNE2")
      # tempDF
    } else {
        if(input$data == "tsneCSbr") {
          tempDF <- as.data.frame(MetaDF)[, c("tsne1_CS_br", "tsne2_CS_br")]; colnames(tempDF) <- c("tSNE1", "tSNE2")
         } else {
            
        
         }
      # tempDF
    }
    
    rownames(tempDF) <- MetaDF$barcode
    

    tempDF$GeneExpr <- rep(0, nrow(tempDF))

    return(tempDF)
  })
  

    observeEvent(input$NextSDA, {
      
    Val <- as.character(min(c(ncol(CellScores), as.numeric(input$ComponentNtext)+1)))
    
    updateTextInput(session, "ComponentNtext", value = Val)
  })
  
  observeEvent(input$PrevSDA, {
    Val <- as.character(max(c(1, as.numeric(input$ComponentNtext)-1)))
    
    updateTextInput(session, "ComponentNtext", value = Val)
  })
  
  observeEvent(input$C2Cpos, {
 
    
    Out1 <- print_gene_list(as.numeric(input$ComponentNtext), PosOnly = T, GeneLoadings = GeneLoadings) %>%
      as.data.frame() %>%
      head(as.numeric(input$NoOfGenes)) 
    Out1 <- Out1$Gene.Name
    
    clipr::write_clip(Out1)
    
  })
  observeEvent(input$C2Cneg, {
  
    
    Out2 <- print_gene_list(as.numeric(input$ComponentNtext), NegOnly = T, GeneLoadings = GeneLoadings) %>%

      as.data.frame() %>%
      head(as.numeric(input$NoOfGenes)) 
    Out2 <- Out2$Gene.Name
    
    clipr::write_clip(Out2)
    
  })
  
  

  
  output$packageTablePos <- renderTable({
    print_gene_list(as.numeric(input$ComponentNtext), PosOnly = T, GeneLoadings = GeneLoadings) %>%
      as.data.frame() %>%
      head(as.numeric(input$NoOfGenes))
  }, digits = 1)
  
  output$packageTableNeg <- renderTable({
    print_gene_list(as.numeric(input$ComponentNtext), NegOnly = T, GeneLoadings = GeneLoadings) %>%
      as.data.frame() %>%
      head(as.numeric(input$NoOfGenes))
  }, digits = 1)
  
  
  
  output$cellinfo1 <- renderValueBox({
     valueBox(
      value = StatFac[paste0("SDA_", input$ComponentNtext, sep=""),2], #format(Sys.time(), "%a %b %d %X %Y %Z"),
      subtitle = StatFac[paste0("SDA_", input$ComponentNtext, sep=""),6],
      icon = icon("area-chart"),
      color = "yellow" #if (downloadRate >= input$rateThreshold) "yellow" else "aqua"
    )
  })
  
  output$GeneName <- renderValueBox({
    if(input$Genetext %in% colnames(GeneLoadings)){
      # GeneLoadings[,"PRM1"]
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
    (ggplot(cbind(tempDF, SDAComp=CellScores[,paste0("SDA_", input$ComponentNtext, sep="")]), 
            aes(tSNE1, tSNE2, color=cut(asinh(SDAComp), breaks = c(-Inf, -1, -.5, 0, .5, 1, Inf)))) +
       geom_point(size=0.1) + theme_bw() +
       scale_color_manual("CS", values = rev(c("red", "orange", "yellow", "lightblue", "dodgerblue", "blue")) ) + 
       guides(colour = guide_legend(override.aes = list(size=2, alpha=1))) +
       theme(legend.position = "bottom", aspect.ratio=1,
             panel.background = element_rect(fill = "black",
                                             colour = "black",
                                             size = 0.5, linetype = "solid")) + 
       ggtitle(paste0("SDA_", input$ComponentNtext, " \n", 
                      StatFac[paste0("SDA_", input$ComponentNtext, sep=""),2], sep="")) + 
       simplify2 + coord_cartesian(xlim = NULL, ylim = NULL, expand = FALSE))
    
  })
  
  output$plot2 <- renderPlot({
    tempDF <- tDF()
    
    if(input$metaselect == "SubjectId") {
      MetaFac <- (MetaDF$SubjectId)
    } else {
        if(input$metaselect == "SampleDate"){
          MetaFac <- (MetaDF$SampleDate)
        } else {
          if(input$metaselect == "BarcodePrefix"){
            MetaFac <- (MetaDF$BarcodePrefix)
          } else {
            if(input$metaselect == "SingleR_Labels"){
              MetaFac <- (MetaDF$SingleR_Labels)
            }  else {
                  
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
      MetaFac <- (MetaDF$SubjectId)
    } else {
      if(input$metaselect == "SampleDate"){
        MetaFac <- (MetaDF$SampleDate)
      } else {
        if(input$metaselect == "BarcodePrefix"){
          MetaFac <- (MetaDF$BarcodePrefix)
        } else {
          if(input$metaselect == "SingleR_Labels"){
            MetaFac <- (MetaDF$SingleR_Labels)
          }  else {
                
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
    
    if(input$Genetext %in% colnames(GeneLoadings)){
      # GeneLoadings[,"PRM1"]
      GeneExpr <- CellScores %*% GeneLoadings[,as.character(input$Genetext)]
    } else {
      GeneExpr <- CellScores %*% rep(0, nrow(GeneLoadings))

    }
    #ggplotly

    LoadOrdVal <- round(GeneLoadings[,as.character(input$Genetext)][order(abs(GeneLoadings[,as.character(input$Genetext)]), decreasing = T)], 3)
    

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
      pgl <- genome_loadings(GeneLoadings[as.numeric(input$ComponentNtext),], 
                             label_both = T, 
                             max.items = 10, 
                             gene_locations =   gene_locations,
                             chromosome_lengths = chromosome.lengths, GeneLoadings = GeneLoadings)+ theme(aspect.ratio = .5)
      print(pgl)
      
    }
    
  })
  
  
  
  output$CellScoresAcross <- renderPlot({
    
    DTt <- ScoreAcrossDT()
    ComponentN <- as.numeric(input$ComponentNtext)
    # CellScores <- CellScores
    
    
    
    if(! (as.numeric(input$ComponentNtext) %in% 1:40)){
      print("No Comp")
    } else {
      lims <- quantile(asinh(CellScores[, paste0("SDA_", ComponentN)]^3), c(0.01, 0.99))
      
      pgl <- ggplot((DTt), 
                    aes(cell_index, asinh(score^3), colour = factor(ColFac))) + 
        geom_point(size = 0.7, stroke = 0) + 
        xlab("Cell Index") + ylab("asinh(Score^3)") + 
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
        ggtitle(paste0("SDA_", ComponentN)) +
        ylim(lims[1], lims[2])
      
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
    GeneSet <- GeneSet[GeneSet %in% colnames(GeneLoadings[,])]
    print("length of your genes in this dataset:")
    print(length(GeneSet))
    
    plotEnrich(GeneSetsDF=SDA_TopNpos, 
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

    GeneSet <- GeneSet[GeneSet %in% colnames(GeneLoadings[,])]

    
    plotEnrich(GeneSetsDF=SDA_TopNneg, 
               GeneVec = GeneSet, 
               plotTitle="Gene-set enrichment\n SDA top 40 neg loadings\n Cust. Input. genes \n Hypergeometric test: * adj.p < 0.01",
               xLab = "SDA Comps",
               N=N,
               k=k)
    
    
  })
  
  
}


shinyApp(ui, server)

