#Download the necessary libraries

library(shiny)
library(WGCNA)
library(ggplot2)

# The fluidPage contains sidebarPanel and mainPanel    
ui <- fluidPage(
  titlePanel("AHR IL4I1"),
  
  sidebarLayout(
    # sidebarPanel contains 3 main sections with possibility to choice object
    sidebarPanel(
      style = "overflow-y: scroll;position:fixed;width:inherit;height:100%;",
                 h5("Weighted correlation network analysis"),
                 fluidRow(column(selectize = T, width = 6,
                                 selectInput("varnet", 
                                             label="TCGA_TOMS", 
                                             choices = names(TCGA_TOMS_bicor),
                                             selected = "TCGA_ACC"))), 
                 h5("AHR_signature_development"),
                 fluidRow(column(selectize = T, width = 8,
                                 selectInput("GO", 
                                             label="processes", 
                                             choices = go_to_plot_BP2$GO_groups,
                                             selected = "Angiogenesis"))),
                 (h5("AHR-associated modules in 32 TCGA tumors")),
                 fluidRow(column(selectize = T, width = 12,
                   splitLayout(cellWidths = c("50%", "50%"), 
                                 radioButtons("circos1", 
                                             label="TCGA_AHR", 
                                             choices = df_img$img_path,
                                             selected = "TCGA_ACC"),
                                 radioButtons("circos2", 
                                             label="TCGA_IDO1/TDO2", 
                                             choices = df_img$img_path,
                                             selected = "TCGA_ACC"))
                 ))
      ),
    
    # mainPanel contains Voom plot to estimate counts in TCGA samples
    mainPanel("Generating a normalized counts DGElist using TMM normalization and voom",
              imageOutput("voom"),
             p("The eigengene dendrogram and heatmap identify groups of correlated
eigengenes - meta-modules"),
             # eigengen dendrogram and heatmap made by WGCNA package
              plotOutput("heatmap", height = 1000, width = 800),
             "AHR_GOs_circular_plot",
             # circos plot for estimatation of the process representation
             p("Gene ontology groups enriched in the pan-tissue AHR signature. The inner most circle represents the color code of each ontology group. Each bar represents
a significantly enriched ontology term. The bars are ordered in a descending order of highest significance in a clockwise fashion. The colors of each bar
correspond to the significance of enrichment. The length of each bar and the numbers in the outer circle represent the number of genes from the AHR signature
sharing the same ontology term"), 
             plotOutput("plot"),
              tableOutput("data"),
             # two plots on the same row for better comparison
             "Circos plots showing the weighted gene coexpression network analysis (WGCNA) modules positively associated with AHR activity. The left side
             plot shows connection between the enzymes and modules for the seven TCE encoding genes (IL4I1, IDO1, IDO2, TDO2, Trp hydroxylase 1 (TPH1), TPH2 or Dopa decarboxylase (DDC))
             The right side plot shows connection between the enzymes and modules for IDO1 or TDO2. The size of each module corresponds to the number of genes within the module.",
             fluidRow(
               splitLayout(cellWidths = c("50%", "50%"),imageOutput("AHRimage"),
                           imageOutput("IDOimage")))
              )
    ))
  
  
# Server part contains function for rendering images
server <- function(input, output) {
 #output for voom image
   output$voom <- renderImage({
    imgname <- normalizePath(file.path(paste(getwd(), "/www/", "voom.png", sep = "")))
    list(src = imgname)  
  }, deleteFile = FALSE)
    # this output uses plotEigengeneNetworks command to build heatmap and dendrogramm based on WGCNA topological overlap matrix 
  output$heatmap <- renderPlot({
    plotEigengeneNetworks(TCGA_TOMS_bicor[[input$varnet]]$MEs,"Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                          plotDendrograms = T, xLabelsAngle = 90)
    
  })  
  # this output for plot and table uses the script built by A.Sadik 
  output$plot <- renderPlot({
    ggplot(go_to_plot_BP2, aes(x=pid, y=Count, fill=logpv, group=GO_groups))+
      geom_bar(stat="identity", position = position_stack(reverse = TRUE), alpha=0.5)+ ylim(c(-30, 40))+
      theme_minimal()+ xlab("BP")+
      scale_fill_gradientn(colours = colorRampPalette(c("#1b9e77","#d95f02","#7570b3"))(99))+
      theme(axis.text = element_blank(),
            axis.title = element_blank(),
            panel.grid = element_blank())+ coord_polar()+
      geom_text(data=label_data, aes(x=pid, y=36, label=Count, hjust=hjust), color="black", fontface="bold",alpha=0.6, angle= label_data$angle, inherit.aes = FALSE )+
      geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5, color=GO_groups), alpha=0.8, size=2.6 , inherit.aes = FALSE, show.legend = TRUE) + 
      scale_colour_brewer(palette = "Dark2", guide=guide_legend(title = "Biological processes GOs"))
  })
  
  output$data <- renderTable({
    na.omit(go_to_plot_BP2[which(go_to_plot_BP2$GO_groups == input$GO), ])
  })
  # the images are the result of AHR signature ontology analysis mady by A.Sadik
  output$AHRimage <- renderImage({
    filename <- normalizePath(file.path(paste(getwd(), "/www/WGCNA_circos_AHR/", input$circos1, ".jpg", sep = "")))
    list(src = filename,
         width = 300,
         height= 300)  
  }, deleteFile = FALSE)
  output$IDOimage <- renderImage({
    filename <- normalizePath(file.path(paste(getwd(), "/www/WGCNA_circos_IDO1_TDO2/", input$circos2, ".jpg", sep = "")))
    list(src = filename,
         width = 300,
         height= 300)  
  }, deleteFile = FALSE)
  
}
shinyApp(ui, server)