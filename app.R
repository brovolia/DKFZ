library(shiny)

base_data <- readRDS("./objects_shiny/base_data.RDS")
go_to_plot_BP2 <- readRDS("./go_to_plot_BP2.RDS")
label_data <- readRDS("./objects_shiny/label_data.RDS")
overlapping_genes <- readRDS("./objects_shiny/overlapping_genes.RDS")
to_add <- readRDS("./objects_shiny/to_add.RDS")

ui <- fluidPage(
    titlePanel("AHR IL4I1"),
    
    sidebarLayout(
        sidebarPanel(
            style = "position:fixed;width:inherit;height:120%;overflow-y: scroll;",
            h3("TCGA_data"),
            h3("AHR_signature_development"),
            fluidRow(column(selectize = T, width = 8,
                            selectInput("GO", 
                                        label="processes", 
                                        choices = go_to_plot_BP2$GO_groups,
                                        selected = "Drug metabolism")))),
        mainPanel("Generating a normalized counts DGElist using TMM normalization and voom",
                  img(src = "voom.png", height = 372, width = 572),
                  "AHR_GOs_circular_plot",
                  plotOutput("plot"),
                  tableOutput("data")
        )
    ))





server <- function(input, output) {
    
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
}
shinyApp(ui, server)