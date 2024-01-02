library(glue)
library(ActivePathways)
library(stringr)
library("RColorBrewer") # https://plotly.com/r/builtin-colorscales/
library(shiny)
library(shinythemes)
library(DT)
library(plotly)
library(assertthat)
library(readr)
library(filesstrings)
linebreaks <- function(n){HTML(strrep(br(), n))}


######################################
############ DIRECTORIES ############# 
######################################
getwd()

# 01 - Read directory for ATAC peak assigned to gene datasets 
peakAnno_df_dir <- "./results/peakAnno_df"
# 02 - Datasets with logFC by gene: gene level
logFC_by_gene_df_dir <- "./results/logFC_by_gene_df"
# 03 - Datasets with Avg ATAC logFC, # sig peaks and genes by GS (Gene Set): pathway level
logFC_by_GS_df_dir <- "./results/logFC_by_GS_df"
# Read geneset collection metadata directory
geneset_collection_metadata_dir <- "./data/geneset_collection_metadata/"

# Experiment metadata
metadata_df_path <-  './data/experiment_metadata_df.csv'
metadata_df <- read.table(metadata_df_path, header=TRUE, sep =',', check.names = F)
metadata_df <- metadata_df[complete.cases(metadata_df), ]

# Geneset collection metadata  
read_GSEA_collection_mapping <- function(geneset_collection_metadata_dir){  # Read GSEA on RNA DGE logFC
  # Add family pathway to GSEA report 
  GSEA_collection_mapping_path <- paste0(geneset_collection_metadata_dir, "msigdb_fht_other_genesets_v1_metadata_r38109x8.txt.gz")
  print(paste('Assert that',  GSEA_collection_mapping_path, 'exists:', assert_that(file.exists(GSEA_collection_mapping_path))))
  GSEA_collection_mapping <- read_delim(GSEA_collection_mapping_path, delim ='\t')
  # GSEA_collection_mapping <- GSEA_collection_mapping[, -which(names(GSEA_collection_mapping) %in% c("genes", "geneset_size"))]
  return(GSEA_collection_mapping)
}
GSEA_collection_mapping <- read_GSEA_collection_mapping(geneset_collection_metadata_dir)
read_genesets_boolean_matrix <- function(geneset_collection_metadata_dir){
  memory.limit() 
  memory.limit(size = 35000)
  genesets_boolean_matrix_path <- paste0(geneset_collection_metadata_dir, "msigdb_fht_other_genesets_v1_boolean_matrix_r42673x38109.rds")
  print(paste('Assert that',  genesets_boolean_matrix_path, 'exists:', assert_that(file.exists(genesets_boolean_matrix_path))))
  genesets_boolean_matrix <- readRDS(genesets_boolean_matrix_path)
  return(genesets_boolean_matrix)
}


######################################
################# UI ################# 
######################################
ui <- fluidPage(theme = shinytheme("flatly"),
                  # theme = "cerulean",  # <--- To use a theme, uncomment this
                  titlePanel("XChrom: Linking gene expression to chromatin profiling"),
                  sidebarLayout(
                    sidebarPanel(
                             h3("Select your options for the scatter plot"),
                             linebreaks(3),
                             selectInput("experiment", "Select an experiment",
                                          levels(as.factor(metadata_df$title)),
                                         selected = levels(as.factor(metadata_df$title))[1]),
                             br(),
                             strong("Experiment Info"),
                             htmlOutput("Experiment_Info"),
                             # br(),
                             # strong("RNA experiment Info"),
                             # htmlOutput("RNA_description"),
                             linebreaks(5),
                             # selectizeInput("geneset","Select a pathway collection",  choices=NULL),
                             selectizeInput("geneset","Select a gene set collection",
                                         choices = unique(GSEA_collection_mapping$collection_display_name),
                                         selected = unique(GSEA_collection_mapping$collection_display_name)[1]),
                             
                             # h5(HTML("Project <b>Description</b>"), 
                             #    style="text-align:center"))),
                             # strong("Genset Collection Info"),
                             htmlOutput("geneset_description"),
                             linebreaks(5),
                             selectizeInput("pathway","Select a gene set from the geneset collection from above", choices=NULL),
                             # strong("Genset Info"),
                             htmlOutput("pathway_description"),
                             width = 3,
                             height = 19
                    
                           ), # sidebarPanel
                    
                    
                           mainPanel(
                             tabsetPanel(type = "tabs",
                             tabPanel("Gene", 
                                      h3("Gene level: expression vs. chromatin dynamics"),
                                      br(),
                                      htmlOutput("full_gene_summary_description"),
                                      br(),
                                      DT::dataTableOutput('full_gene_summary_table',  width = "1400px"),
                                      linebreaks(3),
                                      
                                      
                                      htmlOutput("chosen_gene_description"),
                                      br(),
                                      selectizeInput("gene", "", choices=NULL),
                                      br(),
                                      plotlyOutput(outputId = "gene_scatter_plot", width = "1400px", height = "600px"),
                                      linebreaks(3),
                                      
                                      
                                      htmlOutput("chosen_gene_peak_summary_description"),
                                      br(),
                                      radioButtons("chosen_gene_filter_sig_peaks", "Filter significant peaks",
                                                   c("YES",
                                                     "NO"),
                                                   inline = TRUE),
                                      br(),
                                      DT::dataTableOutput('chosen_gene_peak_summary_table',  width = "1400px"),
                                      linebreaks(3),
                                      
                                      htmlOutput("chosen_gene_peak_id"),
                                      linebreaks(3),
                                      
                                      
                                      
                                      htmlOutput("display_geneset_table_description"),
                                      br(),
                                      radioButtons("chosen_gene_pathway_list_display", "Display table?",
                                                   c("NO",
                                                     "YES"),
                                                   inline = TRUE),
                                      br(),
                                      DT::dataTableOutput('chosen_gene_geneset_table',  width = "1400px"),
                                      linebreaks(3)
                                      ), # Gene Tab
                             
                             tabPanel("Pathway",
                                      
                                      h3("Pathway level: expression vs. chromatin dynamics"),
                                      
                                      div("The x axis represents the change in open chromatin. It is the average logFC from any significant ATAC peaks
proximal to pathway genes weighted by the number of significant genes, i.e., avg(ATAC logFC) * log10(nb of sig genes)."),
                                      div("The y axis represents the change in gene expression. It is the weighted p-value NES, i.e. -log10(pval) * NES."),
                                      div("Each point represents a geneset from the selected gene set collection. Click on one of them to generate a summary table below."),
                                      plotlyOutput(outputId = "scatter_plot", width = "1400px", height = "600px"),
                                      linebreaks(3),
                                      
                                      htmlOutput("gene_summary_description"),
                                      br(),
                                      DT::dataTableOutput('gene_summary_table',  width = "1400px"), #, height = "400px"),
                                      linebreaks(3),
                                      
                                      htmlOutput("peak_summary_description"),
                                      br(),
                                      radioButtons("filter_sig_peaks", "Filter significant peaks",
                                                   c("YES",
                                                     "NO"),
                                                   inline = TRUE),
                                      br(),
                                      DT::dataTableOutput('peak_summary_table',  width = "1400px"),
                                      linebreaks(3),
                                      
                                      htmlOutput("peak_id"),
                                      linebreaks(3)
                                      
                                      # htmlOutput('gene_browser')
                             ), # Pathway Tab
                             
                             ) # tabsetPanel
                           ) # mainPanel
                  ) # sidebarPanel
                  
                  # ) tabPanel

) # fluidPage


######################################
############### SERVER ############### 
######################################
server <- function(input, output,session) {


  # https://stackoverflow.com/questions/56273889/errors-with-plotly-click-in-r-plotly-version-4-9-0-is-there-a-bug-in-the-new-ve
  reactiveData <- reactiveValues()

  toListen <- reactive({
    list(input$experiment, input$geneset)
  }) # toListen
  #

  observeEvent(toListen(),{
    reactiveData$experiment <- input$experiment
    print(paste('reactiveData$experiment:', reactiveData$experiment))

    ### Read peakAnno_df: ATAC peak assigned to gene datasets (Table2) ###
    peakAnno_df_name <- metadata_df[metadata_df$title == reactiveData$experiment,]$peakAnno_df
    peakAnno_df_path <- paste0(glue("{peakAnno_df_dir}/"), peakAnno_df_name)
    print(paste('peakAnno_df path:', peakAnno_df_path))
    peakAnno_df <- read.table(peakAnno_df_path, header=TRUE, sep =',')
    reactiveData$peakAnno_df <- peakAnno_df 
    
    ### Read logFC_by_gene_df: logFC by gene (Table 1) ###
    logFC_by_gene_df_name <- metadata_df[metadata_df$title == reactiveData$experiment,]$logFC_by_gene_df
    logFC_by_gene_df_path <- paste0(glue("{logFC_by_gene_df_dir}/"), logFC_by_gene_df_name)
    print(paste('logFC_by_gene_df path:', logFC_by_gene_df_path))
    logFC_by_gene_df <- read.table(logFC_by_gene_df_path, header=TRUE, sep =',')
    reactiveData$logFC_by_gene_df <- logFC_by_gene_df

    ### Read logFC_by_GS_df: Avg ATAC logFC, # sig peaks and genes by GS (plot) ###
    logFC_by_GS_df_name <- metadata_df[metadata_df$title == reactiveData$experiment,]$logFC_by_GS_df
    logFC_by_GS_df_path <- paste0(glue("{logFC_by_GS_df_dir}/"), logFC_by_GS_df_name)
    print(paste('logFC_by_GS_df path:', logFC_by_gene_df_path))
    logFC_by_GS_df <- read.table(logFC_by_GS_df_path, header=TRUE, sep =',')

    reactiveData$geneset <- input$geneset
    print(paste('reactiveData$geneset:', reactiveData$geneset))
    reactiveData$logFC_by_GS_df <- logFC_by_GS_df[logFC_by_GS_df$collection_display_name == reactiveData$geneset, ] 
    updateSelectizeInput(session = session,
                         "pathway", glue("Select a gene set from ", reactiveData$geneset),
                         unique(reactiveData$logFC_by_GS_df$geneset))
    
    
    updateSelectizeInput(session = session,
                         "gene", "Select a gene to highlight it in the plot below: ",
                         unique(reactiveData$logFC_by_gene_df$RNA_SYMBOL))
  }) # observeEvent 
  ################################################################

  # Listen to several events for the same variable
  # https://groups.google.com/g/shiny-discuss/c/0V4sR4LjAuc
  observe({
    if(!is.null(input$pathway)){
      reactiveData$chosen_pathway <- isolate(input$pathway)
    }
  })

  observe({
    if(!is.null(event_data("plotly_click", source='scatter_plot'))){

      d <- event_data("plotly_click", source='scatter_plot')
      if (is.null(d)) return(shiny::showNotification("No data", type = "error"))
      reactiveData$chosen_pathway  <- isolate(d$customdata)
    }
  })
  
  
  observe({
    if(!is.null(input$gene)){
      reactiveData$chosen_gene <- isolate(input$gene)
    }
  })
  
  observe({
    if(!is.null(input$full_gene_summary_table_cell_clicked$value)){
      reactiveData$chosen_gene <- isolate(input$full_gene_summary_table_cell_clicked$value)
    }
  })
  
  observe({
    if(!is.null(event_data("plotly_click", source='gene_scatter_plot'))){
      
      d <- event_data("plotly_click", source='gene_scatter_plot')
      if (is.null(d)) return(shiny::showNotification("No data", type = "error"))
      reactiveData$chosen_gene  <- isolate(d$customdata)
    }
  })
  
  ################ experiment_description ##########################
  output$Experiment_Info <- renderUI({
    title <- metadata_df[metadata_df$title == reactiveData$experiment,]$title
    firewall_status <- metadata_df[metadata_df$title == reactiveData$experiment,]$firewall_status
    bio_context_id <- metadata_df[metadata_df$title == reactiveData$experiment,]$bio_context_id
    lineage <- metadata_df[metadata_df$title == reactiveData$experiment,]$lineage
    pert_type <- metadata_df[metadata_df$title == reactiveData$experiment,]$pert_type
    time_point <- metadata_df[metadata_df$title == reactiveData$experiment,]$time_point
    pert_dose <- metadata_df[metadata_df$title == reactiveData$experiment,]$pert_dose
    RNA_experiment_id <- metadata_df[metadata_df$title == reactiveData$experiment,]$RNA_experiment_id
    ATAC_experiment_id <- metadata_df[metadata_df$title == reactiveData$experiment,]$ATAC_experiment_id
    Genome <- metadata_df[metadata_df$title == reactiveData$experiment,]$Genome
    filtered_txdb <- metadata_df[metadata_df$title == reactiveData$experiment,]$filtered_txdb
    Contrast <- metadata_df[metadata_df$title == reactiveData$experiment,]$Contrast
    description <- metadata_df[metadata_df$title == reactiveData$experiment,]$description

    HTML(paste('<b>Title: </b>', title,  "<br>",
               '<b>Firewall status: </b>', firewall_status,  "<br>",
               '<b>bio_context_id: </b>', bio_context_id,  "<br>", 
               '<b>Lineage: </b>', lineage,  "<br>", 
               '<b>pert_type: </b>',  pert_type,  "<br>", 
               '<b>time_point: </b>', time_point,  "<br>", 
               '<b>pert_dose: </b>', pert_dose,  "<br>", 
               '<b>RNA_experiment_id: </b>', RNA_experiment_id,  "<br>", 
               '<b>ATAC_experiment_id: </b>', ATAC_experiment_id,  "<br>", 
               '<b>Genome: </b>', Genome,  "<br>", 
               '<b>Filtered assigned peak to gene: </b>', filtered_txdb, "<br>", 
               '<b>Contrast: </b>', Contrast, "<br>",
               '<b>Description: </b>', description))
  })
  ################################################################


  ################ geneset_description ###########################
  output$geneset_description <- renderUI({
    HTML(paste(unique(reactiveData$logFC_by_GS_df[reactiveData$logFC_by_GS_df$collection_display_name== reactiveData$geneset,]$collection_description),
               paste0("<a href=", unique(reactiveData$logFC_by_GS_df[reactiveData$logFC_by_GS_df$collection_display_name== reactiveData$geneset,]$collection_source_url), ">", paste(reactiveData$geneset, "- Collection source URL"), "</a>"),
               sep = "<br/>"))
  })
  ################################################################


  ################# pathway_description ###########################
  output$pathway_description <- renderUI({
    HTML(paste(paste0("<a href=", unique(reactiveData$logFC_by_GS_df[reactiveData$logFC_by_GS_df$geneset== reactiveData$chosen_pathway,]$description), ">", paste(reactiveData$chosen_pathway, "- Pathway source URL"), "</a>"),
               sep = "<br/>"))
  })
  #################################################################

  
  
  
  
  
  
  
  
  
  ################ TAB 1: GENE LEVEL #############################
  
  
  ################ gene_summary table ############################
  full_gene_summary <- reactive({
    
    full_gene_summary <- req(reactiveData$logFC_by_gene_df)
    
    # # Subset columns
    full_gene_summary <- full_gene_summary[, c("RNA_SYMBOL", "RNA_logFC", "RNA_P_Value", 'ATAC_SYMBOL', "Avg_ATAC_logFC", "nb_sig_peak", "GENENAME", "BAF.gene", "TF.gene")]
    # Remove duplicate rows
    full_gene_summary <- unique(full_gene_summary)
    
    # convert values as numeric
    full_gene_summary$RNA_logFC <- as.numeric(full_gene_summary$RNA_logFC)
    full_gene_summary$RNA_P_Value <- as.numeric(full_gene_summary$RNA_P_Value)
    full_gene_summary$Avg_ATAC_logFC <- as.numeric(full_gene_summary$Avg_ATAC_logFC)
    full_gene_summary$nb_sig_peak <- as.numeric(full_gene_summary$nb_sig_peak)
    # round number of decimal
    full_gene_summary$RNA_logFC <- round(full_gene_summary$RNA_logFC, 5)
    full_gene_summary$RNA_P_Value <- round(full_gene_summary$RNA_P_Value, 5)
    full_gene_summary$Avg_ATAC_logFC <- round(full_gene_summary$Avg_ATAC_logFC, 5)
    # convert values as factor
    full_gene_summary$BAF.gene <- as.factor(full_gene_summary$BAF.gene)
    full_gene_summary$TF.gene <- as.factor(full_gene_summary$TF.gene)
    
    full_gene_summary <- full_gene_summary[order(full_gene_summary$RNA_SYMBOL),]
    full_gene_summary
    
  })
  # full_gene_summary()
  ################################################################
  
  
  ################ full_gene_summary_description #################
  output$full_gene_summary_description <- renderUI({
    req(input$experiment)
    chosen_gene = reactiveData$chosen_gene
    HTML(paste(paste('<font size="+2">Summary table of the genes in', '<b>', input$experiment, '</b> .</font>'),
               'Click on one of the gene names in the table to highlight it in the plot below.',
               sep = "<br/>"))
  })
  ################################################################
  
  
  ################ full_gene_summary output table ################
  output$full_gene_summary_table <- DT::renderDT(full_gene_summary(), selection = "single", filter = 'top',
                                                 options = list(
                                                   autoWidth = TRUE
                                                   # selection = "single",
                                                   # initComplete = JS(
                                                   #   "function(settings, json) {",
                                                   #   "$(this.api().table().header()).css({'background-color': '#000', 'color': '#fff'});",
                                                   #   "}")
                                                 )) # output$gene_summary_table
  ################################################################
  
  
  ################ chosen_gene_description #######################
  output$chosen_gene_description <- renderUI({
    req(reactiveData$chosen_gene)
    chosen_gene = reactiveData$chosen_gene
    HTML(paste(paste('<font size="+2">Selected gene from table above or dropdown is:', '<b>', chosen_gene, '</b>.</font>'),
               'The x axis represents the change in open chromatin. It is the average logFC from any significant ATAC peaks
proximal to pathway genes, i.e., avg(ATAC logFC).',
               'The y axis represents the change in gene expression, i.e. RNA logFC.',
               'Each point represents a gene listed in the table above from the selected experiment. Click on one of the points to generate the table below.',
               sep = "<br/>"))
  })
  ################################################################
  
  
  ################ gene_scatter_plot #############################
  output$gene_scatter_plot <- renderPlotly({
    req(reactiveData$logFC_by_gene_df)
    req(reactiveData$chosen_gene)
    logFC_by_gene_df <- reactiveData$logFC_by_gene_df
    chosen_gene <- reactiveData$chosen_gene
    
    fig_gene <- plot_ly(source = "gene_scatter_plot",
                        name = "all genes",
                        data = logFC_by_gene_df,
                        customdata = ~RNA_SYMBOL,
                        x = ~Avg_ATAC_logFC,
                        y = ~RNA_logFC,
                        color = ~log10(nb_sig_peak + 1),
                        colors = "YlGnBu",
                        text = ~paste("<b>Gene</b>: ", RNA_SYMBOL,
                                      '\n <b>RNA logFC</b>: ', RNA_logFC,
                                      '\n <b>Avg ATAC logFC</b>: ', Avg_ATAC_logFC,
                                      '\n <b># sig peaks in ATAC</b>: ', nb_sig_peak
                        ),
                        hoverinfo = 'text')%>%
      colorbar(title = "Nb sig peaks")%>%
      layout(xaxis = list(title = 'change in open chromatin'),
             yaxis = list(title = 'change in gene expression'))  %>%
      event_register("plotly_click")
    
    
    fig_gene <- fig_gene %>% add_trace(
      inherit = FALSE,
      name = chosen_gene,
      x = c(logFC_by_gene_df[logFC_by_gene_df$RNA_SYMBOL == chosen_gene,]$Avg_ATAC_logFC),
      y = c(logFC_by_gene_df[logFC_by_gene_df$RNA_SYMBOL == chosen_gene,]$RNA_logFC),
      text = ~paste("<b>Gene</b>: ", logFC_by_gene_df[logFC_by_gene_df$RNA_SYMBOL == chosen_gene,]$RNA_SYMBOL,
                    '\n <b>RNA logFC</b>: ', logFC_by_gene_df[logFC_by_gene_df$RNA_SYMBOL == chosen_gene,]$RNA_logFC,
                    '\n <b>Avg ATAC logFC</b>: ', logFC_by_gene_df[logFC_by_gene_df$RNA_SYMBOL == chosen_gene,]$Avg_ATAC_logFC,
                    '\n <b># sig peaks in ATAC</b>: ', logFC_by_gene_df[logFC_by_gene_df$RNA_SYMBOL == chosen_gene,]$nb_sig_peak),
      hoverinfo = 'text',
      marker=list(
        showscale=FALSE,
        line=list(
          color='orange',
          width=2,
          opacity = 1
        ))
    )%>%
      layout(showslegend=FALSE)
    
    fig_gene %>%
      event_register("plotly_click")
  }) # output$gene_scatter_plot
  ################################################################
  
  
  
  ################ peak_summary table ############################
  chosen_gene_peak_summary <- reactive({
    req(reactiveData$chosen_gene)
    chosen_gene <- reactiveData$chosen_gene
    req(reactiveData$peakAnno_df)
    peakAnno_df <- reactiveData$peakAnno_df
    
    filter_sig_peaks <- input$chosen_gene_filter_sig_peaks
    if (filter_sig_peaks == 'YES'){filter = 0.05}
    else { filter = 1}
    
    subset_gene_df <- peakAnno_df[peakAnno_df$SYMBOL==chosen_gene & !is.na(peakAnno_df$SYMBOL) & peakAnno_df$P_Value < filter,
                                  c('logFC', 'P_Value', 'distanceToTSS', 'peak_id', 'annotation')]
    
    # convert values as numeric
    subset_gene_df$logFC <- as.numeric(subset_gene_df$logFC)
    subset_gene_df$P_Value <- as.numeric(subset_gene_df$P_Value)
    subset_gene_df$distanceToTSS <- as.numeric(subset_gene_df$distanceToTSS)
    # round number of decimal
    subset_gene_df$logFC <- round(subset_gene_df$logFC, 5)
    subset_gene_df$P_Value <- round(subset_gene_df$P_Value, 5)
    
    subset_gene_df
  })
  ################################################################
  
  
  
  ################ peak_summary_description ######################
  output$chosen_gene_peak_summary_description <- renderUI({
    
    req(reactiveData$chosen_gene)
    chosen_gene <- reactiveData$chosen_gene
    
    HTML(paste(paste('<font size="+2">Summary table of the peaks from ', '<b>', chosen_gene, '</b></font>'),
               'Click on one of the peak_id in the table below to view the genome borowser in that region.',
               sep = "<br/>"))
  })
  ################################################################
  
  
  
  ################ peak_summary output table ######################
  output$chosen_gene_peak_summary_table <- DT::renderDT(chosen_gene_peak_summary(), selection = "single", filter = 'top',
                                                        options = list(
                                                          autoWidth = TRUE
                                                        )) # output$peak_summary_table
  ################################################################
  
  
  
  ################ peak_id from the peak_summary table ###########
  output$chosen_gene_peak_id <- renderUI({
    
    
    req(reactiveData$chosen_gene)
    chosen_gene <- reactiveData$chosen_gene
    peak_id = input$chosen_gene_peak_summary_table_cell_clicked
    URL = paste0("<a href=http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=", peak_id$value,">",
                 paste0('Open the UCSC Gene browser for the gene ', chosen_gene, ' at the region ', peak_id$value),"</a>")
    HTML(paste(paste('<font size="+2">UCSC Genome Browser on Human (GRCh38/hg38) for ', '<b>', peak_id$value, '</b></font>'),
               URL,
               sep = "<br/>"))
  })
  ################################################################
  
  
  ################ chosen_gene_geneset_list ######################
  chosen_gene_geneset_list <- reactive({
    req(input$chosen_gene_pathway_list_display)
    req(reactiveData$chosen_gene)
    chosen_gene_pathway_list_display <- input$chosen_gene_pathway_list_display
    chosen_gene = reactiveData$chosen_gene
    
    # print('chosen_gene_geneset_question: ', chosen_gene_pathway_list_display)
    if (chosen_gene_pathway_list_display == 'NO'){
      chosen_gene_df = data.frame(matrix(ncol = 3, nrow = 0))
      colnames(chosen_gene_df) <- c('gene set collection', 'gene set' , 'description')
    }
    else {
      if (exists("all_genesets_boolean_matrix")){
        chosen_gene_pathway_list = names(all_genesets_boolean_matrix[chosen_gene, ][all_genesets_boolean_matrix[chosen_gene, ]==1])
        chosen_gene_pathway_df = data.frame(chosen_gene_pathway_list)
        chosen_gene_df = merge(x = GSEA_collection_mapping, y = chosen_gene_pathway_df,
                               by.x ='cleaned_name', by.y = "chosen_gene_pathway_list", all.y=T)[, c("collection_display_name", "cleaned_name", "collection_description")]
        colnames(chosen_gene_df) <- c('gene set collection', 'gene set' , 'description')
      }
      else{
        all_genesets_boolean_matrix <<- read_genesets_boolean_matrix(geneset_collection_metadata_dir)
        chosen_gene_pathway_list = names(all_genesets_boolean_matrix[chosen_gene, ][all_genesets_boolean_matrix[chosen_gene, ]==1])
        chosen_gene_pathway_df = data.frame(chosen_gene_pathway_list)
        chosen_gene_df = merge(x = GSEA_collection_mapping, y = chosen_gene_pathway_df,
                               by.x ='cleaned_name', by.y = "chosen_gene_pathway_list", all.y=T)[, c("collection_display_name", "cleaned_name", "collection_description")]
        colnames(chosen_gene_df) <- c('gene set collection', 'gene set' , 'description')
      }
    }
    chosen_gene_df
  }) #chosen_gene_genset_list
  ################################################################
  
  
  ################ display_geneset_table_description #############
  output$display_geneset_table_description <- renderUI({
    req(reactiveData$chosen_gene)
    chosen_gene = reactiveData$chosen_gene
    HTML(paste(paste('<font size="+2">Summary table of the gene set collections that include', '<b>', chosen_gene, '</b>?</font>'),
               'It might take a few seconds to load.',
               sep = "<br/>"))
  })
  ################################################################
  
  
  ################ chosen_gene_geneset_list output table #############################
  output$chosen_gene_geneset_table <- DT::renderDT(chosen_gene_geneset_list(), selection = "single", filter = 'top',
                                                   options = list(
                                                     autoWidth = TRUE
                                                     # selection = "single",
                                                     # initComplete = JS(
                                                     #   "function(settings, json) {",
                                                     #   "$(this.api().table().header()).css({'background-color': '#000', 'color': '#fff'});",
                                                     #   "}")
                                                   )) # output$gene_summary_table
  ################################################################
  
  
  
  
  
  
  
  
  ################ TAB 2: GENE SET/PATHWAY LEVEL ############
  
  

  ################ scatter_plot #############################
  output$scatter_plot <- renderPlotly({
    req(reactiveData$logFC_by_GS_df)
    req(reactiveData$chosen_pathway)
    logFC_by_GS_df <- reactiveData$logFC_by_GS_df
    chosen_pathway <- reactiveData$chosen_pathway

    fig <- plot_ly(source = "scatter_plot",
                   name = "all pathways",
                   data = logFC_by_GS_df,
                   x = ~Avg_ATAC_logFC * log10(nb_sig_gene + 1e-6),
                   # x = ~sign(Avg_ATAC_logFC) * nb_sig_gene * 100 / geneset_size,
                   y = ~NES_RNA_score,
                   customdata = ~geneset,
                   color = ~log10(nb_sig_peak+1),
                   # reversescale=T,
                   colors = "YlGnBu",# "YlOrRd",
                   text = ~paste("<b>Pathway</b>: ", geneset,
                                 # '\n <b># matched genes</b>: ', matched_size,
                                 '\n <b>gene set size</b>: ', geneset_size,

                                 "\n\n <b>GSEA pval (RNA)</b>: ", round(pval, 5),
                                 '\n <b>GSEA NES (RNA)</b>: ', round(nes, 5),

                                 '\n\n <b># sig genes in ATAC</b>:', nb_sig_gene,
                                 '\n <b># sig peaks in ATAC</b>:', nb_sig_peak
                   ),
                   hoverinfo = 'text')%>%
      colorbar(title = "Nb sig peaks")%>%
      layout(xaxis = list(title = 'weighted change in open chromatin'),
             yaxis = list(title = 'weighted change in pathway expression'))  %>%
      event_register("plotly_click")

    hline <- function(y = 0, color = "red") {
      list(
        type = "line",
        x0 = 0,
        x1 = 1,
        xref = "paper",
        y0 = y,
        y1 = y,
        line = list(color = color, dash="dot")
      )
    }

    vline <- function(x = 0, color = "green") {
      list(
        type = "line",
        y0 = 0,
        y1 = 1,
        yref = "paper",
        x0 = x,
        x1 = x,
        line = list(color = color, dash="dot")
      )
    }

    fig <- fig %>%
      layout(plot_bgcolor = "#e5ecf6", shapes = list(hline(-2), hline(2)))

    # https://plotly.com/python/marker-style/

    fig <- fig %>% add_trace(
      inherit = FALSE,
      name = chosen_pathway,
      x = c(logFC_by_GS_df[logFC_by_GS_df$geneset  == chosen_pathway,]$Avg_ATAC_logFC * log10(logFC_by_GS_df[logFC_by_GS_df$geneset == chosen_pathway,]$nb_sig_gene)),
      y = c(logFC_by_GS_df[logFC_by_GS_df$geneset  == chosen_pathway,]$NES_RNA_score),
      text = ~paste("<b>Pathway</b>: ", logFC_by_GS_df[logFC_by_GS_df$geneset  == chosen_pathway,]$geneset,
                    # '\n <b># matched genes</b>: ', logFC_by_GS_df[logFC_by_GS_df$geneset  == chosen_pathway,]$matched_size,
                    # '\n <b>gene set size</b>: ', logFC_by_GS_df[logFC_by_GS_df$geneset  == chosen_pathway,]$geneset_size,

                    "\n\n <b>GSEA pval (RNA)</b>: ", round(logFC_by_GS_df[logFC_by_GS_df$geneset  == chosen_pathway,]$pval, 5),
                    '\n <b>GSEA NES (RNA)</b>: ', round(logFC_by_GS_df[logFC_by_GS_df$geneset  == chosen_pathway,]$nes, 5),

                    '\n\n <b># sig genes in ATAC</b>:', logFC_by_GS_df[logFC_by_GS_df$geneset  == chosen_pathway,]$nb_sig_gene,
                    '\n <b># sig peaks in ATAC</b>:', logFC_by_GS_df[logFC_by_GS_df$geneset  == chosen_pathway,]$nb_sig_peak
      ),
      hoverinfo = 'text',
      marker=list(
        showscale=FALSE,
        line=list(
          color='orange',
          width=2,
          opacity = 1
        ))
    )%>%
      layout(showslegend=FALSE)

    fig %>%
      event_register("plotly_click")
  }) # output$scatter_plot
  ################################################################

  
  ################ gene_summary table #############################
  gene_summary_table <- reactive({
    if(is.null(reactiveData$chosen_pathway)){
      print(paste('reactiveData$chosen_pathway: ', reactiveData$chosen_pathway))
      pathway_summary_table <- NULL
    } else {
      req(reactiveData$logFC_by_GS_df)
      req(reactiveData$logFC_by_gene_df)
      req(reactiveData$chosen_pathway)
      logFC_by_GS_df <- reactiveData$logFC_by_GS_df
      logFC_by_gene_df <- reactiveData$logFC_by_gene_df
      chosen_pathway<- reactiveData$chosen_pathway

      # Get gene list from chosen pathway
      gene_chosen_pathway <- logFC_by_GS_df[logFC_by_GS_df$geneset == chosen_pathway,]$genes
      list_gene_chosen_pathway <- strsplit(gene_chosen_pathway, split=';')[[1]]
      # Subset genes in logFC_by_gene_df
      pathway_summary_table <- logFC_by_gene_df[logFC_by_gene_df$RNA_SYMBOL %in% list_gene_chosen_pathway, ]
      # # Subset columns
      pathway_summary_table <- pathway_summary_table[, c("RNA_SYMBOL", "RNA_logFC", "RNA_P_Value", 'ATAC_SYMBOL', "Avg_ATAC_logFC", "nb_sig_peak", "GENENAME", "BAF.gene", "TF.gene")]
      names(pathway_summary_table)[names(pathway_summary_table) == 'RNA_SYMBOL'] <- 'Gene'
      # Remove duplicate rows
      pathway_summary_table <- unique(pathway_summary_table)

      # convert values as numeric
      pathway_summary_table$RNA_logFC <- as.numeric(pathway_summary_table$RNA_logFC)
      pathway_summary_table$RNA_P_Value <- as.numeric(pathway_summary_table$RNA_P_Value)
      pathway_summary_table$Avg_ATAC_logFC <- as.numeric(pathway_summary_table$Avg_ATAC_logFC)
      pathway_summary_table$nb_sig_peak <- as.numeric(pathway_summary_table$nb_sig_peak)
      # round number of decimal
      pathway_summary_table$RNA_logFC <- round(pathway_summary_table$RNA_logFC, 5)
      pathway_summary_table$RNA_P_Value <- round(pathway_summary_table$RNA_P_Value, 5)
      pathway_summary_table$Avg_ATAC_logFC <- round(pathway_summary_table$Avg_ATAC_logFC, 5)
      # convert values as factor
      pathway_summary_table$BAF.gene <- as.factor(pathway_summary_table$BAF.gene)
      pathway_summary_table$TF.gene <- as.factor(pathway_summary_table$TF.gene)
     
    }
    validate(
      need(!is.null(pathway_summary_table), "No data")
    )
    return(pathway_summary_table)
  })
  # gene_summary_table()
  ################################################################

  
  ################ gene_summary_description ######################
  output$gene_summary_description <- renderUI({
    HTML(paste(paste('<font size="+2">Summary table of the genes in ', '<b>', reactiveData$chosen_pathway, '</b></font>'),
               'Click on one of the gene names in the table below to view the peaks asssociated to this gene.',
               sep = "<br/>"))
  })
  ################################################################


  ################ gene_summary output table #####################
  output$gene_summary_table <- DT::renderDT(gene_summary_table(), selection = "single", filter = 'top',
                                            options = list(
                                              autoWidth = TRUE
                                              # selection = "single",
                                              # initComplete = JS(
                                              #   "function(settings, json) {",
                                              #   "$(this.api().table().header()).css({'background-color': '#000', 'color': '#fff'});",
                                              #   "}")
                                            )) # output$gene_summary_table
  ################################################################


  ################ peak_summary table ############################
  peak_summary <- reactive({
    req(reactiveData$peakAnno_df)
    peakAnno_df <- reactiveData$peakAnno_df

    filter_sig_peaks <- input$filter_sig_peaks
    if (filter_sig_peaks == 'YES'){filter = 0.05}
    else { filter = 1}

    row = input$gene_summary_table_cell_clicked
    subset_gene_df <- peakAnno_df[peakAnno_df$SYMBOL==row$value & !is.na(peakAnno_df$SYMBOL) & peakAnno_df$P_Value < filter,
                                               c('logFC', 'P_Value', 'distanceToTSS', 'peak_id', 'annotation')]
    
    # convert values as numeric
    subset_gene_df$logFC <- as.numeric(subset_gene_df$logFC)
    subset_gene_df$P_Value <- as.numeric(subset_gene_df$P_Value)
    subset_gene_df$distanceToTSS <- as.numeric(subset_gene_df$distanceToTSS)
    # round number of decimal
    subset_gene_df$logFC <- round(subset_gene_df$logFC, 5)
    subset_gene_df$P_Value <- round(subset_gene_df$P_Value, 5)

    subset_gene_df
  })
  ################################################################


  ################ peak_summary_description ######################
  output$peak_summary_description <- renderUI({

    gene = input$gene_summary_table_cell_clicked
    print(gene)

    HTML(paste(paste('<font size="+2">Summary table of the peaks from ', '<b>', gene$value, '</b></font>'),
         'Click on one of the peak_id in the table below to view the genome borowser in that region.',
         sep = "<br/>"))
  })
  ################################################################


  ################ peak_summary output table ######################
  output$peak_summary_table <- DT::renderDT(peak_summary(), selection = "single", filter = 'top',
                                            options = list(
                                              autoWidth = TRUE
                                            )) # output$peak_summary_table
  ################################################################


  ################ peak_id from the peak_summary table #############################
  output$peak_id <- renderUI({

    gene = input$gene_summary_table_cell_clicked
    peak_id = input$peak_summary_table_cell_clicked
    URL = paste0("<a href=http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=", peak_id$value,">",
                 paste0('Open the UCSC Gene browser for the gene ', gene$value, ' at the region ', peak_id$value),"</a>")
    HTML(paste(paste('<font size="+2">UCSC Genome Browser on Human (GRCh38/hg38) for ', '<b>', peak_id$value, '</b></font>'),
               URL,
         sep = "<br/>"))
  })
  ################################################################

  
  
  

  
  
} # server


# Create Shiny object
shinyApp(ui = ui, server = server, options = list(width=200, height=200))

# shinyApp(ui = ui,server=server)

