# Rshiny application for exploring various sources of unwanted in TCGA RNA-seq data
## The formatted datasets can be downloaded here

## These function can be found the scripts folder
source('Libraries_HelperFunctions_ForRshiny.R', local = TRUE)

cancer_types <- c(
  "ACC (Adrenocortical carcinoma)",
  "BLCA (Bladder Urothelial Carcinoma)",
  "BRCA (Breast invasive carcinoma)",
  "CESC (Cervical squamous cell carcinoma and endocervical adenocarcinoma)",
  "CHOL (Cholangiocarcinoma)",
  "COAD (Colon adenocarcinoma)",
  "DLBC (Lymphoid Neoplasm Diffuse Large B-cell Lymphoma)",
  "ESCA (Esophageal carcinoma)",
  "GBM (Glioblastoma multiforme)",
  "HNSC (Head and Neck squamous cell carcinoma)",
  "KICH (Kidney Chromophobe)",
  "KIRC (Kidney renal clear cell carcinoma)",
  "KIRP (Kidney renal papillary cell carcinoma)",
  "LAML (Acute Myeloid Leukemia)",
  "LGG (Brain Lower Grade Glioma)",
  "LIHC (Liver hepatocellular carcinoma)",
  "LUAD (Lung adenocarcinoma)",
  "LUSC (Lung squamous cell carcinoma)",
  "MESO (Mesothelioma)",
  "OV (Ovarian serous cystadenocarcinoma)",
  "PAAD (Pancreatic adenocarcinoma)",
  "PCPG (Pheochromocytoma and Paraganglioma)",
  "PRAD (Prostate adenocarcinoma)",
  "READ (Rectum adenocarcinoma)",
  "SARC (Sarcoma)",
  "SKCM (Skin Cutaneous Melanoma)",
  "STAD (Stomach adenocarcinom)",
  "TGCT (Testicular Germ Cell Tumors)",
  "THCA (Thyroid carcinoma)",
  "THYM (Thymoma)",
  "UCEC (Uterine Corpus Endometrial Carcinoma)",
  "UCS (Uterine Carcinosarcoma)",
  "UVM (Uveal Melanoma)"
)

# UI ----------------------------------------------------------------------
myImgResources <- paste0("imgResources/Figure_Rshiny", seq_len(2), ".png")
addResourcePath(prefix = "imgResources", directoryPath = "image")
ui <- navbarPage(
  title = "TCGA RNA-seq unwanted variation", 
  id = "tabs",
  selected = "Launchpad",
  theme = shinythemes::shinytheme("flatly"),
  
  ## TAB 0 ----------------
  tabPanel("Welcome",
           value = "Launchpad",
           fluidPage(
             splitLayout(
               cellWidths = c("50%", "50%"),
               tags$img(src = myImgResources[1], width = "600px", height = "700px"),
               tags$img(src = myImgResources[2], width = "600px", height = "600px")
             )

           )),
  
  ## TAB 1 --------
  tabPanel("Study",
           
           # Application title
           titlePanel("TCGA RNA-seq studies"),
           # Sidebar with a slider input for number of bins
           sidebarLayout(
             sidebarPanel(
               selectInput(
                 inputId = "sample_type",
                 label = "Select cancet type",
                 choices = cancer_types,
                 selected = "READ (Rectum adenocarcinoma)"
               )
             ),
             
             # Show a plot of the generated distribution
             mainPanel(plotOutput("study_design_plot"),
                       downloadButton("study_design_plot_dl", "Download"),
                       verbatimTextOutput("summary_stat"))
           )),
  
  ## TAB 2 --------
  tabPanel("Samples&genes filtering",
           # Sidebar with a slider input for number of bins
           sidebarLayout(
             sidebarPanel(
               h2("Filtering genes and samples"),
               br(),
               selectizeInput(
                 inputId = "biotype",
                 label = "Select Gene Biology type",
                 choices = NULL
               ),
               selectizeInput(
                 inputId = "tissue",
                 label = "Select Tissue",
                 choices = NULL
               ),
               br(),
               sliderInput(
                 inputId = "plate_range",
                 label = "Select Plate range",
                 min = 0,
                 max = 100,
                 value = c(1,20)
               ),
               br(),
               sliderInput(
                 inputId = "min_gene_count",
                 label = "Select Lowest Gene count threshold",
                 min = 0,
                 max = 100,
                 value = 20
               ),
               sliderInput(
                 inputId = "min_sample_count",
                 label = "Select Smallest sample threshold",
                 min = 0,
                 max = 1000,
                 value = 200
               ),
               br(),
               sliderInput(
                 inputId = "lib_size",
                 label = "Select Library size threshold",
                 min = 0,
                 max = 30,
                 value = 24.8
               ),
               br(),
               sliderInput(
                 inputId = "purity",
                 label = "Select purity threshold",
                 min = .5,
                 max = 1,
                 value = .5
               )
             ),
             
             # Show a plot of the generated distribution
             mainPanel(fluidRow(
               splitLayout(
                 cellWidths = c("50%", "50%"),
                 plotOutput("lib_thresh_plot"),
                 plotOutput("lib_thresh_plot_2")
               )
             ),
             plotOutput("purity_plot"),
             verbatimTextOutput("filter_sum"))
           )
  ),
  
  ## TAB 3 --------
  tabPanel("Global level",
           # Sidebar with a slider input for number of bins
           sidebarLayout(
             sidebarPanel(
               p("Press the button below to run PCA on the filtered data."),
               actionButton(inputId = "pca_trigger",
                            label = "Run PCA"),
               p(
                 "Note: please re-run pca after every filterations in the previous step"
               ),
               radioButtons(
                 inputId = "pca_var",
                 label = "Select the variable of interest",
                 choices = c("Tissues", "Year", "Month",
                             "Plates", "TSS", "Center"),
                 selected = "Year"
               ),
               selectInput(
                 inputId = "pca_plot_type",
                 label = "Plot type",
                 choices = c("Scatterplot", "Boxplot"),
                 selected = "Scatterplot"
               ),
               fluidRow(
                 column(width = 3,
                        numericInput("firstPC", "PC #1", 1, 1, 50)),
                 column(width = 3,
                        numericInput("secondPC", "PC #2", 2, 1, 50)),
                 column(width = 3,
                        numericInput("thirdPC", "PC #3", 3, 1, 50))
               ),
               br(),
               br(),
               selectInput(
                 inputId = "pca_corr",
                 label = "Association between PCs (cumulatively) and variables",
                 choices = c("Library_Size", "Purity_singscore", "Year"),
                 selected = character(0)
               ),
               sliderInput(
                 inputId = "nPCs",
                 label = "Number of PCs",
                 min = 2,
                 max = 10,
                 value = 10
               )
             ),
             
             # Show a plot of the generated distribution
             mainPanel(tabsetPanel(
               type = "tabs",
               tabPanel(
                 "PCA plots",
                 plotOutput("pca_plot1"),
                 plotOutput("pca_plot2"),
                 plotOutput("pca_plot3"),
                 plotOutput("pca_plot4")
               ),
               tabPanel("Association between PCs and variable",
                        plotOutput("corr_plot"))
             ))
           )),
  
  ## TAB 4 --------
  tabPanel("Gene level",
           # Sidebar with a slider input for number of bins
           sidebarLayout(
             sidebarPanel(
               selectInput(
                 inputId = "gene_corr",
                 label = "Association between genes and variables",
                 choices = c("Library Size", "Purity",
                             "Year", "Plate"),
                 selected = "Library Size"
               ),
               selectizeInput(
                 inputId = "gene_name",
                 label = "Select gene",
                 selected = "ACTB",
                 choices = NULL
               ),
               checkboxGroupInput(inputId = "gene_dat",
                                  label = "Select counts",
                                  choices = c("Raw", "FPKM",
                                              "FPKM.UQ", "RUV-III"),
                                  selected = "Raw")
             ),
             
             # Show a plot of the generated distribution
             mainPanel(
               plotOutput("gene_corr_plot_raw"),
               plotOutput("gene_corr_plot_fpkm"),
               plotOutput("gene_corr_plot_fpkm_uq"),
               plotOutput("gene_corr_plot_ruv")
             )
           )),
  
  ## TAB 5 --------
  tabPanel("RUV-III normalization",
           tabsetPanel(
             type = "tabs",
             tabPanel(
               "Negative control genes selection",
               sidebarLayout(
                 sidebarPanel(
                   h3("NCG Selection"),
                   selectInput(
                     inputId = "ncg_select",
                     label = "Various sets of housekeeping genes",
                     choices = c('scRNAseq_HK',
                                 'RNAseq_HK',
                                 'Microrray_HK',
                                 'NanostringPanCancer_HK',
                                 'sinscorePanCancer_HK'),
                     selected = "scRNAseq_HK"
                   ),
                   h3("Assessment of selected negative control genes"),
                   p("Press the button below to run PCA on the selected negative control genes"),
                   actionButton(inputId = "pca_trigger2",
                                label = "Run PCA"),
                   radioButtons(
                     inputId = "pca_var2",
                     label = "Select the variable of interest",
                     choices = c( "Tissues", "Year", "Month",
                                  "Plates", "TSS", "Center"),
                     selected = "Year"
                   ),
                   selectInput(
                     inputId = "pca_plot_type2",
                     label = "Plot type",
                     choices = c("Scatterplot", "Boxplot"),
                     selected = "Scatterplot"
                   ),
                   fluidRow(
                     column(width = 3,
                            numericInput("firstPC2", "PC #1", 1, 1, 50)),
                     column(width = 3,
                            numericInput("secondPC2", "PC #2", 2, 1, 50)),
                     column(width = 3,
                            numericInput("thirdPC2", "PC #3", 3, 1, 50))
                   ),
                   br(),
                   br(),
                   selectInput(
                     inputId = "pca_corr2",
                     label = "Select variable for Correlation",
                     choices = c("Library_Size", "Purity_singscore", "Year"),
                     selected = character(0)
                   ),
                   sliderInput(
                     inputId = "nPCs2",
                     label = "Select no. of PCs",
                     min = 2,
                     max = 10,
                     value = 10
                   )
                 ),
                 
                 # Show a plot of the generated distribution
                 mainPanel(tabsetPanel(
                   type = "tabs",
                   tabPanel(
                     "PCA plots",
                     plotOutput("pca_plot1_2")#,
                     # plotOutput("pca_plot2_2"),
                     # plotOutput("pca_plot3_2")
                   ),
                   tabPanel("Association between PCs (cumulatively) and variables",
                            plotOutput("corr_plot_2"))
                 )))
             ),
             tabPanel(
               "PRPS",
               # Sidebar with a slider input for number of bins
               sidebarLayout(
                 sidebarPanel(
                   p("Pseudo-replicates of pseudo-samples"),
                   checkboxGroupInput(inputId = "prps_batch",
                                      label = "Select batch variables",
                                      choices = c("Year", "Plates"),
                                      selected = "Year"),
                   uiOutput("biology_con"),
                   # selectInput(
                   #   inputId = "prps_bio",
                   #   label = "Select biology",
                   #   choices = "biology",
                   #   selected = "biology"
                   # ),
                   checkboxInput(
                     inputId = "include.ls",
                     label = "Include Library Size",
                     value = TRUE
                   ),
                   checkboxInput(
                     inputId = "include.purity",
                     label = "Include Purity",
                     value = TRUE
                   ),
                   numericInput(inputId = "minSamplesPerBatchPS",
                                label = "Number of samples generate PS for batch effects",
                                value = 3,
                                min = 0,
                                max = 10),
                   numericInput(inputId = "minSamplesForPurityPS",
                                label = "Number of samples generate PS for purity effects",
                                value = 3,
                                min = 0,
                                max = 10),
                   numericInput(inputId = "minSamplesForLibrarySizePS",
                                label = "Number of samples generate PS for library size effects",
                                value = 3,
                                min = 0,
                                max = 30),
                   numericInput(inputId = "minSamplesForPurityPerBiology",
                                label = "Minimun number of samples per biology to generate PRPS for purity effects",
                                value = 12,
                                min = 0,
                                max = 15),
                   numericInput(inputId = "minSamplesForLibrarySizePerBatch",
                                label = "Minimun number of samples per plate to generate PRPS for library size effects",
                                value = 10,
                                min = 0,
                                max = 30),
                   actionButton(inputId = "prps_gen",
                                label = "Generate PRPS"),
                   textOutput("prps_gen_conf")
                 ),
                 
                 # Show a plot of the generated distribution
                 mainPanel(plotOutput("PRPS_map"),
                           verbatimTextOutput("prps_summary")
                 )
               )),
             tabPanel(
               "RUV-III-PRPS normalization",
               mainPanel(p("RUV-III-PRPS"),
                         numericInput(
                           inputId = "select_k",
                           label = "Select k for RUV-III",
                           value = 1,
                           min = 0,
                           max = 20),
                         actionButton(inputId = "ruv_gen",
                                      label = "Run RUV-III"),
                         textOutput("ruv_gen_conf")
               ))
           ))
)


# Server ------------------------------------------------------------------
server <- function(input, output, session) {
  
  observe({
    req(input$sample_type)
    data <-
      base::readRDS(paste0(
        input$sample_type,
        ".rds")
        )
    gene.annot <<- as.data.frame(SummarizedExperiment::rowData(data))
    sample.info <<-
      as.data.frame(SummarizedExperiment::colData(data))
    raw.count <<-
      as.data.frame(SummarizedExperiment::assay(data, 'HTseq_counts'))
    Library_Size <<- log2(colSums(raw.count))
    tcga.harmonized <<- names(SummarizedExperiment::assays(data))
    data.set.names <<- tcga.harmonized
    floor_dec <-
      function(x, level = 1)
        round(x - 5 * 10 ^ (-level - 1), level)
    ceiling_dec <-
      function(x, level = 1)
        round(x + 5 * 10 ^ (-level - 1), level)
    
    # Tab 2
    updateSelectizeInput(session,
                         inputId = "biotype",
                         label = "Select gene biotypes",
                         choices = unique(gene.annot$Gene_BioType),
                         selected = "protein.coding",
                         options = list(maxItems = 6))
    updateSelectizeInput(session,
                         inputId = "tissue",
                         label = "Select Tissues",
                         choices = unique(sample.info$Tissues),
                         selected = unique(sample.info$Tissues),
                         options = list(maxItems = 6))
    updateSliderInput(session,
                      inputId = "plate_range",
                      label = "Minimum and maximum number of samples per plate",
                      min = min(dplyr::add_count(sample.info, Plates)$n),
                      max = max(dplyr::add_count(sample.info, Plates)$n),
                      value = c(min(dplyr::add_count(sample.info, Plates)$n), max(dplyr::add_count(sample.info, Plates)$n))
    )
    updateSliderInput(session,
                      inputId = "min_gene_count",
                      label = "Lowest gene count threshold",
                      min = 0,
                      max = 100,
                      value = 20)
    updateSliderInput(session,
                      inputId = "min_sample_count",
                      label = "Smallest sample size threshold",
                      min = 5/100*(dim(data)[[2]]),
                      max = dim(data)[[2]],
                      value = 15/100*(dim(data)[[2]]))
    updateSliderInput(session,
                      inputId = "lib_size",
                      label = "Library size threshold",
                      min = floor(min(Library_Size, na.rm = TRUE)),
                      max = ceiling(max(Library_Size, na.rm = TRUE)),
                      value = floor(min(Library_Size, na.rm = TRUE)))
    updateSliderInput(session,
                      inputId = "purity",
                      label = "Tumour purity threshold",
                      min = floor_dec(min(sample.info$Purity_singscore)),
                      max = ceiling_dec(max(sample.info$Purity_singscore)),
                      value = floor_dec(min(sample.info$Purity_singscore))
    )
    
    # Tab 4
    updateSelectizeInput(session,
                         inputId = "gene_name",
                         label = "Select Gene",
                         choices = unique(gene.annot$hgnc_symbol_BioMart),
                         selected = "ACTB",
                         server = TRUE
    )
  
  })
  
  observe({
    req(input$pca_trigger)
    
    updateActionButton(session, "pca_trigger",
                       label = "Re-run PCA")
  })
  
  observe({
    req(input$pca_trigger2)
    
    updateActionButton(session, "pca_trigger2",
                       label = "Re-run PCA")
  })
  
  # FILTERS ------------------------------------------------------------------
  
  ## TAB 1 ---------
  # Filter Sample Type -- FILTER 1
  selected_sampletype <- reactive({
    data <-
      base::readRDS(
        paste0(
          input$sample_type,
          ".rds"
        )
      )
    values$pca.cancer.tcga.ruv <- list()
    values$final_dat <- SummarizedExperiment()
    data
  })
  
  ## TAB 2 ---------
  # Initial Filters
  init_filters <- reactive({
    data <- selected_sampletype()
    keep.genes <- gene.annot$Gene_BioType %in% input$biotype
    data <- data[keep.genes ,]
    
    sample.info <- as.data.frame(SummarizedExperiment::colData(data))
    keep.samples <- sample.info$Tissues %in% input$tissue
    data <- data[, keep.samples]
    
    sample.info <- as.data.frame(SummarizedExperiment::colData(data))
    keep.samples <- sample.info$Purity_singscore > input$purity
    data <- data[, keep.samples]
    
    sample.info <- as.data.frame(SummarizedExperiment::colData(data))
    keep.samples <- dplyr::add_count(sample.info, Plates)
    keep.samples <- keep.samples$n >= input$plate_range[1] & keep.samples$n <= input$plate_range[2]
    
    raw.count <- as.matrix(SummarizedExperiment::assay(data, 'HTseq_counts'))
    keep.high <- apply(
      raw.count,
      1,
      function(x)
        length(x[x > input$min_gene_count]) > input$min_sample_count)
    data <- data[keep.high, keep.samples]
    data
  })
  
  # Filter using library size
  selected_lib_size <- reactive({
    data <- init_filters()
    raw.count <- as.data.frame(SummarizedExperiment::assay(data, 'HTseq_counts'))
    library_size <- log2(colSums(raw.count))
    keep.samples <- library_size > input$lib_size
    data <- data[, keep.samples]
    
    sample.info <- as.data.frame(SummarizedExperiment::colData(data))
    keep.samples <- dplyr::add_count(sample.info, Plates)
    keep.samples <- keep.samples$n >= input$plate_range[1] & keep.samples$n <= input$plate_range[2]
    #data <- data[, ]
    
    raw.count <- as.matrix(SummarizedExperiment::assay(data, 'HTseq_counts'))
    keep.high <- apply(
      raw.count,
      1,
      function(x)
        length(x[x > input$min_gene_count]) > input$min_sample_count)
    data <- data[keep.high, keep.samples]
    data
  })
  
  output$biology_con <- renderUI({
    data <- selected_lib_size()
    if(length(data$Subtypes)>0){
      selectizeInput(inputId = "prps_bio",
                     label = "Select biology",
                     choices = unique(sample.info$Subtypes),
                     selected = unique(sample.info$Subtypes),
                     options = list(maxItems = 6)
      )
    }
  })
  
  
  ## TAB 3 ---------
  values <- reactiveValues(
    pca.cancer.tcga = list(),
    pca.cancer.tcga2 = list(),
    pca.cancer.tcga.ruv = list(),
    prps = list(),
    ruv_adj = data.frame(),
    final_dat = SummarizedExperiment())
  
  
  ### TAB 3 - PCA ---------------
  observeEvent(input$pca_trigger, {
    # Use filtered dataset
    filter4 <- selected_lib_size()
    #filter4 <- values$final_dat
    
    # PCA
    sample.info <- as.data.frame(SummarizedExperiment::colData(filter4))
    ### we perform pca on raw counts, fpkm and fpkm.uq
    data.set.names <- names(SummarizedExperiment::assays(filter4))
    set.seed(1234)
    pca.cancer.tcga  <- lapply(
      data.set.names,
      function(x) {
        .pca(
          data = as.matrix(SummarizedExperiment::assay(filter4, x)),
          is.log = FALSE)
      })
    names(pca.cancer.tcga) <- data.set.names
    values$pca.cancer.tcga <- pca.cancer.tcga
  })
  
  ### TAB 3 - PCA-RUV ---------------
  observeEvent(input$pca_trigger, {
    # Use filtered dataset
    if(dim(values$final_dat)[1] != 0){
      filter4 <- values$final_dat
      
      # PCA
      sample.info <- as.data.frame(SummarizedExperiment::colData(filter4))
      set.seed(1234)
      ### we perform pca on raw counts, fpkm and fpkm.uq
      #data.set.names <- names(SummarizedExperiment::assays(filter4))
      pca.cancer.tcga.ruv  <-
        .pca(
          data = as.matrix(SummarizedExperiment::assay(filter4, "RUV_III")),
          is.log = TRUE)
      #names(pca.cancer.tcga) <- data.set.names
      values$pca.cancer.tcga.ruv <- pca.cancer.tcga.ruv
    }
  })
  
  
  ### TAB 3 - CORRELATION ------
  selected_corr_var <- reactive({
    # Use filtered dataset
    filter4 <- selected_lib_size()
    #filter4 <- values$final_dat
    sample.info <- as.data.frame(SummarizedExperiment::colData(filter4))
    raw.count <- as.data.frame(SummarizedExperiment::assay(filter4, 'HTseq_counts'))
    sample.info$Library_Size <- log2(colSums(raw.count))
    
    if (is.numeric(sample.info[,input$pca_corr])) {
      corr_var <- sample.info[,input$pca_corr]
    }
    else {
      corr_var <- fastDummies::dummy_cols(sample.info[,input$pca_corr])
      corr_var <- corr_var[, c(2:ncol(corr_var))]
    }
    corr_var
  })
  
  ## TAB 4 ----------
  gene_corr_sp <- reactive({
    filter4 <- selected_lib_size()
    raw.count <-
      as.data.frame(SummarizedExperiment::assay(filter4, 'HTseq_counts'))
    sample.info <-
      as.data.frame(SummarizedExperiment::colData(filter4))
    
    raw.count <- as.matrix(raw.count)
    
    
    if (input$gene_corr == "Library Size") {
      corr_var <- log2(colSums(raw.count))
    }
    else if (input$gene_corr == "Purity") {
      corr_var <- sample.info$Purity_singscore
    }
    else if (input$gene_corr == "Year") {
      corr_var <- sample.info$Year
    }
    else if (input$gene_corr == "Plate") {
      corr_var <- sample.info$Plates
    }
    
    if (is.numeric(corr_var)) {
      plot_hist_raw <- .corr.gene.variable(
        expr.data = raw.count,
        is.log = FALSE,
        variable = corr_var,
        method = 'spearman',
        n.cores = 5
      )}
    
    else{
      plot_hist_raw <- .Ftest(
        data = raw.count,
        variable = corr_var,
        is.log = F,
        n.cores = 5
      )}
    
    plot_hist_raw
    
  })
  
  gene_corr_sp_fpkm <- reactive({
    filter4 <- selected_lib_size()
    raw.count <-
      as.data.frame(SummarizedExperiment::assay(filter4, 'HTseq_counts'))
    fpkm <-
      as.data.frame(SummarizedExperiment::assay(filter4, 'HTseq_FPKM'))
    sample.info <-
      as.data.frame(SummarizedExperiment::colData(filter4))
    
    raw.count <- as.matrix(raw.count)
    fpkm <- as.matrix(fpkm)
    
    if (input$gene_corr == "Library Size") {
      corr_var <- log2(colSums(raw.count))
    }
    else if (input$gene_corr == "Purity") {
      corr_var <- sample.info$Purity_singscore
    }
    else if (input$gene_corr == "Year") {
      corr_var <- sample.info$Year
    }
    else if (input$gene_corr == "Plate") {
      corr_var <- sample.info$Plates
    }
    
    if (is.numeric(corr_var)) {
      
      plot_hist_fpkm <- .corr.gene.variable(
        expr.data = fpkm,
        is.log = FALSE,
        variable = corr_var,
        method = 'spearman',
        n.cores = 5
      )}
    
    else{
      plot_hist_fpkm <- .Ftest(
        data = fpkm,
        variable = corr_var,
        is.log = F,
        n.cores = 5
      )}
    
    plot_hist_fpkm
    
  })
  
  gene_corr_sp_ruv <- reactive({
    filter4 <- values$final_dat
    raw.count <-
      as.data.frame(SummarizedExperiment::assay(filter4, 'HTseq_counts'))
    ruv.iii <-
      as.data.frame(SummarizedExperiment::assay(filter4, 'RUV_III'))
    sample.info <-
      as.data.frame(SummarizedExperiment::colData(filter4))
    
    raw.count <- as.matrix(raw.count)
    ruv.iii <- as.matrix(ruv.iii)
    
    
    if (input$gene_corr == "Library Size") {
      corr_var <- log2(colSums(raw.count))
    }
    else if (input$gene_corr == "Purity") {
      corr_var <- sample.info$Purity_singscore
    }
    else if (input$gene_corr == "Year") {
      corr_var <- sample.info$Year
    }
    else if (input$gene_corr == "Plate") {
      corr_var <- sample.info$Plates
    }
    
    if (is.numeric(corr_var)) {
      plot_hist_ruv <- .corr.gene.variable(
        expr.data = ruv.iii,
        is.log = TRUE,
        variable = corr_var,
        method = 'spearman',
        n.cores = 5
      )}
    
    else{
      plot_hist_ruv <- .Ftest(
        data = ruv.iii,
        variable = corr_var,
        is.log = T,
        n.cores = 5
      )}
    
    plot_hist_ruv
    
  })
  
  gene_corr_sp_fpkm_uq <- reactive({
    filter4 <- selected_lib_size()
    raw.count <-
      as.data.frame(SummarizedExperiment::assay(filter4, 'HTseq_counts'))
    fpkm.uq <-
      as.data.frame(SummarizedExperiment::assay(filter4, 'HTseq_FPKM.UQ'))
    sample.info <-
      as.data.frame(SummarizedExperiment::colData(filter4))
    
    raw.count <- as.matrix(raw.count)
    fpkm.uq <- as.matrix(fpkm.uq)
    
    
    if (input$gene_corr == "Library Size") {
      corr_var <- log2(colSums(raw.count))
    }
    else if (input$gene_corr == "Purity") {
      corr_var <- sample.info$Purity_singscore
    }
    else if (input$gene_corr == "Year") {
      corr_var <- sample.info$Year
    }
    else if (input$gene_corr == "Plate") {
      corr_var <- sample.info$Plates
    }
    
    if (is.numeric(corr_var)) {
      plot_hist_fpkm_uq <- .corr.gene.variable(
        expr.data = fpkm.uq,
        is.log = FALSE,
        variable = corr_var,
        method = 'spearman',
        n.cores = 5
      )}
    
    else{
      plot_hist_fpkm_uq <- .Ftest(
        data = fpkm.uq,
        variable = corr_var,
        is.log = F,
        n.cores = 5
      )}
    
    plot_hist_fpkm_uq
    
  })
  
  ## TAB 5 ----------
  
  selected_ncg <- reactive({
    data <- selected_lib_size()
    gene.annot <-  as.data.frame(SummarizedExperiment::rowData(data))
    keep.genes <- as.vector(gene.annot[input$ncg_select] == "yes")
    data <- data[keep.genes, ]
    data
  })
  
  ### NCG - PCA ------
  observeEvent(input$pca_trigger2, {
    # Use filtered dataset
    data <- selected_ncg()
    
    # PCA
    sample.info <- as.data.frame(SummarizedExperiment::colData(data))
    ### we perform pca on raw counts, fpkm and fpkm.uq
    data.set.names <- names(SummarizedExperiment::assays(data))
    set.seed(1234)
    pca.cancer.tcga  <- lapply(
      data.set.names,
      function(x) {
        .pca(
          data = as.matrix(SummarizedExperiment::assay(data, x)),
          is.log = FALSE)
      })
    names(pca.cancer.tcga) <- data.set.names
    values$pca.cancer.tcga2 <- pca.cancer.tcga
  })
  
  ### NCG - CORRELATION ------
  selected_corr_var2 <- reactive({
    # Use filtered dataset
    data <- selected_ncg()
    sample.info <- as.data.frame(SummarizedExperiment::colData(data))
    raw.count <- as.data.frame(SummarizedExperiment::assay(data, 'HTseq_counts'))
    sample.info$Library_Size <- log2(colSums(raw.count))
    
    if (is.numeric(sample.info[,input$pca_corr2])) {
      corr_var <- sample.info[,input$pca_corr2]
    }
    else {
      corr_var <- fastDummies::dummy_cols(sample.info[,input$pca_corr2])
      corr_var <- corr_var[, c(2:ncol(corr_var))]
    }
    corr_var
  })
  ### RUV-III ------
  observeEvent(input$prps_gen, {
    req(input$sample_type %in% c("BRCA (Breast invasive carcinoma)",
                                 "COAD (Colon adenocarcinoma)",
                                 "LUAD (Lung adenocarcinoma)",
                                 "READ (Rectum adenocarcinoma)"))
    data <- selected_lib_size()
    
    sample.info <- as.data.frame(SummarizedExperiment::colData(data))
    keep.samples <- sample.info$Subtypes %in% input$prps_bio
    data <- data[, keep.samples]
    sample.info <- as.data.frame(SummarizedExperiment::colData(data))
    raw.count <- as.matrix(SummarizedExperiment::assay(data, 'HTseq_counts'))
    sample.info$ls <- log2(colSums(raw.count))
    
    prps <-
      .CreatePseudoSamplesForLsPurityBatch(
        expr.data = raw.count,
        sample.info = sample.info,
        librarySize = "ls",
        batch = input$prps_batch,
        biology = "Subtypes",
        purity = "Purity_singscore",
        include.ls = input$include.ls,
        include.purity = input$include.purity,
        minSamplesPerBatchPS = input$minSamplesPerBatchPS,
        minSamplesForPurityPS = input$minSamplesForPurityPS,
        minSamplesForLibrarySizePS = input$minSamplesForLibrarySizePS,
        minSamplesForPurityPerBiology = input$minSamplesForPurityPerBiology,
        minSamplesForLibrarySizePerBatch = input$minSamplesForLibrarySizePerBatch
      )
    
    values$prps <- prps
    output$prps_gen_conf <- renderText("PRPS Generation Completed.")
    
    
  })
  
  observeEvent(input$ruv_gen, {
    req(input$sample_type %in% c("BRCA (Breast invasive carcinoma)",
                                 "COAD (Colon adenocarcinoma)",
                                 "LUAD (Lung adenocarcinoma)",
                                 "READ (Rectum adenocarcinoma)"))
    data <- selected_lib_size()
    raw.count <- as.matrix(SummarizedExperiment::assay(data, 'HTseq_counts'))
    fpkm <- as.matrix(SummarizedExperiment::assay(data, 'HTseq_FPKM'))
    fpkm.uq <- as.matrix(SummarizedExperiment::assay(data, 'HTseq_FPKM.UQ'))
    gene.annot <- as.data.frame(SummarizedExperiment::rowData(data))
    sample.info <- as.data.frame(SummarizedExperiment::colData(data))
    prps <- values$prps
    
    prps.batch <- prps$ps.batch
    colnames(prps.batch) <- unlist(lapply(
      colnames(prps.batch),
      function(x) strsplit(x, '_')[[1]][1]
    ))
    prps.ls <- prps$ps.ls
    ps.purity <- prps$ps.purity
    ### data input
    if(is.list(prps.ls) & is.list(ps.purity) ){
    ruv.data <- cbind(raw.count, prps.batch)
    ruv.data <- t(log2(ruv.data + 1))
    } else if(!is.list(prps.ls) & !is.list(ps.purity)){
      ruv.data <- cbind(raw.count, prps.batch, prps.ls, ps.purity)
      ruv.data <- t(log2(ruv.data + 1))
    } else if(!is.list(ps.purity)){
      ruv.data <- cbind(raw.count, prps.batch, ps.purity)
      ruv.data <- t(log2(ruv.data + 1))
    } else if(!is.list(prps.ls)){
      ruv.data <- cbind(raw.count, prps.batch, prps.ls)
      ruv.data <- t(log2(ruv.data + 1))
    }
  
    ### replicate matrix
    ruv.rep <- ruv::replicate.matrix(row.names(ruv.data))
    
    ### NCG sets
    message('RUV-III rrr')
    ncg.set <- colnames(ruv.data) %in% gene.annot$hgnc_symbol_BioMart[as.vector(gene.annot[input$ncg_select] == "yes")]
    message(sum(ncg.set))
    message(dim(ruv.data))
    message(sum(ruv.rep))
    #sum(ncg.set)
    message('RUV-III hhhh')
    ruv.iii.nor <-
      .fastRUVIII(
        Y = ruv.data,
        M = ruv.rep,
        ctl = ncg.set,
        k = input$select_k,
        BSPARAM = BiocSingular::bsparam(),
        return.info = TRUE
      )
    
    ruv.iii.adj <- t(ruv.iii.nor$newY[1:ncol(raw.count) , ])
    
    values$final_dat <- SummarizedExperiment(
      assays = list(
        HTseq_counts = raw.count,
        HTseq_FPKM = fpkm,
        HTseq_FPKM.UQ = fpkm.uq,
        RUV_III = ruv.iii.adj
      ),
      colData = sample.info,
      rowData = gene.annot)
    output$ruv_gen_conf <- renderText(
      paste("The RUV-III normalization is completed. k is ", 
            input$select_k)
      )
  
  })
  
  # PLOTS -------------------------------------------------------------------
  ## TAB 1 ---------
  # Study Design plot
  study_design_plot_input <- reactive({
    data <- selected_sampletype()
    data$ls <- log2(colSums(SummarizedExperiment::assay(data, 'HTseq_counts')))
    sample.info <- as.data.frame(SummarizedExperiment::colData(data))
    
    currentCols <-  c(
      RColorBrewer::brewer.pal(8, "Dark2")[-5],
      RColorBrewer::brewer.pal(10, "Paired"),
      RColorBrewer::brewer.pal(12, "Set3"),
      RColorBrewer::brewer.pal(9, "Blues")[c(8, 3, 7, 4, 6, 9, 5)],
      RColorBrewer::brewer.pal(9, "Oranges")[c(8, 3, 7, 4, 6, 9, 5)],
      RColorBrewer::brewer.pal(9, "Greens")[c(8, 3, 7, 4, 6, 9, 5)],
      RColorBrewer::brewer.pal(9, "Purples")[c(8, 3, 7, 4, 6, 9, 5)],
      RColorBrewer::brewer.pal(9, "Reds")[c(8, 3, 7, 4, 6, 9, 5)],
      RColorBrewer::brewer.pal(9, "Greys")[c(8, 3, 7, 4, 6, 9, 5)],
      RColorBrewer::brewer.pal(9, "BuGn")[c(8, 3, 7, 4, 6, 9, 5)],
      RColorBrewer::brewer.pal(9, "PuRd")[c(8, 3, 7, 4, 6, 9, 5)],
      RColorBrewer::brewer.pal(9, "BuPu")[c(8, 3, 7, 4, 6, 9, 5)],
      RColorBrewer::brewer.pal(9, "YlGn")[c(8, 3, 7, 4, 6, 9, 5)],
      RColorBrewer::brewer.pal(10, "Paired")
    )
    
    cols <- c(
      'Year',
      'Plates',
      'TSS',
      'Tissues',
      'Center',
      'ls',
      'Purity_singscore',
      "Tumor.stage"
    )
    
    sample.info <- sample.info[ , cols]
    years.colors <- currentCols[1:length(unique(sample.info$Year))]
    ### Year
    H.time <- ComplexHeatmap::Heatmap(
      rev(sample.info$Year),
      cluster_columns  = FALSE,
      column_names_gp = grid::gpar(fontsize = 18),
      col =  years.colors,
      name = 'Time (years)',
      heatmap_legend_param = list(
        color_bar = "discrete" ,
        ncol = 2,
        title_gp = grid::gpar(fontsize = 18)
      )
    )
    ### Plates
    n.plates <- length(unique(sample.info$Plates)) # 38
    colfunc <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(11, 'PRGn')[-6])
    color.plates <- colfunc(n.plates)
    H.plate <- ComplexHeatmap::Heatmap(
      rev(sample.info$Plates),
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      column_names_gp = grid::gpar(fontsize = 18),
      col = color.plates,
      name = 'Plates',
      heatmap_legend_param = list(
        color_bar = "discrete" ,
        ncol = 4,
        title_gp = grid::gpar(fontsize = 18)
      )
    )
    ### TSS
    n.tss <- length(unique(sample.info$TSS)) # 40
    colfunc <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(11, 'BrBG')[-6])
    color.tss <- colfunc(n.tss)
    H.tss <- ComplexHeatmap::Heatmap(
      rev(sample.info$TSS),
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      column_names_gp = grid::gpar(fontsize = 18),
      col = color.tss,
      name = 'Tissue source sites',
      heatmap_legend_param = list(
        color_bar = "discrete" ,
        ncol = 4,
        title_gp = grid::gpar(fontsize = 18)
      )
    )
    ### Tissue
    # H.tissue <- ComplexHeatmap::Heatmap(
    #   rev(sample.info$Tissues),
    #   cluster_rows = FALSE,
    #   column_names_gp = grid::gpar(fontsize = 18),
    #   col = c("#252525", "#D9D9D9"),
    #   name = 'Tissues',
    #   heatmap_legend_param = list(
    #     color_bar = "discrete" ,
    #     direction = "vertical",
    #     ncol = 1,
    #     title_gp = grid::gpar(fontsize = 18),
    #     labels = unique(sample.info$Tissues)
    #   )
    # )
    ### Purity
    H.purity <- ComplexHeatmap::Heatmap(
      rev(sample.info$Purity_singscore),
      column_names_gp = grid::gpar(fontsize = 18),
      cluster_rows = FALSE,
      name = 'Tumor purity score',
      col = viridis::plasma(n = 10),
      heatmap_legend_param = list(
        title_gp = grid::gpar(fontsize = 18)
      )
    )
    ### library size
    H.ls <- ComplexHeatmap::Heatmap(
      rev(sample.info$ls),
      cluster_rows = FALSE,
      name = 'Library size',
      column_names_gp = grid::gpar(fontsize = 18),
      col = viridis::viridis(n = 10),
      heatmap_legend_param = list(
        title_gp = grid::gpar(fontsize = 18)
      )
    )
    ## Tumor stage
    n.ts <- length(unique(sample.info$Tumor.stage)) # 40
    colfunc <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(11, 'BrBG')[-6])
    color.ts <- colfunc(n.ts)
    H.ts <- ComplexHeatmap::Heatmap(
      rev(sample.info$Tumor.stage),
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      column_names_gp = grid::gpar(fontsize = 18),
      col = color.ts,
      name = 'Tumor stage',
      heatmap_legend_param = list(
        color_bar = "discrete" ,
        ncol = 4,
        title_gp = grid::gpar(fontsize = 18)
      )
    )
    
    ComplexHeatmap::draw(
      H.time +
        H.plate +
        H.tss +
        #H.tissue +
        H.ts +
        H.purity +
        H.ls ,
      merge_legends = FALSE,
      heatmap_legend_side = 'left'
    )
  })
  
  output$study_design_plot <- renderPlot({
    study_design_plot_input()
  })
  
  output$study_design_plot_dl <- downloadHandler(
    filename = function(){
      paste("study_design_plot", "png", sep = ".")
    },
    content = function(file){
      ggsave(study_design_plot_input(), filename = file)
      # png(file = file)
      # study_design_plot_input()
      # dev.off()
    }
  )
  
  output$summary_stat <- renderPrint({
    data <- selected_lib_size()
    # data
    cat('Summary', "\n")
    cat('Number of samples: ', ncol(data), "\n")
    cat('Number of genes: ', nrow(data), "\n")
    cat('Tissue: ', paste0(names(table(data$Tissues)), '(', unname(table(data$Tissues)), '),'), "\n")
    cat('Years: ', paste0(names(table(data$Year)), '(', unname(table(data$Year)), '),'), "\n")
    cat('Plates: ', paste0(names(table(data$Plates)), '(', unname(table(data$Plates)), '),'), "\n")
    cat('Tissue source sites: ', paste0(names(table(data$TSS)), '(', unname(table(data$TSS)), '),'))
  })
  
  ## TAB 2 ---------
  # Library threshold plot
  output$lib_thresh_plot <- renderPlot({
    filter3 <- init_filters()
    raw.count_f3 <- as.data.frame(SummarizedExperiment::assay(filter3, 'HTseq_counts'))
    library_size_f3 <- log2(colSums(raw.count_f3))
    # plot(library_size_f3, xlab = 'sample', ylab = 'log2 library size (total counts)')
    to.plot.ls <- data.frame(ls = library_size_f3, samples = 1:length(library_size_f3))
    ggplot2::ggplot(to.plot.ls, aes(x = samples, y = ls)) +
      geom_point() +
      geom_hline(yintercept = input$lib_size) +
      ylab('Log2 library size (totall counts)') +
      xlab('Samples') +
      ggtitle('Library size') +
      theme(
        panel.background = element_blank(),
        axis.line = element_line(colour = 'black', size = 1),
        plot.title = element_text(size = 22),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 10),
        strip.text = element_text(size = 14))
    # abline(h = input$lib_size)
  })
  
  # Library threshold plot 2
  output$lib_thresh_plot_2 <- renderPlot({
    filter3 <- selected_lib_size()
    raw.count <- as.data.frame(SummarizedExperiment::assay(filter3, 'HTseq_counts'))
    library_size <- log2(colSums(raw.count))
    to.plot.ls.2 <- data.frame(ls = library_size, samples = 1:length(library_size))
    ggplot2::ggplot(to.plot.ls.2, aes(x = samples, y = ls)) +
      geom_point() +
      ylab('Log2 library size (totall counts)') +
      xlab('Samples') +
      ggtitle('Library size (removed samples)') +
      theme(
        panel.background = element_blank(),
        plot.title = element_text(size = 22), 
        axis.line = element_line(colour = 'black', size = 1),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 10),
        strip.text = element_text(size = 14))
  })
  
  # Purity Plot
  output$purity_plot <- renderPlot({
    filter3 <- selected_lib_size()
    sample.info <- as.data.frame(SummarizedExperiment::colData(filter3))
    # plot(sample.info$Purity_singscore)
    to.plot.ls.2 <- data.frame(
      purity = sample.info$Purity_singscore, 
      samples = 1:length(sample.info$Purity_singscore)
    )
    ggplot2::ggplot(to.plot.ls.2, aes(x = samples, y = purity)) +
      geom_point() +
      ylab('Tumour purity scores') +
      xlab('Samples') +
      ggtitle('Tumour purity') +
      theme(
        panel.background = element_blank(),
        plot.title = element_text(size = 22), 
        axis.line = element_line(colour = 'black', size = 1),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 10),
        strip.text = element_text(size = 14))
  })
  
  # Summary
  output$filter_sum <- renderPrint({
    data <- selected_lib_size()
    # cat("    Genes   Samples", "\n")
    # dim(data)
    cat('After above-filterarion', "\n")
    cat('Number of samples: ', ncol(data), "\n")
    cat('Number of genes: ', nrow(data), "\n")
    cat('Tissue: ', paste0(names(table(data$Tissues)), '(', unname(table(data$Tissues)), '),'), "\n")
    cat('Years: ', paste0(names(table(data$Year)), '(', unname(table(data$Year)), '),'), "\n")
    cat('Plates: ', paste0(names(table(data$Plates)), '(', unname(table(data$Plates)), '),'), "\n")
    cat('Tissue source sites: ', paste0(names(table(data$TSS)), '(', unname(table(data$TSS)), '),'))
  })
  
  ## TAB 3 ---------
  output$pca_plot1 <- renderPlot({
    validate(need(
      input$pca_trigger, 
      message = 'Please Start PCA for the plots to appear')
    )
    
    filter4 <- selected_lib_size()
    sample.info <- as.data.frame(SummarizedExperiment::colData(filter4))
    pca.cancer.tcga <- values$pca.cancer.tcga
    
    if (input$pca_plot_type == "Scatterplot") {
      pcs <- pca.cancer.tcga[[tcga.harmonized[1]]]
      p <- .scatter.density.pc(
        pcs = pcs$sing.val$u[, c(input$firstPC, input$secondPC, input$thirdPC)],
        pc.var = pcs$var,
        pcs.no = c(input$firstPC, input$secondPC, input$thirdPC),
        group.name = as.character(input$pca_var),
        group = sample.info[,input$pca_var],
        color = as.factor(unique(sample.info[,input$pca_var])),
        strokeSize = .2,
        pointSize = 3,
        strokeColor = 'gray30',
        alpha = .6
      )
      cowplot::plot_grid(
        plotlist = p, 
        ncol = 4, 
        rel_widths = c(.3, .3, .3, .1), 
        labels = 'Raw counts', 
        label_size = 18
      )
      
    }
    else {
      to.plot.pc <- data.frame(
        PC1 = pca.cancer.tcga$HTseq_counts$sing.val$u[, input$firstPC],
        PC2 = pca.cancer.tcga$HTseq_counts$sing.val$u[, input$secondPC],
        PC3 = pca.cancer.tcga$HTseq_counts$sing.val$u[, input$thirdPC],
        variable = sample.info[,input$pca_var]
      )
      colnames(to.plot.pc)[1] <- paste0('PC', input$firstPC)
      colnames(to.plot.pc)[2] <- paste0('PC', input$secondPC)
      colnames(to.plot.pc)[3] <- paste0('PC', input$thirdPC)
      colnames(to.plot.pc)[4] <- input$pca_var
      to.plot.pc <- to.plot.pc %>% 
        tidyr::pivot_longer(-input$pca_var, names_to = 'pcs', values_to = 'var') %>% 
        data.frame(.)
      ggplot(to.plot.pc, aes(x = to.plot.pc[,1], y = var)) +
        geom_boxplot()+
        facet_wrap(~pcs, scales = 'free') +
        ylab('PC') +
        xlab(input$pca_var) +
        ggtitle('Raw counts') +
        theme(
          panel.background = element_blank(),
          plot.title = element_text(size = 22), 
          axis.line = element_line(colour = 'black', size = 1),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14),
          strip.text.x = element_text(size = 20),
          strip.text = element_text(size = 22))
    }
  })
  
  output$pca_plot2 <- renderPlot({
    req(input$pca_trigger)
    
    filter4 <- selected_lib_size()
    sample.info <- as.data.frame(SummarizedExperiment::colData(filter4))
    pca.cancer.tcga <- values$pca.cancer.tcga
    
    if (input$pca_plot_type == "Scatterplot") {
      pcs <- pca.cancer.tcga[[tcga.harmonized[2]]]
      p <- .scatter.density.pc(
        pcs = pcs$sing.val$u[, c(input$firstPC, input$secondPC, input$thirdPC)],
        pc.var = pcs$var,
        pcs.no = c(input$firstPC, input$secondPC, input$thirdPC),
        group.name = as.character(input$pca_var),
        group = sample.info[,input$pca_var],
        color = as.factor(unique(sample.info[,input$pca_var])),
        strokeSize = .2,
        pointSize = 3,
        strokeColor = 'gray30',
        alpha = .6
      )
      cowplot::plot_grid(
        plotlist = p, 
        ncol = 4, 
        rel_widths = c(.3, .3, .3, .1), 
        labels = 'FPKM',
        label_size = 18
      )
    }
    else {
      to.plot.pc <- data.frame(
        PC1 = pca.cancer.tcga$HTseq_FPKM$sing.val$u[, input$firstPC],
        PC2 = pca.cancer.tcga$HTseq_FPKM$sing.val$u[, input$secondPC],
        PC3 = pca.cancer.tcga$HTseq_FPKM$sing.val$u[, input$thirdPC],
        variable = sample.info[,input$pca_var]
      )
      colnames(to.plot.pc)[1] <- paste0('PC', input$firstPC)
      colnames(to.plot.pc)[2] <- paste0('PC', input$secondPC)
      colnames(to.plot.pc)[3] <- paste0('PC', input$thirdPC)
      colnames(to.plot.pc)[4] <- input$pca_var
      to.plot.pc <- to.plot.pc %>% 
        tidyr::pivot_longer(-input$pca_var, names_to = 'pcs', values_to = 'var') %>% 
        data.frame(.)
      ggplot(to.plot.pc, aes(x = to.plot.pc[,1], y = var)) +
        geom_boxplot()+
        facet_wrap(~pcs, scales = 'free') +
        ylab('PC') +
        xlab(input$pca_var) +
        ggtitle('FPKM') +
        theme(
          panel.background = element_blank(),
          plot.title = element_text(size = 22), 
          axis.line = element_line(colour = 'black', size = 1),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14),
          strip.text.x = element_text(size = 20),
          strip.text = element_text(size = 22))
    }
  })
  
  output$pca_plot3 <- renderPlot({
    req(input$pca_trigger)
    
    filter4 <- selected_lib_size()
    sample.info <-
      as.data.frame(SummarizedExperiment::colData(filter4))
    pca.cancer.tcga <- values$pca.cancer.tcga
    
    if (input$pca_plot_type == "Scatterplot") {
      pcs <- pca.cancer.tcga[[tcga.harmonized[3]]]
      p <- .scatter.density.pc(
        pcs = pcs$sing.val$u[, c(input$firstPC, input$secondPC, input$thirdPC)],
        pc.var = pcs$var,
        pcs.no = c(input$firstPC, input$secondPC, input$thirdPC),
        group.name = as.character(input$pca_var),
        group = sample.info[,input$pca_var],
        color = as.factor(unique(sample.info[,input$pca_var])),
        strokeSize = .2,
        pointSize = 3,
        strokeColor = 'gray30',
        alpha = .6
      )
      cowplot::plot_grid(
        plotlist = p, 
        ncol = 4, 
        rel_widths = c(.3, .3, .3, .1), 
        labels = 'FPKM.UQ',
        label_size = 18
      )
    }
    else {
      to.plot.pc <- data.frame(
        PC1 = pca.cancer.tcga$HTseq_FPKM.UQ$sing.val$u[, input$firstPC],
        PC2 = pca.cancer.tcga$HTseq_FPKM.UQ$sing.val$u[, input$secondPC],
        PC3 = pca.cancer.tcga$HTseq_FPKM.UQ$sing.val$u[, input$thirdPC],
        variable = sample.info[,input$pca_var]
      )
      colnames(to.plot.pc)[1] <- paste0('PC', input$firstPC)
      colnames(to.plot.pc)[2] <- paste0('PC', input$secondPC)
      colnames(to.plot.pc)[3] <- paste0('PC', input$thirdPC)
      colnames(to.plot.pc)[4] <- input$pca_var
      to.plot.pc <- to.plot.pc %>% 
        tidyr::pivot_longer(-input$pca_var, names_to = 'pcs', values_to = 'var') %>% 
        data.frame(.)
      ggplot(to.plot.pc, aes(x = to.plot.pc[,1], y = var)) +
        geom_boxplot()+
        facet_wrap(~pcs, scales = 'free') +
        ylab('PC') +
        xlab(input$pca_var)+
        ggtitle('FPKM-UQ') +
        theme(
          panel.background = element_blank(),
          plot.title = element_text(size = 22), 
          axis.line = element_line(colour = 'black', size = 1),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14),
          strip.text.x = element_text(size = 20),
          strip.text = element_text(size = 22))
    }
  })
  
  output$pca_plot4 <- renderPlot({
    req(input$pca_trigger)
    req(input$ruv_gen)
    
    filter4 <- values$final_dat
    sample.info <- as.data.frame(SummarizedExperiment::colData(filter4))
    pca.cancer.tcga.ruv <- values$pca.cancer.tcga.ruv
    
    if (input$pca_plot_type == "Scatterplot") {
      pcs <- pca.cancer.tcga.ruv
      p <- .scatter.density.pc(
        pcs = pcs$sing.val$u[, c(input$firstPC, input$secondPC, input$thirdPC)],
        pc.var = pcs$var,
        pcs.no = c(input$firstPC, input$secondPC, input$thirdPC),
        group.name = as.character(input$pca_var),
        group = sample.info[,input$pca_var],
        color = as.factor(unique(sample.info[,input$pca_var])),
        strokeSize = .2,
        pointSize = 3,
        strokeColor = 'gray30',
        alpha = .6
      )
      cowplot::plot_grid(
        plotlist = p, 
        ncol = 4, 
        rel_widths = c(.3, .3, .3, .1), 
        labels = 'RUV-III',
        label_size = 18
      )
    }
    else {
      par(mfrow = c(1, 3))
      to.plot.pc <- data.frame(
        PC1 = pca.cancer.tcga.ruv$sing.val$u[, input$firstPC],
        PC2 = pca.cancer.tcga.ruv$sing.val$u[, input$secondPC],
        PC3 = pca.cancer.tcga.ruv$sing.val$u[, input$thirdPC],
        variable = sample.info[,input$pca_var]
      )
      colnames(to.plot.pc)[1] <- paste0('PC', input$firstPC)
      colnames(to.plot.pc)[2] <- paste0('PC', input$secondPC)
      colnames(to.plot.pc)[3] <- paste0('PC', input$thirdPC)
      colnames(to.plot.pc)[4] <- input$pca_var
      to.plot.pc <- to.plot.pc %>% 
        tidyr::pivot_longer(-input$pca_var, names_to = 'pcs', values_to = 'var') %>% 
        data.frame(.)
      ggplot(to.plot.pc, aes(x = to.plot.pc[,1], y = var)) +
        geom_boxplot()+
        facet_wrap(~pcs, scales = 'free') +
        ylab('PC') +
        xlab(input$pca_var) +
        ggtitle('RUV-III') +
        theme(
          panel.background = element_blank(),
          plot.title = element_text(size = 22), 
          axis.line = element_line(colour = 'black', size = 1),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14),
          strip.text.x = element_text(size = 20),
          strip.text = element_text(size = 22))
      
    }
  })
  
  output$corr_plot <- renderPlot({
    req(input$pca_corr)
    # Color Pallette Set
    dataSets.colors <- wesanderson::wes_palette(n = 4,
                                                name = "GrandBudapest1")[c(1, 2, 4, 3)]
    names(dataSets.colors) <- c('Raw counts',
                                'FPKM',
                                'FPKM.UQ',
                                'RUV.III')
    nPCs <- input$nPCs
    pca.cancer.tcga <- values$pca.cancer.tcga
    corr_var <- selected_corr_var()
    final_dat <- values$final_dat
    
    if (is.numeric(corr_var)) {
      if (dim(final_dat)[1] != 0) {
        pca.cancer.tcga$RUV_III <- values$pca.cancer.tcga.ruv
        tcga.harmonized <- names(pca.cancer.tcga)
        lreg.ls.cancer.tcga <- lapply(tcga.harmonized,
                                      function(x) {
                                        pcs <- pca.cancer.tcga[[x]]$sing.val$u
                                        tcga.ls.rSquared <- sapply(1:nPCs,
                                                                   function(y) {
                                                                     lm.ls <- summary(lm(corr_var ~ pcs[, 1:y]))$r.squared
                                                                   })
                                      })
        names(lreg.ls.cancer.tcga) <- tcga.harmonized
        
        # Create Dataframe for corrplot
        ls.lnreg.normAssess <- data.frame(
          Raw.counts = lreg.ls.cancer.tcga$HTseq_counts,
          FPKM = lreg.ls.cancer.tcga$HTseq_FPKM,
          FPKM.UQ = lreg.ls.cancer.tcga$HTseq_FPKM.UQ,
          RUV.III = lreg.ls.cancer.tcga$RUV_III,
          pcs = c(1:nPCs)
        )
        
        
        # Modify DF for corrplot
        ls.lnreg.normAssess <- ls.lnreg.normAssess %>%
          tidyr::pivot_longer(-pcs,
                              names_to = 'Datasets',
                              values_to = 'ls') %>%
          dplyr::mutate(Datasets = replace(Datasets,
                                           Datasets == 'Raw.counts', 'Raw counts')) %>%
          dplyr::mutate(Datasets = factor(
            x = Datasets,
            levels = c('Raw counts', 'FPKM', 'FPKM.UQ', 'RUV.III')
          )) %>%
          data.frame(.)
        
        # Visualise
        ggplot(ls.lnreg.normAssess, aes(x = pcs, y = ls, group = Datasets)) +
          geom_line(aes(color = Datasets), size = 1) +
          geom_point(aes(color = Datasets), size = 3) +
          xlab('PCs') + ylab (expression("R" ^ "2")) +
          scale_color_manual(
            values = c(dataSets.colors[1:4]),
            labels = c('Raw counts', 'FPKM', 'FPKM.UQ', 'RUV.III')
          ) +
          scale_x_continuous(breaks = (1:nPCs),
                             labels = c('PC1', paste0('PC1:', 2:nPCs))) +
          scale_y_continuous(breaks = scales::pretty_breaks(n = 5),
                             limits = c(0, 1)) +
          theme(
            panel.background = element_blank(),
            axis.line = element_line(colour = 'black', size = 1),
            axis.title.x = element_text(size = 14),
            axis.title.y = element_text(size = 14),
            axis.text.x = element_text(
              size = 12,
              angle = 45,
              hjust = 1
            ),
            axis.text.y = element_text(size = 12),
            legend.text = element_text(size = 10),
            legend.title = element_text(size = 14),
            strip.text.x = element_text(size = 10)
          ) +
          ggtitle(input$pca_corr)
      }
      else{
        lreg.ls.cancer.tcga <- lapply(tcga.harmonized,
                                      function(x) {
                                        pcs <- pca.cancer.tcga[[x]]$sing.val$u
                                        tcga.ls.rSquared <- sapply(1:nPCs,
                                                                   function(y) {
                                                                     lm.ls <- summary(lm(corr_var ~ pcs[, 1:y]))$r.squared
                                                                   })
                                      })
        names(lreg.ls.cancer.tcga) <- tcga.harmonized
        
        # Create Dataframe for corrplot
        ls.lnreg.normAssess <- data.frame(
          Raw.counts = lreg.ls.cancer.tcga$HTseq_counts,
          FPKM = lreg.ls.cancer.tcga$HTseq_FPKM,
          FPKM.UQ = lreg.ls.cancer.tcga$HTseq_FPKM.UQ,
          pcs = c(1:nPCs)
        )
        
        
        # Modify DF for corrplot
        ls.lnreg.normAssess <- ls.lnreg.normAssess %>%
          tidyr::pivot_longer(-pcs,
                              names_to = 'Datasets',
                              values_to = 'ls') %>%
          dplyr::mutate(Datasets = replace(Datasets,
                                           Datasets == 'Raw.counts', 'Raw counts')) %>%
          dplyr::mutate(Datasets = factor(
            x = Datasets,
            levels = c('Raw counts', 'FPKM', 'FPKM.UQ')
          )) %>%
          data.frame(.)
        
        # Visualise
        ggplot(ls.lnreg.normAssess, aes(x = pcs, y = ls, group = Datasets)) +
          geom_line(aes(color = Datasets), size = 1) +
          geom_point(aes(color = Datasets), size = 3) +
          xlab('PCs') + ylab (expression("R" ^ "2")) +
          scale_color_manual(
            values = c(dataSets.colors[1:3]),
            labels = c('Raw counts', 'FPKM', 'FPKM.UQ')
          ) +
          scale_x_continuous(breaks = (1:nPCs),
                             labels = c('PC1', paste0('PC1:', 2:nPCs))) +
          scale_y_continuous(breaks = scales::pretty_breaks(n = 5),
                             limits = c(0, 1)) +
          theme(
            panel.background = element_blank(),
            axis.line = element_line(colour = 'black', size = 1),
            axis.title.x = element_text(size = 14),
            axis.title.y = element_text(size = 14),
            axis.text.x = element_text(
              size = 12,
              angle = 45,
              hjust = 1
            ),
            axis.text.y = element_text(size = 12),
            legend.text = element_text(size = 10),
            legend.title = element_text(size = 14),
            strip.text.x = element_text(size = 10)
          ) +
          ggtitle(input$pca_corr)
      }
    }
    
    else {
      if(dim(final_dat)[1] != 0){
        pca.cancer.tcga$RUV_III <- values$pca.cancer.tcga.ruv
        tcga.harmonized <- names(pca.cancer.tcga)
        cca.time.year <-
          lapply(tcga.harmonized,
                 function(x) {
                   sapply(1:nPCs,
                          function(y) {
                            cca.pam50 <- stats::cancor(x = pca.cancer.tcga[[x]]$sing.val$u[, 1:y, drop = FALSE],
                                                       y = corr_var)
                            1 - prod(1 - cca.pam50$cor ^ 2)
                          })
                   
                 })
        names(cca.time.year) <- tcga.harmonized
        
        # Create DF for corrplot
        cca.time.year <- data.frame(
          Raw.counts = cca.time.year$HTseq_counts,
          FPKM = cca.time.year$HTseq_FPKM,
          FPKM.UQ = cca.time.year$HTseq_FPKM.UQ,
          RUV.III = cca.time.year$RUV_III,
          pcs = c(1:nPCs)
        )
        
        # Modify DF for corrplot
        cca.time.year <- cca.time.year %>%
          tidyr::pivot_longer(
            -pcs,
            names_to = 'Datasets',
            values_to = 'ls') %>% 
          dplyr::mutate(Datasets = replace(
            Datasets,
            Datasets == 'Raw.counts', 'Raw counts')) %>%
          
          dplyr::mutate(Datasets = factor(
            x = Datasets,
            levels = c('Raw counts', 'FPKM', 'FPKM.UQ', 'RUV.III')
          )) %>%
          data.frame(.)
        
        # Visualise
        ggplot(cca.time.year, aes(x = pcs, y = ls, group = Datasets)) +
          geom_line(aes(color = Datasets), size = 1) +
          geom_point(aes(color = Datasets), size = 3) +
          xlab('PCs') + ylab (expression("Vector Correlation")) +
          scale_color_manual(
            values = c(dataSets.colors),
            labels = c('Raw counts', 'FPKM', 'FPKM.UQ', 'RUV.III')
          ) +
          scale_x_continuous(
            breaks = (1:nPCs),
            labels = c('PC1', paste0('PC1:', 2:nPCs))) +
          scale_y_continuous(
            breaks = scales::pretty_breaks(n = 5),
            limits = c(0, 1)) +
          theme(
            panel.background = element_blank(),
            axis.line = element_line(colour = 'black', size = 1),
            axis.title.x = element_text(size = 14),
            axis.title.y = element_text(size = 14),
            axis.text.x = element_text(
              size = 12,
              angle = 45,
              hjust = 1
            ),
            axis.text.y = element_text(size = 12),
            legend.text = element_text(size = 10),
            legend.title = element_text(size = 14),
            strip.text.x = element_text(size = 10)
          ) +
          ggtitle(input$pca_corr)
      }
      else{
        cca.time.year <-
          lapply(data.set.names,
                 function(x) {
                   sapply(1:nPCs,
                          function(y) {
                            cca.pam50 <- stats::cancor(x = pca.cancer.tcga[[x]]$sing.val$u[, 1:y, drop = FALSE],
                                                       y = corr_var)
                            1 - prod(1 - cca.pam50$cor ^ 2)
                          })
                   
                 })
        names(cca.time.year) <- data.set.names
        
        # Create DF for corrplot
        cca.time.year <- data.frame(
          Raw.counts = cca.time.year$HTseq_counts,
          FPKM = cca.time.year$HTseq_FPKM,
          FPKM.UQ = cca.time.year$HTseq_FPKM.UQ,
          pcs = c(1:nPCs)
        )
        dataSets.colors <- wesanderson::wes_palette(
          n = 4,
          name = "GrandBudapest1")[c(1, 2, 4)]
        names(dataSets.colors) <- c(
          'Raw counts',
          'FPKM',
          'FPKM.UQ')
        
        # Modify DF for corrplot
        cca.time.year <- cca.time.year %>%
          tidyr::pivot_longer(
            -pcs,
            names_to = 'Datasets',
            values_to = 'ls') %>% 
          dplyr::mutate(Datasets = replace(
            Datasets,
            Datasets == 'Raw.counts', 'Raw counts')) %>%
          
          dplyr::mutate(Datasets = factor(
            x = Datasets,
            levels = c('Raw counts', 'FPKM', 'FPKM.UQ')
          )) %>%
          data.frame(.)
        
        # Visualise
        ggplot(cca.time.year, aes(x = pcs, y = ls, group = Datasets)) +
          geom_line(aes(color = Datasets), size = 1) +
          geom_point(aes(color = Datasets), size = 3) +
          xlab('PCs') + ylab (expression("Vector Correlation")) +
          scale_color_manual(
            values = c(dataSets.colors),
            labels = c('Raw counts', 'FPKM', 'FPKM.UQ')
          ) +
          scale_x_continuous(
            breaks = (1:nPCs),
            labels = c('PC1', paste0('PC1:', 2:nPCs))) +
          scale_y_continuous(
            breaks = scales::pretty_breaks(n = 5),
            limits = c(0, 1)) +
          theme(
            panel.background = element_blank(),
            axis.line = element_line(colour = 'black', size = 1),
            axis.title.x = element_text(size = 14),
            axis.title.y = element_text(size = 14),
            axis.text.x = element_text(
              size = 12,
              angle = 45,
              hjust = 1
            ),
            axis.text.y = element_text(size = 12),
            legend.text = element_text(size = 10),
            legend.title = element_text(size = 14),
            strip.text.x = element_text(size = 10)
          ) +
          ggtitle(input$pca_corr)
      }}
  })
  
  
  
  ## TAB 4 ---------
  output$gene_corr_plot_raw <- renderPlot({
    req("Raw" %in% input$gene_dat)
    filter4 <- selected_lib_size()
    raw.count <-
      as.data.frame(SummarizedExperiment::assay(filter4, 'HTseq_counts'))
    sample.info <-
      as.data.frame(SummarizedExperiment::colData(filter4))
    
    raw.count <- as.matrix(raw.count)
    
    if (input$gene_corr == "Library Size") {
      corr_var <- log2(colSums(raw.count))
    }
    else if (input$gene_corr == "Purity") {
      corr_var <- sample.info$Purity_singscore
    }
    else if (input$gene_corr == "Year") {
      corr_var <- sample.info$Year
    }
    else if (input$gene_corr == "Plate") {
      corr_var <- sample.info$Plates
    }
    
    if (is.numeric(corr_var)){
      genes.ls <- gene_corr_sp()
      to.plot.corr <- data.frame(corr = genes.ls$rho)
      p1 <- ggplot2::ggplot(to.plot.corr, aes(x = corr)) +
        geom_histogram() +
        xlab('Spearman correlation') +
        ylab("Counts") +
        ggtitle('Raw counts') +
        xlim(c(-1, 1)) +
        theme(
          panel.background = element_blank(),
          plot.title = element_text(size = 24),
          axis.line = element_line(colour = 'black', size = 1),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 14),
          strip.text.x = element_text(size = 10)
        )
      to.plot.gene <- data.frame(var = corr_var, gene = log2(raw.count[input$gene_name ,] + 1))
      p2 <- ggplot2::ggplot(to.plot.gene, aes(x = var, y = gene)) +
        geom_point() +
        ylab(input$gene_name) +
        xlab(input$gene_corr) +
        ggpubr::stat_cor(
          aes(label = ..r.label..),
          label.x.npc = .4,
          label.y.npc = 1,
          hjust = 0,
          col = 'red', 
          r.accuracy = .1,
          size = 10,
          cor.coef.name = "rho") +
        geom_smooth(
          method = 'lm',
          formula = y ~ x,
          col = 'red',
          se = T ) +
        theme(
          panel.background = element_blank(),
          axis.line = element_line(colour = 'black', size = 1),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 14),
          strip.text.x = element_text(size = 10)
        )
      gridExtra::grid.arrange(p1, p2, ncol = 2)
      
      
    }
    else {
      time.effects <- gene_corr_sp()
      par(mfrow = c(1, 2))
      hist(log2(time.effects$FValue))
      boxplot(log2(raw.count[input$gene_name ,] + 1) ~ corr_var)
      
    }
    
  })
  
  output$gene_corr_plot_fpkm <- renderPlot({
    req("FPKM" %in% input$gene_dat)
    filter4 <- selected_lib_size()
    raw.count <-
      as.data.frame(SummarizedExperiment::assay(filter4, 'HTseq_counts'))
    fpkm <-
      as.data.frame(SummarizedExperiment::assay(filter4, 'HTseq_FPKM'))
    sample.info <-
      as.data.frame(SummarizedExperiment::colData(filter4))
    
    raw.count <- as.matrix(raw.count)
    fpkm <- as.matrix(fpkm)
    
    if (input$gene_corr == "Library Size") {
      corr_var <- log2(colSums(raw.count))
    }
    else if (input$gene_corr == "Purity") {
      corr_var <- sample.info$Purity_singscore
    }
    else if (input$gene_corr == "Year") {
      corr_var <- sample.info$Year
    }
    else if (input$gene_corr == "Plate") {
      corr_var <- sample.info$Plates
    }
    
    if (is.numeric(corr_var)){
      genes.ls <- gene_corr_sp_fpkm()
      to.plot.corr <- data.frame(corr = genes.ls$rho)
      p1 <- ggplot2::ggplot(to.plot.corr, aes(x = corr)) +
        geom_histogram() +
        xlab('Spearman correlation') +
        ylab("Counts") +
        ggtitle('FPKM') +
        xlim(c(-1, 1)) +
        theme(
          panel.background = element_blank(),
          axis.line = element_line(colour = 'black', size = 1),
          plot.title = element_text(size = 24),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 14),
          strip.text.x = element_text(size = 10)
        )
      to.plot.gene <- data.frame(var = corr_var, gene = log2(fpkm[input$gene_name ,] + 1))
      p2 <- ggplot2::ggplot(to.plot.gene, aes(x = var, y = gene)) +
        geom_point() +
        ylab(input$gene_name) +
        xlab(input$gene_corr) +
        ggpubr::stat_cor(
          aes(label = ..r.label..),
          label.x.npc = .4,
          label.y.npc = 1,
          hjust = 0,
          col = 'red', 
          r.accuracy = .1,
          size = 10,
          cor.coef.name = "rho") +
        geom_smooth(
          method = 'lm',
          formula = y ~ x,
          col = 'red',
          se = T ) +
        theme(
          panel.background = element_blank(),
          axis.line = element_line(colour = 'black', size = 1),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 14),
          strip.text.x = element_text(size = 10)
        )
      gridExtra::grid.arrange(p1, p2, ncol = 2)
      
    }
    else {
      time.effects <- gene_corr_sp_fpkm()
      par(mfrow = c(1, 2))
      hist(log2(time.effects$FValue))
      boxplot(log2(fpkm[input$gene_name ,] + 1) ~ corr_var)
      
    }
    
  })
  
  output$gene_corr_plot_fpkm_uq <- renderPlot({
    req("FPKM.UQ" %in% input$gene_dat)
    filter4 <- selected_lib_size()
    raw.count <-
      as.data.frame(SummarizedExperiment::assay(filter4, 'HTseq_counts'))
    fpkm.uq <-
      as.data.frame(SummarizedExperiment::assay(filter4, 'HTseq_FPKM.UQ'))
    sample.info <-
      as.data.frame(SummarizedExperiment::colData(filter4))
    
    raw.count <- as.matrix(raw.count)
    fpkm.uq <- as.matrix(fpkm.uq)
    
    if (input$gene_corr == "Library Size") {
      corr_var <- log2(colSums(raw.count))
    }
    else if (input$gene_corr == "Purity") {
      corr_var <- sample.info$Purity_singscore
    }
    else if (input$gene_corr == "Year") {
      corr_var <- sample.info$Year
    }
    else if (input$gene_corr == "Plate") {
      corr_var <- sample.info$Plates
    }
    
    if (is.numeric(corr_var)){
      genes.ls <- gene_corr_sp_fpkm_uq()
      to.plot.corr <- data.frame(corr = genes.ls$rho)
      p1 <- ggplot2::ggplot(to.plot.corr, aes(x = corr)) +
        geom_histogram() +
        xlab('Spearman correlation') +
        ylab("Counts") +
        ggtitle('FPKM.UQ') +
        xlim(c(-1, 1)) +
        theme(
          panel.background = element_blank(),
          axis.line = element_line(colour = 'black', size = 1),
          plot.title = element_text(size = 24),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 14),
          strip.text.x = element_text(size = 10)
        )
      to.plot.gene <- data.frame(var = corr_var, gene = log2(fpkm.uq[input$gene_name ,] + 1))
      p2 <- ggplot2::ggplot(to.plot.gene, aes(x = var, y = gene)) +
        geom_point() +
        ylab(input$gene_name) +
        xlab(input$gene_corr) +
        ggpubr::stat_cor(
          aes(label = ..r.label..),
          label.x.npc = .4,
          label.y.npc = 1,
          hjust = 0,
          col = 'red', 
          r.accuracy = .1,
          size = 10,
          cor.coef.name = "rho") +
        geom_smooth(
          method = 'lm',
          formula = y ~ x,
          col = 'red',
          se = T ) +
        theme(
          panel.background = element_blank(),
          axis.line = element_line(colour = 'black', size = 1),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 14),
          strip.text.x = element_text(size = 10)
        )
      gridExtra::grid.arrange(p1, p2, ncol = 2, top = "RUV-III")
      
      
      
    }
    else {
      time.effects <- gene_corr_sp_fpkm_uq()
      par(mfrow = c(1, 2))
      hist(log2(time.effects$FValue))
      boxplot(log2(fpkm.uq[input$gene_name ,] + 1) ~ corr_var)
      
    }
  })
  
  output$gene_corr_plot_ruv <- renderPlot({
    req(input$ruv_gen)
    req("RUV-III" %in% input$gene_dat)
    filter4 <- values$final_dat
    raw.count <-
      as.data.frame(SummarizedExperiment::assay(filter4, 'HTseq_counts'))
    ruv.iii <-
      as.data.frame(SummarizedExperiment::assay(filter4, 'RUV_III'))
    sample.info <-
      as.data.frame(SummarizedExperiment::colData(filter4))
    
    raw.count <- as.matrix(raw.count)
    ruv.iii <- as.matrix(ruv.iii)
    
    if (input$gene_corr == "Library Size") {
      corr_var <- log2(colSums(raw.count))
    }
    else if (input$gene_corr == "Purity") {
      corr_var <- sample.info$Purity_singscore
    }
    else if (input$gene_corr == "Year") {
      corr_var <- sample.info$Year
    }
    else if (input$gene_corr == "Plate") {
      corr_var <- sample.info$Plates
    }
    
    if (is.numeric(corr_var)){
      genes.ls <- gene_corr_sp_ruv()
      to.plot.corr <- data.frame(corr = genes.ls$rho)
      p1 <- ggplot2::ggplot(to.plot.corr, aes(x = corr)) +
        geom_histogram() +
        xlab('Spearman correlation') +
        ylab("RUV-III") +
        xlim(c(-1, 1)) +
        theme(
          panel.background = element_blank(),
          axis.line = element_line(colour = 'black', size = 1),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 14),
          strip.text.x = element_text(size = 10)
        )
      to.plot.gene <- data.frame(var = corr_var, gene = ruv.iii[input$gene_name ,])
      p2 <- ggplot2::ggplot(to.plot.gene, aes(x = var, y = gene)) +
        geom_point() +
        ylab(input$gene_name) +
        xlab(input$gene_corr) +
        theme(
          panel.background = element_blank(),
          axis.line = element_line(colour = 'black', size = 1),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 14),
          strip.text.x = element_text(size = 10)
        )
      gridExtra::grid.arrange(p1, p2, ncol = 2)
      
    }
    else {
      time.effects <- gene_corr_sp_ruv()
      par(mfrow = c(1, 2))
      hist(log2(time.effects$FValue))
      boxplot(ruv.iii[input$gene_name ,] ~ corr_var)
      
    }
    
  })
  
  ## TAB 5 RUV-III---------
  ## TAB 5 - NCG - PCA ------
  output$pca_plot1_2 <- renderPlot({
    validate(need(input$pca_trigger2, message = 'Please Start PCA for the plots to appear'))
    
    filter4 <- selected_ncg()
    sample.info <- as.data.frame(SummarizedExperiment::colData(filter4))
    pca.cancer.tcga <- values$pca.cancer.tcga2
    
    if (input$pca_plot_type2 == "Scatterplot") {
      pcs <- pca.cancer.tcga[[tcga.harmonized[1]]]
      p <- .scatter.density.pc(
        pcs = pcs$sing.val$u[, c(input$firstPC2, input$secondPC2, input$thirdPC2)],
        pc.var = pcs$variation,
        pcs.no =  c(input$firstPC2, input$secondPC2, input$thirdPC2), 
        group.name = as.character(input$pca_var2),
        group = sample.info[,input$pca_var2],
        color = as.factor(unique(sample.info[,input$pca_var2])),
        strokeSize = .2,
        pointSize = 3,
        strokeColor = 'gray30',
        alpha = .6
      )
      # do.call(gridExtra::grid.arrange,
      #         c(p,
      #           ncol = 4, top = tcga.harmonized[1]))
      cowplot::plot_grid(
        plotlist = p, 
        ncol = 4, 
        rel_widths = c(.3, .3, .3, .1), 
        labels = 'Raw counts',
        label_size = 18
      )
    }
    else {
      to.plot.pc <- data.frame(
        PC1 = pca.cancer.tcga$HTseq_counts$sing.val$u[, input$firstPC2],
        PC2 = pca.cancer.tcga$HTseq_counts$sing.val$u[, input$secondPC2],
        PC3 = pca.cancer.tcga$HTseq_counts$sing.val$u[, input$thirdPC2],
        variable = sample.info[,input$pca_var2]
      )
      colnames(to.plot.pc)[1] <- paste0('PC', input$firstPC2)
      colnames(to.plot.pc)[2] <- paste0('PC', input$secondPC2)
      colnames(to.plot.pc)[3] <- paste0('PC', input$thirdPC2)
      colnames(to.plot.pc)[4] <- input$pca_var2
      to.plot.pc <- to.plot.pc %>% 
        tidyr::pivot_longer(-input$pca_var2, names_to = 'pcs', values_to = 'var') %>% 
        data.frame(.)
      ggplot(to.plot.pc, aes(x = to.plot.pc[,1], y = var)) +
        geom_boxplot()+
        facet_wrap(~pcs, scales = 'free') +
        ylab('PC') +
        xlab(input$pca_var2) +
        ggtitle('Raw counts') +
        theme(
          panel.background = element_blank(),
          plot.title = element_text(size = 22), 
          axis.line = element_line(colour = 'black', size = 1),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14),
          strip.text.x = element_text(size = 20),
          strip.text = element_text(size = 22))
      
      
    }
  })
  
  ## TAB 5 - NCG - COR ------
  output$corr_plot_2 <- renderPlot({
    req(input$pca_corr2)
    
    nPCs <- input$nPCs2
    pca.cancer.tcga <- values$pca.cancer.tcga2
    corr_var <- selected_corr_var2()
    # Color Pallette Set
    dataSets.colors <- wesanderson::wes_palette(n = 4,
                                                name = "GrandBudapest1")[c(1, 2, 4)]
    names(dataSets.colors) <- c('Raw counts',
                                'FPKM',
                                'FPKM.UQ')
    
    if (is.numeric(corr_var)) {
      lreg.ls.cancer.tcga <- lapply(tcga.harmonized,
                                    function(x) {
                                      pcs <- pca.cancer.tcga[[x]]$sing.val$u
                                      tcga.ls.rSquared <- sapply(1:nPCs,
                                                                 function(y) {
                                                                   lm.ls <- summary(lm(corr_var ~ pcs[, 1:y]))$r.squared
                                                                 })
                                    })
      names(lreg.ls.cancer.tcga) <- tcga.harmonized
      
      # Create Dataframe for corrplot
      ls.lnreg.normAssess <- data.frame(
        Raw.counts = lreg.ls.cancer.tcga$HTseq_counts,
        FPKM = lreg.ls.cancer.tcga$HTseq_FPKM,
        FPKM.UQ = lreg.ls.cancer.tcga$HTseq_FPKM.UQ,
        pcs = c(1:nPCs)
      )
      
      # Modify DF for corrplot
      ls.lnreg.normAssess <- ls.lnreg.normAssess %>%
        tidyr::pivot_longer(-pcs,
                            names_to = 'Datasets',
                            values_to = 'ls') %>%
        dplyr::mutate(Datasets = replace(Datasets,
                                         Datasets == 'Raw.counts', 'Raw counts')) %>%
        dplyr::mutate(Datasets = factor(
          x = Datasets,
          levels = c('Raw counts', 'FPKM', 'FPKM.UQ')
        )) %>%
        data.frame(.)
      
      # Visualise
      ggplot(ls.lnreg.normAssess, aes(x = pcs, y = ls, group = Datasets)) +
        geom_line(aes(color = Datasets), size = 1) +
        geom_point(aes(color = Datasets), size = 3) +
        xlab('PCs') + ylab (expression("R" ^ "2")) +
        scale_color_manual(
          values = c(dataSets.colors[1:3]),
          labels = c('Raw counts', 'FPKM', 'FPKM.UQ')
        ) +
        scale_x_continuous(breaks = (1:nPCs),
                           labels = c('PC1', paste0('PC1:', 2:nPCs))) +
        scale_y_continuous(breaks = scales::pretty_breaks(n = 5),
                           limits = c(0, 1)) +
        theme(
          panel.background = element_blank(),
          axis.line = element_line(colour = 'black', size = 1),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(
            size = 12,
            angle = 45,
            hjust = 1
          ),
          axis.text.y = element_text(size = 12),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 14),
          strip.text.x = element_text(size = 10)
        ) +
        ggtitle(input$pca_corr)
    }
    
    else {
      cca.time.year <-
        lapply(data.set.names,
               function(x) {
                 sapply(1:nPCs,
                        function(y) {
                          cca.pam50 <- stats::cancor(x = pca.cancer.tcga[[x]]$sing.val$u[, 1:y, drop = FALSE],
                                                     y = corr_var)
                          1 - prod(1 - cca.pam50$cor ^ 2)
                        })
                 
               })
      names(cca.time.year) <- data.set.names
      
      # Create DF for corrplot
      cca.time.year <- data.frame(
        Raw.counts = cca.time.year$HTseq_counts,
        FPKM = cca.time.year$HTseq_FPKM,
        FPKM.UQ = cca.time.year$HTseq_FPKM.UQ,
        pcs = c(1:nPCs)
      )
      
      # Modify DF for corrplot
      cca.time.year <- cca.time.year %>%
        tidyr::pivot_longer(-pcs,
                            names_to = 'Datasets',
                            values_to = 'ls') %>%
        dplyr::mutate(Datasets = replace(Datasets,
                                         Datasets == 'Raw.counts', 'Raw counts')) %>%
        dplyr::mutate(Datasets = factor(
          x = Datasets,
          levels = c('Raw counts', 'FPKM', 'FPKM.UQ')
        )) %>%
        data.frame(.)
      
      # Visualise
      ggplot(cca.time.year, aes(x = pcs, y = ls, group = Datasets)) +
        geom_line(aes(color = Datasets), size = 1) +
        geom_point(aes(color = Datasets), size = 3) +
        xlab('PCs') + ylab (expression("Vector Correlation")) +
        scale_color_manual(
          values = c(dataSets.colors),
          labels = c('Raw counts', 'FPKM', 'FPKM.UQ')
        ) +
        scale_x_continuous(breaks = (1:nPCs),
                           labels = c('PC1', paste0('PC1:', 2:nPCs))) +
        scale_y_continuous(breaks = scales::pretty_breaks(n = 5),
                           limits = c(0, 1)) +
        theme(
          panel.background = element_blank(),
          axis.line = element_line(colour = 'black', size = 1),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(
            size = 12,
            angle = 45,
            hjust = 1
          ),
          axis.text.y = element_text(size = 12),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 14),
          strip.text.x = element_text(size = 10)
        ) +
        ggtitle(input$pca_corr2)
    }
  })
  
  ## TAB 5 - PRPS------
  output$PRPS_map <- renderPlot({
    validate(need(
      input$sample_type %in% c("BRCA (Breast invasive carcinoma)",
                               "COAD (Colon adenocarcinoma)",
                               "LUAD (Lung adenocarcinoma)",
                               "READ (Rectum adenocarcinoma)"), 
      message = 'Please select either BRCA, READ, COAD or LUAD for RUV-III')
    )
    data <- selected_sampletype() 
    sample.info <- as.data.frame(SummarizedExperiment::colData(data))
    keep.samples <- sample.info$Subtypes %in% input$prps_bio
    data <- data[, keep.samples]
    sample.info <- as.data.frame(SummarizedExperiment::colData(data))
    #sample.info$biology <- sample(letters[1:4], nrow(sample.info), replace = TRUE)
    
    if ("Year" %in% input$prps_batch & "Plates" %in% input$prps_batch) {
      sample.info$new.batch <- paste0(
        sample.info$Year,
        '_',
        sample.info$Plates
      )
    }
    else{
      sample.info$new.batch <- sample.info[,input$prps_batch]
    }
    
    df_count <- sample.info %>%
      dplyr::count(new.batch, Subtypes)
    
    df_count$use <- 'Un-selected'
    df_count$use[df_count$n >= input$minSamplesPerBatchPS] <- 'Selected'
    
    ggplot(df_count, aes(x = new.batch, y = Subtypes)) +
      geom_count(aes(color = use)) +
      geom_text(aes(
        label = n,
        hjust = 0.5,
        vjust = 0.5
      )) +
      xlab('Years-plates') +
      ylab('Biological groups') +
      theme_bw()+
      theme(
        axis.line = element_line(colour = 'black', size = .85),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        plot.title = element_text(size = 15),
        axis.text.x = element_text(
          size = 12,
          angle = 45,
          hjust = 1
        ),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 10),
        legend.position = 'none'
      )
  })
  
  ## TAB 5 - GENE ---------
  output$gene_corr_plot_2 <- renderPlot({
    filter4 <- ruv_iii_gen()
    raw.count <-
      as.data.frame(SummarizedExperiment::assay(filter4, 'HTseq_counts'))
    sample.info <-
      as.data.frame(SummarizedExperiment::colData(filter4))
    
    raw.count <- as.matrix(raw.count)
    
    if (input$gene_corr2 == "Library Size") {
      corr_var <- log2(colSums(raw.count))
    }
    else if (input$gene_corr2 == "Purity") {
      corr_var <- sample.info$Purity_singscore
    }
    else if (input$gene_corr2 == "Year") {
      corr_var <- sample.info$Year
    }
    else if (input$gene_corr2 == "Plate") {
      corr_var <- sample.info$Plates
    }
    
    if (is.numeric(corr_var)) {
      genes.ls <- .corr.gene.variable(
        expr.data = raw.count,
        is.log = FALSE,
        variable = corr_var,
        method = 'spearman',
        n.cores = 5
      )
      par(mfrow = c(1, 2))
      hist(genes.ls$rho)
      plot(corr_var, log2(raw.count[input$gene_name2 ,] + 1))
    }
    else {
      time.effects <- .Ftest(
        data = raw.count,
        variable = corr_var,
        is.log = F,
        n.cores = 5
      )
      par(mfrow = c(1, 2))
      hist(log2(time.effects$FValue))
      boxplot(log2(raw.count[input$gene_name2 ,] + 1) ~ corr_var)
    }
    
  })
  
  output$prps_summary <- renderPrint({
    req(input$prps_gen)
    prps <- values$prps
    cat('After above-filterarion', "\n")
    cat('ps.batch: ', dim(prps$ps.batch), "\n")
    cat('ps.ls: ', dim(prps$ps.ls), "\n")
    cat('ps.purity: ', dim(prps$ps.purity), "\n")
  
  })
  
}

# Launch App --------------------------------------------------------------
shinyApp(ui = ui, server = server)

