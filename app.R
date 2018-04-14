library(shiny)
library(DT)
library(readxl)
library(tidyverse)

# read in sheet 1 data and make better names

associations <- read_csv("df_all_annotated.csv")
associations <- associations %>%
  rename(Func = Func.knownGene, Gene = Gene.knownGene) %>%
  rename(Marker = SNP) %>%
  rename(dbSNP = avsnp147)

# create unique chromosomes
unique_chromsomes <- associations %>%
  select(Chr) %>%
  distinct() %>%
  arrange(Chr)

# create unique genes

genes <- associations$Gene
refgenes <- associations$Gene.refGene
unique_genes <- unique(c(genes, refgenes))

# create unique illumina_ids
illumina_id <- associations$Marker
dbSNP <- associations$dbSNP
unique_markers <- unique(c(illumina_id, dbSNP))

# Define UI for application 
ui <- fluidPage(
  titlePanel("Cedars IBD Known Associations"),
  sidebarLayout(
    sidebarPanel(
      
      #p value filter
      sliderInput(inputId = "p_value",
                  label = "Associations with P Value <= 10^-",
                  min = 1,
                  max = 15,
                  value = 1,
                  step = 1),
      
      # selector for search method
      radioButtons(inputId = "searchby",
                   label = "Search By",
                   choices = c("Gene" ="gene", "Position" = "position",
                              "SNP" = "snp", "Upload Gene File" = "gene_file",
                              "Upload SNP File" = "snp_file"),
                   selected = "gene"),
      
      #conditional selector if search method is gene
      conditionalPanel(
        condition = "input.searchby == 'gene'",
        selectizeInput("genelist", "Enter Genes of Interest", choices = NULL, multiple = TRUE, 
                       options = list(placeholder = 'enter gene names',
                                      splitOn = I("(function() { return /[,;]/; })()"),
                                      create = I("function(input, callback) {
                                                              return { 
                                                              value: input,
                                                              text: input
                                                              };
                                                              }")))
        ),
      
      #conditional selector if search method is position
      conditionalPanel(
        condition = "input.searchby == 'position'",
        selectInput(inputId = "chromosome_choice",
                    label = "Chromosome",
                    choices = unique_chromsomes$Chr),
        numericInput("min_bp", "From BP", 0),
        numericInput("max_bp", "To BP", 0)
        ),
      
      #conditional selector if search method is snp
      conditionalPanel(
        condition = "input.searchby == 'snp'",
        selectizeInput("snplist", "Enter SNPs of Interest", choices = NULL, multiple = TRUE, 
                       options = list(placeholder = 'enter snp names',
                                      splitOn = I("(function() { return /[,;]/; })()"),
                                      create = I("function(input, callback) {
                                                              return { 
                                                              value: input,
                                                              text: input
                                                              };
                                                              }")))
      ),
      
      #conditional selector if search method is file
      conditionalPanel(
        condition = "input.searchby == 'gene_file' | input.searchby == 'snp_file'",
        # Input: Checkbox if file has header ----
        checkboxInput("header", "Check if your Column is Named", FALSE),
        
        # Input: Select separator ----
        radioButtons("sep", "Separator",
                     choices = c(Comma = ",",
                                 Semicolon = ";",
                                 Tab = "\t"),
                     selected = ","),
        
        # Input: Select quotes ----
        radioButtons("quote", "Are your Character Strings Quoted?",
                     choices = c(None = "",
                                 "Double Quote" = '"',
                                 "Single Quote" = "'"),
                     selected = ''),
        # Input: Select a file ----
        fileInput("file1", "Choose CSV File",
                  multiple = FALSE,
                  accept = c("text/csv",
                             "text/comma-separated-values,text/plain",
                             ".csv"))
      ),
      # snp location search
      checkboxGroupInput(inputId = "snp_location",
                         label = "Acceptable SNP Location",
                         choices = c(unique(associations$Func)),
                         selected = c(unique(associations$Func)))
    ),
    # Show a table
    mainPanel(
      tabsetPanel(
        tabPanel("By Gene Search", 
                 DT::dataTableOutput("table_gene"),
                 downloadButton("download_table_gene", "Download")),
        tabPanel("By Gene-Pheno Search", DT::dataTableOutput("table_gene_pheno"),
                 downloadButton("download_table_gene_pheno", "Download")),
        tabPanel("By Gene-Pheno-SNP Search", DT::dataTableOutput("table_gene_pheno_snp"),
                 downloadButton("download_table_gene_pheno_snp", "Download")),
        tabPanel("Information", includeMarkdown("README.md")))
    )
  )
)


  

# Define server logic 
server <- function(input, output, session) {
  
  search_choice <- reactive(input$searchby)
  
  associations_filtered <- reactive({
    if(search_choice() == "gene") {
      associations %>%
        filter(P <= 10^-input$p_value) %>%
        filter(Func %in% input$snp_location) %>%
        filter(Gene %in% input$genelist | Gene.refGene %in% input$genelist)
        }
        
    else if(search_choice() == "position") {
      associations %>%
        filter(P <= 10^-input$p_value) %>%
        filter(Func %in% input$snp_location) %>%
        filter(Chr == input$chromosome_choice) %>%
        filter(between(Start, input$min_bp, input$max_bp))
        }
    else if(search_choice() == "snp"){
      associations %>%
        filter(P <= 10^-input$p_value) %>%
        filter(Func %in% input$snp_location) %>%
        filter(Marker %in% input$snplist | dbSNP %in% input$snplist)
    }
    
    else if(search_choice() == "gene_file"){
      
      genes_from_file <- uploaded_file()
      
      associations %>%
        filter(P <= 10^-input$p_value) %>%
        filter(Func %in% input$snp_location) %>%
        filter(Gene %in% genes_from_file[,1] | Gene.refGene %in% genes_from_file[,1])
    }
    
    else if(search_choice() == "snp_file"){
      
      snps_from_file <- uploaded_file()
      
      associations %>%
        filter(P <= 10^-input$p_value) %>%
        filter(Func %in% input$snp_location) %>%
        filter(Marker %in% snps_from_file[,1] | dbSNP %in% snps_from_file[,1])
    }
  }
  )
  
  # table_genes ----------------------------------------------------------------
  output$table_gene<- DT::renderDataTable(associations_filtered() %>%
      group_by(Population, Gene, Gene.refGene) %>%
      summarise(n_phenos = length(unique(PHENOTYPE)),
                phenos = paste(unique(PHENOTYPE), collapse = ", "),
                min_p = min(P),
                max_OR_Z_B = max(OR_Z_B),
                max_OR_Z_B = min(max_OR_Z_B)) %>%
      arrange(desc(n_phenos)), filter = "bottom", options = 
        list(pageLength = 15, lengthMenu = c(15, 30, 60)))
  
  # table_gene_pheno------------------------------------------------------------
  output$table_gene_pheno<- DT::renderDataTable(associations_filtered() %>%
      group_by(Population, Gene, Gene.refGene, PHENOTYPE, Analyst, Year) %>%
      summarise(n_sig_snps = length(unique(Marker)),
                Markers = paste(unique(Marker), collapse = ", "),
                dbSNPs = paste(unique(dbSNP), collapse = ", "),
                min_p = min(P),
                max_OR_Z_B = max(OR_Z_B),
                min_OR_Z_B = min(OR_Z_B)) %>%
      arrange(desc(n_sig_snps)), filter = "bottom", options = 
        list(pageLength = 15, lengthMenu = c(15, 30, 60)))
  
  # table_gene_pheno_snp---------------------------------------------------------
  output$table_gene_pheno_snp <- DT::renderDataTable(
    associations_filtered() %>%
      group_by(Population, Gene, Gene.refGene, PHENOTYPE, Analyst, Year, Marker, dbSNP) %>%
      summarise(p_value = P,
                Chromosome = Chr,
                A1 = A1,
                n_miss = NMISS,
                OR_Z_B = OR_Z_B,
                base_pair = Start,
                snp_location = Func),
    filter = "bottom",
    options = list(pageLength = 15, lengthMenu = c(15, 30, 60)))
  
  # ui_genes-------------------------------------------------------------------
  updateSelectizeInput(session, 'genelist', choices = unique_genes, server = TRUE)
  updateSelectizeInput(session, 'snplist', choices = unique_markers, server = TRUE)
  
  # uplaod list of genes or snps
  uploaded_file <- reactive({
    inFile <- input$file1
    if (is.null(inFile))
      return(NULL)
    data <- read.csv(inFile$datapath, header = input$header, sep = input$sep, quote = input$quote)
    data
  })
    
  
  # Downloadable csv of selected dataset ---------------------------------------
  # download on tab gene
  output$download_table_gene <- downloadHandler(
    filename = function() {
      paste(Sys.Date(), "associations-table-gene", ".csv", sep = "-")
    },
    content = function(file) {
      
      write.csv(associations_filtered() %>%
                  group_by(Population, Gene, Gene.refGene) %>%
                  summarise(n_phenos = length(unique(PHENOTYPE)),
                            phenos = paste(unique(PHENOTYPE), collapse = ", "),
                            min_p = min(P),
                            max_OR_Z_B = max(OR_Z_B),
                            max_OR_Z_B = min(max_OR_Z_B)) %>%
                  arrange(desc(n_phenos)), file, row.names = FALSE)
    }
  )
  
  # download on tab gene-pheno
  output$download_table_gene_pheno <- downloadHandler(
    filename = function() {
      paste(Sys.Date(), "associations-table-gene-pheno", ".csv", sep = "-")
    },
    content = function(file) {
      
      write.csv(associations_filtered() %>%
                  group_by(Population, Gene, Gene.refGene) %>%
                  summarise(n_phenos = length(unique(PHENOTYPE)),
                            phenos = paste(unique(PHENOTYPE), collapse = ", "),
                            min_p = min(P),
                            max_OR_Z_B = max(OR_Z_B),
                            max_OR_Z_B = min(max_OR_Z_B)) %>%
                  arrange(desc(n_phenos)), file, row.names = FALSE)
    }
  )
  
  # download on tab gene-pheno-snp
  output$download_table_gene_pheno_snp <- downloadHandler(
    filename = function() {
      paste(Sys.Date(), "associations-table-gene-pheno-snp", ".csv", sep = "-")
    },
    content = function(file) {
      
      write.csv(associations_filtered() %>%
                  group_by(Population, Gene, Gene.refGene, PHENOTYPE, Analyst, Year, Marker, dbSNP) %>%
                  summarise(p_value = P,
                            Chromosome = Chr,
                            A1 = A1,
                            n_miss = NMISS,
                            OR_Z_B = OR_Z_B,
                            base_pair = Start,
                            snp_location = Func), file, row.names = FALSE)
    }
  )
  
}

# Run the application 
shinyApp(ui = ui, server = server)

