library(shiny)
library(DT)
library(readxl)
library(tidyverse)

# read in sheet 1 data and make better names

associations <- read_excel("data/IBD1_IBD2_IBD9_IBD10_top ten each category.xlsx", 
                           sheet = 1)
colnames(associations) <- make.names(colnames(associations))

# split GENE.LOCATION.CodingStatus
associations <- associations %>%
  separate(GENE.LOCATION.CodingStatus, into = c("gene", "snp_location", "coding_status"),
           sep = "=") %>%
  rename(population = CD_pheno) %>%
  rename()

# create unique chromosomes
unique_chromsomes <- associations %>%
  select(CHR) %>%
  distinct() %>%
  arrange(CHR)

# create unique genes
unique_genes <- associations %>%
  select(gene) %>%
  distinct(gene) 

# create unique rsids
unique_rsids <- associations %>%
  select(RSID) %>%
  distinct()


# Define UI for application 
ui <- fluidPage(
  titlePanel("Cedars IBD Known Associations"),
  sidebarLayout(
    sidebarPanel(
      sliderInput(inputId = "maf",
                  label = "SNP with Minor Allele Frequency >=",
                  min = 0,
                  max = 0.2,
                  value = 0,
                  step = 0.01),
      checkboxGroupInput(inputId = "snp_location",
                         label = "Acceptable SNP Location",
                         choices = c(unique(associations$snp_location)),
                         selected = c(unique(associations$snp_location))),
      sliderInput(inputId = "p_value",
                  label = "Associations with P Value <= 10^-",
                  min = 1,
                  max = 15,
                  value = 1,
                  step = 1),
      radioButtons(inputId = "searchby",
                   label = "Search By",
                   choices = c("Gene" ="gene", "Position" = "position",
                                                         "SNP" = "snp"),
                   selected = "gene"),
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
      conditionalPanel(
        condition = "input.searchby == 'position'",
        selectInput(inputId = "chromosome_choice",
                    label = "Chromosome",
                    choices = unique_chromsomes$CHR),
        numericInput("min_bp", "From BP", 0),
        numericInput("max_bp", "To BP", 0)
        ),
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
      )
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
        filter(snp_location %in% input$snp_location) %>%
        filter(MAF >= input$maf) %>%
        filter(gene %in% input$genelist)
        }
        
    else if(search_choice() == "position") {
      associations %>%
        filter(P <= 10^-input$p_value) %>%
        filter(snp_location %in% input$snp_location) %>%
        filter(MAF >= input$maf) %>%
        filter(CHR == input$chromosome_choice) %>%
        filter(between(BP, input$min_bp, input$max_bp))
        }
    else if(search_choice() == "snp"){
      associations %>%
        filter(P <= 10^-input$p_value) %>%
        filter(snp_location %in% input$snp_location) %>%
        filter(MAF >= input$maf) %>%
        filter(RSID %in% input$snplist)
    }
  }
  )
    
  
  # table_genes ----------------------------------------------------------------
  output$table_gene<- DT::renderDataTable(associations_filtered() %>%
      group_by(population, gene) %>%
      summarise(n_phenos = length(unique(PHENOTYPE)),
                phenos = paste(unique(PHENOTYPE), collapse = ", "),
                min_p = min(P),
                max_OR = max(OR),
                min_OR = min(OR)) %>%
      arrange(desc(n_phenos)), filter = "bottom", options = 
        list(pageLength = 25, lengthMenu = c(25, 50, 100)))
  
  # table_gene_pheno------------------------------------------------------------
  output$table_gene_pheno<- DT::renderDataTable(associations_filtered() %>%
      group_by(population, gene, PHENOTYPE) %>%
      summarise(n_sig_snps = length(unique(RSID)),
                SNPs = paste(unique(RSID), collapse = ", "),
                min_p = min(P),
                max_OR = max(OR),
                min_OR = min(OR)) %>%
      arrange(desc(n_sig_snps)), filter = "bottom", options = 
        list(pageLength = 25, lengthMenu = c(25, 50, 100)))
  
  # table_gene_pheno_snp---------------------------------------------------------
  # talk with talin about rs2066844, non unique SNP
  output$table_gene_pheno_snp <- DT::renderDataTable(
    associations_filtered() %>%
      group_by(population, gene, PHENOTYPE, RSID, SNP) %>%
      summarise(p_value = P,
                chromosome = CHR,
                a1 = A1,
                n_miss = NMISS,
                odds_ratio = OR,
                base_pair = BP,
                snp_location = snp_location),
    filter = "bottom",
    options = list(pageLength = 25, lengthMenu = c(25, 50, 100)))
  
  # ui_genes-------------------------------------------------------------------
  updateSelectizeInput(session, 'genelist', choices = unique_genes$gene, server = TRUE)
  updateSelectizeInput(session, 'snplist', choices = unique_rsids$RSID, server = TRUE)
  
  # Downloadable csv of selected dataset ---------------------------------------
  output$download_table_gene <- downloadHandler(
    filename = function() {
      paste(Sys.Date(), "associations-table-gene", ".csv", sep = "-")
    },
    content = function(file) {
      
      write.csv(associations_filtered() %>%
                  group_by(population, gene) %>%
                  summarise(n_phenos = length(unique(PHENOTYPE)),
                            phenos = paste(unique(PHENOTYPE), collapse = ", "),
                            min_p = min(P),
                            max_OR = max(OR),
                            min_OR = min(OR)) %>%
                  arrange(desc(n_phenos)), file, row.names = FALSE)
    }
  )
  
  output$download_table_gene_pheno <- downloadHandler(
    filename = function() {
      paste(Sys.Date(), "associations-table-gene-pheno", ".csv", sep = "-")
    },
    content = function(file) {
      
      write.csv(associations_filtered() %>%
                  group_by(population, gene, PHENOTYPE) %>%
                  summarise(n_sig_snps = length(unique(RSID)),
                            SNPs = paste(unique(RSID), collapse = ", "),
                            min_p = min(P),
                            max_OR = max(OR),
                            min_OR = min(OR)) %>%
                  arrange(desc(n_sig_snps)), file, row.names = FALSE)
    }
  )
  
  output$download_table_gene_pheno_snp <- downloadHandler(
    filename = function() {
      paste(Sys.Date(), "associations-table-gene-pheno-snp", ".csv", sep = "-")
    },
    content = function(file) {
      
      write.csv(associations_filtered() %>%
                  group_by(population, gene, PHENOTYPE, RSID, SNP) %>%
                  summarise(p_value = P,
                            chromosome = CHR,
                            a1 = A1,
                            n_miss = NMISS, 
                            odds_ratio = OR,
                            basse_pair = BP,
                            snp_location = snp_location), file, row.names = FALSE)
    }
  )
}

# Run the application 
shinyApp(ui = ui, server = server)

