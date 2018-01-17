library(shiny)
library(readxl)
library(tidyverse)

# read in sheet 1 data and make better names

associations <- read_excel("data/IBD1_IBD2_IBD9_IBD10_top ten each category.xlsx", 
                           sheet = 1)
colnames(associations) <- make.names(colnames(associations))

# split GENE.LOCATION.CodingStatus
associations <- associations %>%
  separate(GENE.LOCATION.CodingStatus, into = c("gene", "location", "coding_status"),
           sep = "=") %>%
  rename(sub_type = CD_pheno)

# Define UI for application 
ui <- fluidPage(
  titlePanel("Cedars IBD Known Associations"),
  sidebarLayout(
    sidebarPanel(
      radioButtons(inputId = "disease_class",
                   label = "Select Disease Sub-Type",
                   choices = list("Crohn's Disease" = "CD_pheno",
                                  "Ulcerative_colitits" = "UC_pheno",
                                  "CD and UC" = "both"),
                   selected = "CD_pheno"),
      sliderInput(inputId = "maf",
                  label = "SNP with Minor Allele Frequency >=",
                  min = 0,
                  max = 0.2,
                  value = 0,
                  step = 0.01),
      checkboxGroupInput(inputId = "snp_location",
                         label = "Acceptable SNP Location",
                         choices = c(unique(associations$location)),
                         selected = c(unique(associations$location))),
      sliderInput(inputId = "p_value",
                  label = "Associations with P Value <= 10^-",
                  min = 1,
                  max = 15,
                  value = 1,
                  step = 1),
      uiOutput("genes")
    ),
    # Show a table
    mainPanel(
      tabsetPanel(
        tabPanel("Gene Based", dataTableOutput("table_genes")),
        tabPanel("SNP-Gene Based", dataTableOutput("table_gene_snps")),
        tabPanel("SNP Based", dataTableOutput("table_snps"))
      )
    )
  )
)

  

# Define server logic 
server <- function(input, output, session) {
  associations_filtered <- reactive({
    
    if (input$disease_class != "both") {
      associations %>%
        filter(sub_type == input$disease_class) %>%
        filter(P <= 10^-input$p_value) %>%
        filter(gene %in% input$gene_list) %>%
        filter(location %in% input$snp_location) %>%
        filter(MAF >= input$maf)
  }
      associations %>%
        filter(P <= 10^-input$p_value) %>%
        filter(gene %in% input$gene_list) %>%
        filter(location %in% input$snp_location) %>%
        filter(MAF >= input$maf)
  
  })
  
  gene_options <- reactive({
    
    if (input$disease_class != "both") {
      associations %>%
        filter(sub_type == input$disease_class) %>%
        filter(P <= 10^-input$p_value) %>%
        filter(location %in% input$snp_location) %>%
        filter(MAF >= input$maf) %>%
        select(gene) %>%
        distinct()
    }
      associations %>%
        filter(P <= 10^-input$p_value) %>%
        filter(location %in% input$snp_location) %>%
        filter(MAF >= input$maf) %>%
        select(gene) %>%
        distinct()
    
  })
  
  # table_genes ----------------------------------------------------------------
  output$table_genes<- renderDataTable({
    
    associations_filtered() %>%
      group_by(sub_type, gene) %>%
      summarise(n_phenos = length(unique(PHENOTYPE)),
                phenos = paste(unique(PHENOTYPE), collapse = ", "),
                min_p = min(P),
                max_OR = max(OR),
                min_OR = min(OR)) %>%
      arrange(desc(n_phenos))
    
  })
  
  # table_gene_snps ------------------------------------------------------------
  output$table_gene_snps<- renderDataTable({

    associations_filtered() %>%
      group_by(sub_type, gene, PHENOTYPE) %>%
      summarise(n_sig_snps = length(unique(RSID)),
                SNPs = paste(unique(RSID), collapse = ", "),
                min_p = min(P),
                max_OR = max(OR),
                min_OR = min(OR)) %>%
      arrange(desc(n_sig_snps))
    
  })
  
  # table_snps ----------------------------------------------------------------
  # talk with talin about rs2066844, non unique SNP
  output$table_snps <- renderDataTable({
    associations_filtered() %>%
      group_by(sub_type, gene, PHENOTYPE, RSID) %>%
      summarise(p_value = P,
                chromosome = CHR,
                a1 = A1,
                n_miss = NMISS, 
                odds_ration = OR,
                bass_pair = BP,
                location = location)
  })
  
  # ui_genes-------------------------------------------------------------------
  # test gene list IL1,IL2,IL3,IL4,IL5,IL6,IL7,IL8,IL9,IL10
  
  output$genes <- renderUI({
    
    selectizeInput("gene_list", "Select Genes of Interest", choices = gene_options(), multiple = TRUE,
                   selected = NULL, options = list(placeholder = 'select gene names',
                                                   splitOn = I("(function() { return /[,;]/; })()"),
                                                   create = I("function(input, callback) {
                                                              return { 
                                                              value: input,
                                                              text: input
                                                              };
                                                              }")
                   )
    )
  }
  )
  
  updateSelectizeInput(session, 'gene_list', choices = NULL, server = TRUE)
  
}

# Run the application 
shinyApp(ui = ui, server = server)

