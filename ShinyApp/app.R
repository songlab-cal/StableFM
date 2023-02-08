#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
# This version dated: 1/23/23 
# 
# To (re)deploy:
# Remove rsconnect directory
# Go to https://www.shinyapps.io/admin/#/dashboard and archive + delete
# Follow instructions on launching new app. Last step is to run:
# rsconnect::deployApp('/clusterfs/nilah/alan/Shiny/stability_v0')

library(shiny)
library(shinythemes)
library(shinyjs) # [TESTING]
library(dplyr)
library(ggplot2)
#source('/clusterfs/nilah/alan/Shiny/test/auxiliary.R')
source('auxiliary.R')

chromosomes <- 1:22
gene_names <- vector('list', length = 22)
gene_names <- lapply(1:22, readFile <- function(x) {
  #read.table(paste0('/clusterfs/nilah/alan/peer_phenotypes/chr', x, '_gene_names.txt'))$V1
  #read.table(paste0('peer_phenotypes/chr', x, '_gene_names.txt'))$V1
  read.csv(paste0('peer_phenotypes/chr', x, '_gene_names.csv'))$INPUT
})
names(gene_names) <- 1:22
# List of all functional annotations, group by functional annotation class
func_annots <- c('Ensembl', 'FAVOR', 'FIRE', 'Enformer')
func_annot_names <- lapply(func_annots, readFile <- function(x) {
  #readLines(paste0('/clusterfs/nilah/alan/all-func-annots/', x, '_annot_names.txt'))
  readLines(paste0('all-func-annots/', x, '_annot_names.txt'))
})
names(func_annot_names) <- func_annots
# Links to short descriptions for 201 functional annotations
func_annot_descriptions <- lapply(func_annots, readFile <- function(x) {
  #readLines(paste0('/clusterfs/nilah/alan/all-func-annots/', x, '_annot_descriptions.txt'))
  readLines(paste0('all-func-annots/', x, '_annot_descriptions.txt'))
})
names(func_annot_descriptions) <- func_annots
for (func in func_annots) {
  names(func_annot_descriptions[[func]]) <- func_annot_names[[func]]
}
# Links to references for 201 functional annotations
func_annot_links <- lapply(func_annots, readFile <- function(x) {
  #readLines(paste0('/clusterfs/nilah/alan/all-func-annots/', x, '_annot_links.txt'))
  readLines(paste0('all-func-annots/', x, '_annot_links.txt'))
})
names(func_annot_links) <- func_annots
for (func in func_annots) {
  names(func_annot_links[[func]]) <- func_annot_names[[func]]
}

# Define UI for application that draws a histogram
ui <- fluidPage(
  useShinyjs(), # [TESTING]
  theme = shinytheme("cosmo"),
  # Application title
  titlePanel("Visualizing and Interpreting Statistical Fine Mapping"),
  p("Created by: Alan Aw, Lionel Jin"),
  p("To learn more about our work, check out our preprint: A.J. Aw, L.C. Jin, N.M. Ioannidis and Y.S. Song (2023+) 'The impact of stability considerations on statistical fine-mapping' (in preparation). We recommend browsing this app on a Desktop for optimal user experience."),
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      h3("Pick a Chromosome and Gene"),
      h4("If looking up by gene name, select chromosome in which the gene is transcribed."),
      fluidRow(
        #h3("Pick a Chromosome and Gene"),
        column(3,selectInput(inputId = "table_selector", 
                             label = "Chromosome", 
                             choices = chromosomes, 
                             selected = "1")),
        column(8,selectInput(inputId = "gene_selector", 
                             label = "Gene", 
                             choices = ""))
      ),
      h3("Summary of Gene Metadata"),
      h4("Gene Information"),
      uiOutput("GeneInfoSummary"),
      fluidRow(column(2,uiOutput(outputId = "GTEx")),
               column(2,uiOutput(outputId = "Ensembl"))),
      h4("Expression Distribution"),
      plotOutput(outputId = "geneExpressionPlot"),
      actionButton("geneExpressionShowBtn", "Show Summary"), # [TESTING]
      actionButton("geneExpressionHideBtn", "Hide Summary"), # [TESTING]
      hidden(tableOutput("geneExpressionSumTab"))
    ),
    
    # MAIN PANEL
    mainPanel(
      tabsetPanel(
        # Panel: Single-gene fine-mapping
        tabPanel("Single Gene Comparison",
                 # Let user choose Potential Set
                 fluidRow(
                   h3("Pick a Potential Set"),
                   h4("If unsure, use 1. See Methodology tab for the definition of a Potential Set."),
                   # Input: Select dataset (in the future, build upload functionality)
                   selectInput(inputId = "potential_set",
                               label = "",
                               choices = 1:3)
                 ),
                 # Show table of PICS2 statistics for top vs stable variant
                 fluidRow(
                   tabPanel(title = "PICS2 Results",
                            p("Below summarizes the gene-specific fine-mapping results for the gene selected in left panel and the Potential Set selected above."),
                            fluidRow(column(6,h3("1. Summarize Fine-Mapping Output"),
                                            uiOutput("TopStableMatch"),
                                            tableOutput("PICS2TopVSStable"),
                                            h4("SNP Information (dbSNP)"),
                                            fluidRow(column(3,uiOutput(outputId = "PICS2TopDBSNP")),
                                                     column(3,uiOutput(outputId = "PICS2StableDBSNP")))),
                                     column(6,h3("2. Stable Variant Diversity-Related Statistics"),
                                            h4("Statistics capturing the", 
                                              shiny::tags$i("degree of stability"),
                                              "of the stable variant identified"), 
                                            tableOutput("PICS2StableStats"),
                                            actionButton("StabStatsShowBtn", "Show Details"), # [TESTING]
                                            actionButton("StabStatsHideBtn", "Hide Details"), # [TESTING]
                                            hidden(shiny::tags$ul(id = "stab_stats",
                                              shiny::tags$li(p("The", shiny::tags$b("No. Slices Reporting Variant with Positive Probability"),
                                                               "measures how often the stable variant has a posterior probability larger than zero, when performing fine-mapping across all user-defined slices (five in total for this study).")), 
                                              shiny::tags$li(p("Based on the quantity measured above,", 
                                                               shiny::tags$b("Max Slice-Slice F_ST"), 
                                                               "computes the maximum",
                                                               HTML(paste0(shiny::tags$i('F'),shiny::tags$sub('ST'))),
                                                               "(a measure of genetic differentiaton; see",
                                                               shiny::tags$a("Holsinger and Weir, 2009", href="https://www.nature.com/articles/nrg2611"),
                                                               ") between any pair of slices for which the variant has a positive posterior probability. Larger values imply greater genetic differentiation.")), 
                                              shiny::tags$li(p("Similarly,",
                                                               shiny::tags$b("Max Slice-Slice AF Difference"),
                                                               "computes the maximum",
                                                               shiny::tags$i('mean allele frequency difference'),
                                                               "(a measure of genetic heterogeneity) between any pair of slices for which the variant has a positive posterior probability. Larger values imply greater genetic heterogeneity."))
                                            )))),
                            fluidRow(column(6, h3("3. Compare Functional Annotations"),
                                            fluidRow(column(6, selectInput(inputId = "func_annot_class",
                                                                           label = "Pick an Annotation Class",
                                                                           choices = func_annots)),
                                                     column(6,selectInput(inputId = "func_annot_selector", 
                                                                          label = "Choose Annotation", 
                                                                          choices = ""))),
                                            fluidRow(column(12,
                                                            h4("Annotation Description"),
                                                            textOutput(outputId = "FuncAnnotDescription"),
                                                            uiOutput("EnsemblNote"),
                                                            uiOutput(outputId = "FuncAnnotReference"),
                                                            tableOutput("FuncAnnotTable")))))
                            
                   ))
        ), 
        tabPanel("Genome-Wide Comparison",
                 h3("Make selections below"),
                 fluidRow(column(4, selectInput(inputId = "gw_func_annot_class",
                                                label = "Choose Annotation Class",
                                                choices = func_annots)),
                          column(4,selectInput(inputId = "gw_func_annot_selector", 
                                               label = "Choose Annotation", 
                                               choices = "")),
                          column(4,selectInput(inputId = "gw_ps",
                                               label = "Choose Potential Set",
                                               choice = 1:3))),
                 uiOutput("EnformerFormula"),
                 h4("Annotation Description"),
                 textOutput(outputId = "funcAnnotDescription"),
                 uiOutput("GWEnsemblNote"),
                 uiOutput(outputId = "funcAnnotReference"),
                 h3("Compare Matching and Non-matching Variants"),
                 p("The density plots show the distributions of functional annotation scores of matching variants (i.e., genes where residualization and stability-guided approaches agree) and non-matching variants (either the top or the stable) across all GEUVADIS genes."),
                 actionButton("funcDensityPlotShowBtn", "Show Plot"), # [TESTING]
                 actionButton("funcDensityPlotHideBtn", "Hide Plot"), # [TESTING]
                 hidden(plotOutput(outputId = "funcDensityPlot", # [TESTING]
                            height = "400px")),
                 h3("Compare Top and Stable Variants"),
                 p("The pairwise plot compares functional annotation scores of the top variant and the stable variant across all GEUVADIS genes, whenever the two variants disagree.  (Genes with matching top and stable variants not shown.)"),
                 actionButton("funcAnnotPlotShowBtn", "Show Plot"), # [TESTING]
                 actionButton("funcAnnotPlotHideBtn", "Hide Plot"), # [TESTING]
                 hidden(plotOutput(outputId = "funcAnnotPlot", # [TESTING]
                            height = "400px"))
        ),
        tabPanel("Methodology",
                 h3("Fine-Mapping Approach"),
                 p("We implement PICS2, a non-parametric statistical fine-mapping approach introduced in", 
                   a("Taylor et al. (2021)", href = "https://doi.org/10.1093/bioinformatics/btab122"),
                   " and ",
                   a("Farh et al. (2015)", href = "https://www.nature.com/articles/nature13835"),
                   ". PICS2 performs fine mapping by reporting the posterior probability of each input variant, based on a permutation procedure that preserves both LD relationships between the input SNPs and the marginal association signal of that input variant. The figure below provides a visualization."),
                 br(),
                 HTML('<center><img src="pics2_tikz.jpg" width="771.25"></center>'),
                 #imageOutput("PICS2Graphic"),
                 p("Our results are generated by running PICS2 on GEUVADIS PEER-factor-normalized gene expression data, and performing fine mapping on variants within 1Mb upstream and downstream of the gene's canonical TSS ('cis fine mapping')."),
                 h3("Top versus Stable"),
                 p("We investigate two approaches to running PICS2 (see accompanying paper for details). One accounts for population structure by residualizing the phenotype against principal components (",
                   shiny::tags$b("Top"),
                   "); the other uses available external annotations and runs PICS2 on each slice, while prioritizing variants that appear in as many slices and with the largest posterior probability (",
                   shiny::tags$b("Stable"),
                   ")."),
                 h3("Reporting Posterior Probabilities"),
                 p("Our tables report two posterior probabilities corresponding to the top and the stable variant. For the top variant, this quantity is the output of PICS2 on the pooled sample", 
                   shiny::tags$i("after residualization"), 
                   "(using the top PCs) is performed on the phenotype. For the stable variant, this quantity is the output of PICS2 on the pooled sample, with the", 
                   shiny::tags$i("choice"), 
                   "of causal variant depending on PICS2 outputs for the slices."),
                 h3("Potential Set"),
                 p("In our paper, we define a potential set as a set of variants with positive posterior probability after a number of iterations of the permutation algorithm. In our implementation, we iteratively run a permutation algorithm that generates increasingly",
                   shiny::tags$i('less reliable'), "potential sets of variants. Thus, Potential Set 1 contains variants that have LD at least 0.5 with the lead variant, and can be treated as the set containing the putatively causal variant. Subsequent potential sets are obtained by removing variants from an earlier potential set and running the same permutation algorithm. Notably, subsequent potential sets do not overlap with previous ones, and are less marginally correlated with the gene expression phenotype.")
        )
      )
    )
  )
)

# Define server logic 
server <- function(input, output, session) {
  
  # Enable selection of gene matched to a selected chromosome
  observeEvent(input$table_selector, {
    message("Table event observed -- chromosome selection")
    choices <- gene_names[[input$table_selector]]
    updateSelectInput(session, "gene_selector", choices = choices)
  })
  
  # Enable selection of annotation matched to an annotation class
  observeEvent(input$func_annot_class, {
    message("Table event observed -- functional annotation selection")
    choices <- func_annot_names[[input$func_annot_class]]
    updateSelectInput(session, "func_annot_selector", choices = choices)
  })
  
  # Enable selection of annotation matched to an annotation class
  # for genome-wide visualization
  observeEvent(input$gw_func_annot_class, {
    message("Table event observed -- genome-wide functional annotation selection")
    choices <- func_annot_names[[input$gw_func_annot_class]]
    updateSelectInput(session, "gw_func_annot_selector", choices = choices)
  })
  
  # Enable hiding and revealing of functional plots (Genome-wide tab) # [TESTING]
  observeEvent(input$funcAnnotPlotShowBtn,
               {show("funcAnnotPlot")})
  observeEvent(input$funcAnnotPlotHideBtn,
               {hide("funcAnnotPlot")})
  
  observeEvent(input$funcDensityPlotShowBtn,
               {show("funcDensityPlot")})
  observeEvent(input$funcDensityPlotHideBtn,
               {hide("funcDensityPlot")})
  
  # Enable hiding and showing of stability statistics details (Single gene tab) # [TESTING]
  observeEvent(input$StabStatsShowBtn,
               {show("stab_stats")})
  observeEvent(input$StabStatsHideBtn,
               {hide("stab_stats")})
  
  # Enable hiding and showing of gene expression table (Side panel) # [TESTING]
  observeEvent(input$geneExpressionShowBtn,
               {show("geneExpressionSumTab")})
  observeEvent(input$geneExpressionHideBtn,
               {hide("geneExpressionSumTab")})
  
  # Gene Expression metadata
  output$geneExpressionPlot <- renderPlot(req(try(summarizeGeneExp(chr=input$table_selector, 
                                                                   input.gene.name=input$gene_selector)[['PLOT']])))
  output$geneExpressionSumTab <- renderTable(req(try(summarizeGeneExp(chr=input$table_selector, 
                                                               input.gene.name=input$gene_selector)[['SUM_TABLE']])), 
                                             digits = 4, rownames = TRUE)
  output$GTEx <-renderUI(a(href=paste0('https://www.gtexportal.org/home/gene/', 
                                       stringr::str_split(input$gene_selector,"[.]")[[1]][1]),
                           "GTEx",target="_blank"))
  output$Ensembl <- renderUI(a(href=paste0('http://grch37.ensembl.org/Human/Search/Results?q=', 
                                           stringr::str_split(input$gene_selector,"[.]")[[1]][1],
                                           ';site=ensembl;facet_species=Human'), 
                               "Ensembl", target="_blank"))
  
  
  # PICS comparison tables and dbSNP URLs
  PICSResult <- reactive({
    summarizePICSResult(chr=input$table_selector,
                        input.gene.name=input$gene_selector,
                        ps=input$potential_set)})
  
  output$PICS2TopVSStable <- renderTable(req(try(PICSResult()$TOP_VS_STABLE)), rownames = TRUE)
  
  output$TopStableMatch <- renderUI({
    req(try(top_variant <- PICSResult()$TOP_VS_STABLE["rsID", "Top"]))
    req(try(stable_variant <- PICSResult()$TOP_VS_STABLE["rsID", "Stable"]))
    if (top_variant == stable_variant){
      h4("Top and stable variants", strong("match"), paste0("(", top_variant, ")"))
    } else {
      h4("Top and stable variants are", strong("different"), paste0("(", top_variant, " vs ", stable_variant, ")"))
    }
  })
  
  # Create a paragraph-sentence summarizing user-selected gene
  summarizedGeneExp <- reactive({
    req(input$table_selector, input$gene_selector)
    summarizeGeneExp(chr=input$table_selector, input.gene.name=input$gene_selector)
  })
  
  output$GeneInfoSummary <- renderUI({
    gene_symbol <- summarizedGeneExp()[['GENE_SYMBOL']]
    req(try(ensembl_gene <- summarizedGeneExp()[['ENSEMBL_GENE']]))
    biotype <- summarizedGeneExp()[['BIOTYPE']]
    if (is.na(gene_symbol)) {
      p("The selected gene is ",
        strong(ensembl_gene), 
        ", with no gene name available on Biomart. Click on links below for more information.")
    } else {
      p("The selected gene is ",
        strong(ensembl_gene),
        " (Gene Name: ",
        strong(gene_symbol),
        "). Its biotype classification is: ",
        strong(biotype),
        ". Click on links below for more information.")
    }
  })
  
  output$PICS2StableStats <- renderTable(req(try(PICSResult()$STABLE_STATS)), 
                                         rownames = TRUE,
                                         digits = 5)
  output$PICS2TopDBSNP <- renderUI(a(href=req(try(PICSResult()$TOP_DBSNP)), 
                                     "Top Variant", 
                                     target="_blank"))
  output$PICS2StableDBSNP <- renderUI(a(href=req(try(PICSResult()$STAB_DBSNP)), 
                                        "Stable Variant", 
                                        target="_blank"))
  
  # PICS functional annotation 
  output$FuncAnnotTable <- renderTable(req(try(getFuncAnnotTable(chr=input$table_selector,
                                                                 input.gene.name=input$gene_selector, 
                                                                 ps=input$potential_set, 
                                                                 func.class=input$func_annot_class,
                                                                 func.annot=input$func_annot_selector))), 
                                       rownames = TRUE,
                                       digits = 5)
  output$FuncAnnotDescription <- renderText(func_annot_descriptions[[input$func_annot_class]][input$func_annot_selector])
  output$FuncAnnotReference <- renderUI(a(href=func_annot_links[[input$func_annot_class]][input$func_annot_selector], 
                                          "Reference", 
                                          target="_blank"))
  
  # Genome-wide pairwise plots of functional annotation
  output$funcDensityPlot <- renderPlot({
    # If functional annotation input has not refreshed, pick the first annotation in the correct list
    func.class=input$gw_func_annot_class
    func.annot=input$gw_func_annot_selector
    ps=input$gw_ps
    func_annot_list <- readLines(paste0('all-func-annots/', func.class, '_annot_names.txt'))
    req(func.annot %in% func_annot_list)
    # Draw density plot
    getDensityPlot(func.class, func.annot, ps)})
  
  output$funcAnnotPlot <- renderPlot({
    # If functional annotation input has not refreshed, pick the first annotation in the correct list
    func.class=input$gw_func_annot_class
    func.annot=input$gw_func_annot_selector
    ps=input$gw_ps
    func_annot_list <- readLines(paste0('all-func-annots/', func.class, '_annot_names.txt'))
    req(func.annot %in% func_annot_list)
    # Draw pairwise plot
    getPairwisePlot(func.class, func.annot, ps)})
  
  output$funcAnnotDescription <- renderText(func_annot_descriptions[[input$gw_func_annot_class]][input$gw_func_annot_selector])
  output$funcAnnotReference <- renderUI(a(href=func_annot_links[[input$gw_func_annot_class]][input$gw_func_annot_selector], 
                                          "Reference", 
                                          target="_blank"))
  
  # Explain Enformer log-root scale
  output$EnformerFormula <- renderUI({
    withMathJax(
      sprintf("Note: For Enformer annotations, original perturbations are transformed by the function \\(f(x) = \\log(x^\\frac{1}{4} + 1)\\) (log-root scale) to improve visualization.")
    )
  })
  
  # Explain enrichment scores for Ensembl (Single gene tab)
  output$EnsemblNote <- renderUI({
    req(try(annot <- input$func_annot_selector))
    if (annot %in% c('CTCF Binding Site', 
                     'Enhancer', 'Open chromatin', 
                     'Promoter Flanking Region', 
                     'Promoter', 'TF binding site')) {
      p("Note: Enrichments are calculated by allocating a score of 1 or 0 for each input variant based on its inclusion in a given annotated region, and then taking the ratio of the fine-mapped variant (top or stable) and the mean score.")
    }
  })
  
  # Explain enrichment scores for Ensembl (Genome-wide tab)
  output$GWEnsemblNote <- renderUI({
    req(try(annot <- input$gw_func_annot_selector))
    if (annot %in% c('CTCF Binding Site', 
                     'Enhancer', 'Open chromatin', 
                     'Promoter Flanking Region', 
                     'Promoter', 'TF binding site')) {
      p("Note: Enrichments are calculated by allocating a score of 1 or 0 for each input variant based on its inclusion in a given annotated region, and then taking the ratio of the fine-mapped variant (top or stable) and the mean score.")
    }
  })
  
  # Graphic summarizing PICS2
  output$PICS2Graphic <- renderImage({
    list(
      src = file.path("www/pics2_tikz.jpg"),
      #contentType = "image/jpeg",
      width = 3085/4,
      height = 1568/4
    )
  }, deleteFile = FALSE)
}

# Run the application 
shinyApp(ui = ui, 
         server = server)