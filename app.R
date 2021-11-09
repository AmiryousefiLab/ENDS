#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
# 

# Load functions that work data and make plots
# library(rstudioapi)
# currentpath = rstudioapi::getActiveDocumentContext()$path
# setwd(dirname(currentpath))

library(markdown)
library(shinyWidgets)
library(periscope)
library(gridExtra)
library(shiny)

source('PreliminaryFunctions.R')
source('NPDS.R')
source('IsotonicRegressionFit.R')
source('SigmoidFit.R') 
source('MH_AMspline.R')

df = openxlsx::read.xlsx('data/Drug_response_S8.xlsx', sheet = 1)
df_list  = read_excel_allsheets('data/Drug_response_S8.xlsx')
# block = extract_dose_block(df_list, drug = '5FU', patient = 'P1', treatment = 'T1', sample = 1)
# block2 = preprocess_data(block)
df_example = read.csv('data/Example1.csv')



# Define UI for application that draws a histogram
ui <- fluidPage(
  
    withMathJax(),

    shinyjs::useShinyjs(),
    
    
    tags$header(
        HTML('<style type="text/css">
      p {color:#551414;}
      #title{

        right:16px;
        left:16px;
        height:178px;
        background-color:#FEF2DC;
        margin-top:50px;
        margin-bottom:50px;
        padding:20px;
        border-radius:300px;
        border:20px solid;
        border-color:#BD3C3C;
      }
      #tool_name{
        font-family:Copperplate Gothic Bold;
        letter-spacing:8px;
        text-align:center;
      }
      #tool_exp{
        font-family:Charlesworth;
        text-align:center;
        color:#117704;
      }
      #tool{
        margin-left: 100px;
      }
      #gbfile{
        margin-bottom:-15px;
      }
      #remove_file{
        width:100%; 
        margin-top: 25px; 
        margin-left: 0px;
      }
      #remove_layer1, #remove_layer2, #remove_layer3, #layer1_example, #layer2_example, #layer3_example{
        margin-top: 0px;
        margin-bottom: 35px;
      }
      label[for="layer1"]{
        font-weight: normal;
      }
      label[for="layer2"] {
        font-weight: normal;
      }
      label[for="layer3"] {
        font-weight: normal;
      }

      label[for="layer1_type"]{
        font-weight: normal;
      }
      label[for="layer2_type"]{
        font-weight: normal;
      }
      </style>')
    ),
    
    tags$div(id = "title",
             tags$h1(id = "tool_name", "The ENDS", style = "color:#BD3C3C;"),
             tags$h4(id = "tool_exp", tags$em(" An online suite for the Epistemic Nonparametric Drug response Scoring"),
                     style = "color:#551414;")
    ),
    
    # Color of letters in tags 
    tags$head(tags$style(type = "text/css", "a{color: #551414;}")),
    tags$div(id = "tool",
             tabsetPanel(
               tabPanel(title = strong("Home"), 
                        br(),
                        tags$div(includeMarkdown("documents/home.md"), style = "max-width:800px;"),
                        ),
               tabPanel(strong("Fit"), br(),
                        ####----------------------------------------------------------------
                        # Sidebar with inputs for plot parameters (to be included in tab)
                        sidebarLayout(
                          sidebarPanel( width = 4, br(),
                            
                            tabsetPanel(
                              tabPanel(strong("Upload"), br(),
                                       fluidRow(
                                         column(width=8, 
                                                 fileInput("file1", "Choose CSV File",multiple = FALSE),
                                                 ),
                                         column(width=4, actionButton("remove_file", "Remove", icon = icon("trash")))
                                       ),
                                       checkboxInput('header', 'Header', TRUE),
                                       # Plot Uploaded data
                                       actionButton("add_graph0", "Plot", icon = icon("paint-brush")),
                                       
                                       # Horizontal line ----
                                       tags$hr(),
                                       p("Please upload a", em(".csv") ,"file with drug-response data following the format described in",  
                                       strong("Help,"), "
                                        or download the example data set found in",
                                       downloadLink('downloadData2', strong(' here.') )
                                         )
                                        
                              ),
                              tabPanel(strong("Options"), br(),
                                       p("By default Nonparametric spline is selected, select other round checkmarks for other models. 
                                       Note that the nonparametric Bayesian model takes around 10 seconds to fit.
                                         The checkboxes add layers to the spline plot. The switches control the processing of the input data. 
                                         A complete explanation of each option is found in ",strong('Help.')),
                                       fluidRow(column(10,
                                                       # Check box for Sigmoid Fitting
                                                       strong('Nonparametric spline model'),
                                                       prettyCheckbox(inputId = "SplinePlot",
                                                                      label = "Nonparametric spline",
                                                                      value = TRUE, 
                                                                      bigger = TRUE, 
                                                                      shape = 'round',
                                                                      thick = TRUE,
                                                                      status = 'primary',
                                                                      animation = 'tada',
                                                                      icon = icon('check')
                                                       ),
                                                       prettyCheckboxGroup(
                                                         inputId = "checkgroup1",
                                                         label = "Plot Layers",
                                                         choices = c("Point samples",
                                                                     "Spline fit",
                                                                     "Min-max bands",
                                                                     "Empirical viability bands",
                                                                     "Drug span gradient",
                                                                     "Absolute doses",
                                                                     "Relative doses"                                                                     ),
                                                         icon = icon("check"),
                                                         animation = "tada",
                                                         status = "info",
                                                         selected = c("Point samples", "Spline fit")
                                                       ),
                                                       strong('Extra Models'),
                                                       # Check box for Monotone Fit
                                                       prettyCheckbox(inputId = "MonotonePlot",
                                                                      label = "Nonparamatric monotonic",
                                                                      value = FALSE, 
                                                                      bigger = TRUE, 
                                                                      shape = 'round',
                                                                      thick = TRUE,
                                                                      status = 'success',
                                                                      animation = 'tada',
                                                                      icon = icon('check')
                                                       ),
                                                       # Check box for Sigmoid Fitting
                                                       prettyCheckbox(inputId = "SigmoidPlot",
                                                                      label = "Parametric logistic",
                                                                      value = FALSE, 
                                                                      bigger = TRUE, 
                                                                      shape = 'round',
                                                                      thick = TRUE,
                                                                      status = 'danger',
                                                                      animation = 'tada',
                                                                      icon = icon('check')
                                                       ),
                                                       # Check box for NPB Fitting
                                                       prettyCheckbox(inputId = "NPBPlot",
                                                                      label = "Nonparametric Bayesian",
                                                                      value = FALSE, 
                                                                      bigger = TRUE, 
                                                                      shape = 'round',
                                                                      thick = TRUE,
                                                                      status = 'warning',
                                                                      animation = 'tada',
                                                                      icon = icon('check')
                                                       ),
                                                       strong('Data Processing'),
                                                       materialSwitch(inputId = "mean_switch", 
                                                                      label = "Median / Mean",
                                                                      status = "info",
                                                                      right=TRUE,
                                                                      value = TRUE ),
                                                       materialSwitch(inputId = "outlier_switch", 
                                                                      label = "Outliers Kept",
                                                                      status = "info",
                                                                      right=TRUE,
                                                                      value = TRUE ),
                                                       materialSwitch(inputId = "onehunda_switch", 
                                                                      label = "Viability over 100",
                                                                      status = "info",
                                                                      right=TRUE,
                                                                      value = TRUE ),
                                                       materialSwitch(inputId = "dosedep_auc", 
                                                                      label = "Dose dependent AUC",
                                                                      status = "danger",
                                                                      right=TRUE,
                                                                      value = TRUE )
                                       )
                                       )
                                       ),
                            
                              tabPanel(strong("Download"), br(),
                                       textInput(inputId  = "file_name", label = "File name"),
                                       radioButtons(inputId  = "file_type", label = "Select format for plot downloaded.",
                                                    choices = c("pdf", "png", "jpeg")),
                                       em("NB! If 'pdf' format was selected, the resolution setting will not work."),
                                       numericInput(inputId = "resolution", label = "Resolution (DPI)", value = 300, min = 30, max = 1000, step = 10),
                                       downloadButton(outputId = "download", label = "Download")
                                       
                              )
                              
                            ),
              

                          ),
                          # Show a plots of the generated distributions
                          mainPanel(
                            fluidRow(
                              column( 12,
                                     plotOutput("Plot1")
                                      ),
                              column( 12, 
                                      plotOutput("Plot2")
                                      ),
                              column( 12, 
                                      plotOutput("Plot3")
                              ),
                              column( 12, 
                                      plotOutput("Plot3.1")
                              ),
                            )
                          ))),
               tabPanel(strong("Help"),
                        tabsetPanel(
                          tabPanel(
                            strong("Data Instuction and Plot Options"),
                            br(),
                            tags$div(includeMarkdown("documents/help1.md"), style = "max-width:800px;"),
                            p("An Example dataset can be found in ",
                              downloadLink('downloadData', strong('here ')) ),
                            tags$div(includeMarkdown("documents/help2.md"), style = "max-width:800px;")    
                          ),
                          tabPanel(
                            strong("Model Specifications"),
                            tabsetPanel(
                              tabPanel(
                                strong('Nonparametric Spline'),
                                br(),
                                tags$div(includeMarkdown("documents/models1.md"), style = "max-width:800px;")
                              ),
                              tabPanel(
                                strong('Parametric Logistic'),
                                br(),
                                tags$div(includeMarkdown("documents/models2.md"), style = "max-width:800px;")
                              ),
                              tabPanel(
                                strong('Nonparametric Monotonic'),
                                br(),
                                tags$div(includeMarkdown("documents/models3.md"), style = "max-width:800px;")
                              ),
                              tabPanel(
                                strong('Nonparametric Bayesian'),
                                br(),
                                tags$div(includeMarkdown("documents/models4.md"), style = "max-width:800px;")
                              )
                            ),
                            
                          ),
                          tabPanel(
                            strong("Workshop"),
                            br(),
                            sidebarLayout(
                              sidebarPanel( width = 4, br(),
                                            titlePanel("Desired Drug-Patient Characteristics"),
                                            p('Drugs, patients, treatments and samples are options from the data found in the paper 
                                              "Intra-tumour diversification in colorectal cancer at the single-cell level" (2018), selecting different options will 
                                              result in a different model fit, note that not all combinations exist in the data.' ),
                                            prettyRadioButtons(
                                              inputId = "Drug",
                                              label = "Select Drug:",
                                              choices = c('5FU', 'afatinib', 'akt', 'doxorubicin', 'irinotecan', 'MEK12', 'nutlin3a'),
                                              selected = "5FU",
                                              animation = "pulse",
                                              status = "info",
                                              fill=TRUE
                                            ),
                                            # Select which Division(s) to plot
                                            prettyRadioButtons(inputId = "Patient",
                                                               label = "Select Patient:",
                                                               choices = c("P1", "P2", "P3"),
                                                               selected = "P1",
                                                               animation = "pulse",
                                                               status = "info",
                                                               fill = TRUE,
                                                               inline = TRUE
                                            ),
                                            prettyRadioButtons(inputId = "Treatment",
                                                               label = "Select Treatment:",
                                                               choices = c("T1",'T2','T3','T4','T5'),
                                                               selected = "T1",
                                                               animation = "pulse",
                                                               status = "info",
                                                               fill = TRUE,
                                                               inline = TRUE),
                                            prettyRadioButtons(inputId = "Sample",
                                                               label = "Select Sample:",
                                                               choices = 1:5,
                                                               selected = 1,
                                                               animation = "pulse",
                                                               status = "info",
                                                               fill = TRUE,
                                                               inline = TRUE),
                                            tags$hr(),
                                            actionButton("add_graph1", "Plot", icon = icon("paint-brush")),
                                            tags$hr(),
                                            strong('Nonparametric spline model'),
                                            prettyCheckbox(inputId = "SplinePlot_",
                                                           label = "Nonparametric spline",
                                                           value = TRUE, 
                                                           bigger = TRUE, 
                                                           shape = 'round',
                                                           thick = TRUE,
                                                           status = 'primary',
                                                           animation = 'tada',
                                                           icon = icon('check')
                                            ),
                                            prettyCheckboxGroup(
                                              inputId = "checkgroup1_",
                                              label = "Plot Layers",
                                              choices = c("Point samples",
                                                          "Spline fit",
                                                          "Min-max bands",
                                                          "Empirical viability bands",
                                                          "Drug span gradient",
                                                          "Absolute doses",
                                                          "Relative doses"                                                          ),
                                              icon = icon("check"),
                                              animation = "tada",
                                              status = "info",
                                              selected = c("Point samples","Spline fit")
                                            ),
                                            # Check box for Monotone Fit
                                            strong('Extra Models'),
                                            prettyCheckbox(inputId = "MonotonePlot_",
                                                           label = "Nonparametric monotonic",
                                                           value = FALSE, 
                                                           bigger = TRUE, 
                                                           shape = 'round',
                                                           thick = TRUE,
                                                           status = 'success',
                                                           animation = 'tada',
                                                           icon = icon('check')
                                            ),
                                            # Check box for Sigmoid Fitting
                                            prettyCheckbox(inputId = "SigmoidPlot_",
                                                           label = "Parametric logistic",
                                                           value = FALSE, 
                                                           bigger = TRUE, 
                                                           shape = 'round',
                                                           thick = TRUE,
                                                           status = 'danger',
                                                           animation = 'tada',
                                                           icon = icon('check')
                                            ),
                                            # Check box for Nonparametric Bayesian Fitting
                                            prettyCheckbox(inputId = "NPBPlot_",
                                                           label = "Nonparametric Bayesian",
                                                           value = FALSE, 
                                                           bigger = TRUE, 
                                                           shape = 'round',
                                                           thick = TRUE,
                                                           status = 'warning',
                                                           animation = 'tada',
                                                           icon = icon('check')
                                            ),
                                            strong('Data Processing'),
                                            materialSwitch(inputId = "mean_switch_", 
                                                           label = "Median / Mean",
                                                           status = "info",
                                                           right=TRUE,
                                                           value = TRUE ),
                                            materialSwitch(inputId = "outlier_switch_", 
                                                           label = "Outliers Kept",
                                                           status = "info",
                                                           right=TRUE,
                                                           value = TRUE ),
                                            materialSwitch(inputId = "onehunda_switch_", 
                                                           label = "Viability over 100",
                                                           status = "info",
                                                           right=TRUE,
                                                           value = TRUE ),
                                            materialSwitch(inputId = "dosedep_auc_", 
                                                           label = "Dose dependent AUC",
                                                           status = "danger",
                                                           right=TRUE,
                                                           value = TRUE )
                              ),
                              mainPanel(
                                fluidRow(
                                  column( 12,
                                          plotOutput("Plot4")
                                  ),
                                  column( 12, 
                                          plotOutput("Plot5")
                                  ),
                                  column( 12, 
                                          plotOutput("Plot6")
                                  ),
                                  column( 12, 
                                          plotOutput("Plot6.1")
                                  ),
                                )
                              )
                            )
                          )
                        )
               ),
               tabPanel(strong("FAQ"), br(),
                        tags$div(includeMarkdown("documents/questions.md"), style = "max-width:800px;")
               ),
               tabPanel(strong("Contact"), br(),
                        tags$div(includeMarkdown("documents/contact.md"), style = "max-width:800px;")
               )
             )
    )
        
)


server <- function(input, output) {
    
    cb0 = reactive(input$MonotonePlot)
    cb1 = reactive(input$SigmoidPlot)
    cb2 = reactive(input$SplinePlot)
    cb3 = reactive(input$NPBPlot)
    cbs = reactive(input$checkgroup1)
    ms = reactive(input$mean_switch)
    os = reactive(input$outlier_switch)
    hs = reactive(input$onehunda_switch)
    dd = reactive(input$dosedep_auc)
    
    cb0_ = reactive(input$MonotonePlot_)
    cb1_ = reactive(input$SigmoidPlot_)
    cb2_ = reactive(input$SplinePlot_)
    cb3_ = reactive(input$NPBPlot_)
    cbs_ = reactive(input$checkgroup1_)
    ms_ = reactive(input$mean_switch_)
    os_ = reactive(input$outlier_switch_)
    hs_ = reactive(input$onehunda_switch_)
    dd_ = reactive(input$dosedep_auc_)
    
    
    d1 = reactive(input$Drug)
    p1 = reactive(input$Patient)
    t1 = reactive(input$Treatment)
    s1 =  reactive(input$Sample)
    kp = reactive(input$Plots)
    
    fn = reactive(input$file_name)
    ft = reactive(input$file_type)
    rs = reactive(input$resolution)
    
    
    
    mydata <- reactive({
      
      inFile <- input$file1
      
      if (is.null(inFile))
        return(NULL)
      
      tbl <- read.csv(inFile$datapath, header=input$header)
      
      return(tbl)
    })
    
    
    
    observeEvent(input$add_graph1, {
    
      output$Plot4 <- renderPlot({
        showNotification("Generating plot...", duration = NULL, id = "message")
        
          drug = d1()
          patient = p1()
          treatment = t1()
          samp = s1()
          key_plot = kp()
          
          check_boxes  = cbs_()
          check_box0 = cb0_()
          check_box1 = cb1_()
          check_box2 = cb2_()
          check_box3 = cb3_()
          mean_switch = ms_()
          outlier_switch = os_()
          onehunda_switch = hs_()
          dosedependent_auc = dd_()
          
          
          p1 = NULL
          if( check_box2 ){
            # Note that it does not exist for all parameter combinations
            # add safepoint to function such that it returns null if its not present
            block = extract_dose_block(df_list, drug, patient, treatment, samp)
            if(is.null(block)){
              showNotification('Please select other drug-patient characteristics')
              return('')
            }
            block2 = preprocess_data(block, mean_switch, outlier_switch, onehunda_switch )
            p1 = PlotOverlay(block2, check_boxes, dosedependent_auc)
          }
          if( (!check_box2) & check_box0){
            block = extract_dose_block(df_list, drug, patient, treatment, samp)
            block2 = preprocess_data(block, mean_switch, outlier_switch, onehunda_switch )
            if(is.null(block)){
              showNotification('Please select other drug-patient characteristics')
              return('')
            }
            p1 = plot_monotoneFit(block2, dosedependent_auc )
          }
          if( (!check_box2) & (!check_box0) & check_box1 ){
            block = extract_dose_block(df_list, drug, patient, treatment, samp)
            block2 = preprocess_data(block, mean_switch, outlier_switch, onehunda_switch )
            if(is.null(block)){
              showNotification('Please select other drug-patient characteristics')
              return('')
            }
            p1 = plot_sigmodiFit(block2, dosedependent_auc)
          }
          if( (!check_box2) & (!check_box0) & (!check_box1) &  check_box3 ){
            block = extract_dose_block(df_list, drug, patient, treatment, samp)
            block2 = preprocess_data(block, mean_switch, outlier_switch, onehunda_switch, drop_values=F )
            if(is.null(block)){
              showNotification('Please select other drug-patient characteristics')
              return('')
            }
            p1 = plot_npbFit(block2, dosedependent_auc)
          }
          
          p1
          
          })
      
      output$Plot5 <- renderPlot({
        drug = d1()
        patient = p1()
        treatment = t1()
        samp = s1()
        key_plot = kp()
        
        check_box0 = cb0_()
        check_box1 = cb1_()
        check_box2 = cb2_()
        check_box3 = cb3_()
        mean_switch = ms_()
        outlier_switch = os_()
        onehunda_switch = hs_()
        dosedependent_auc = dd_()
        
        p2 = NULL
        
        if(check_box2 & check_box0 ){
          block = extract_dose_block(df_list, drug, patient, treatment, samp)
          block2 = preprocess_data(block, mean_switch, outlier_switch, onehunda_switch )
          if(is.null(block)){
            showNotification('Please select other drug-patient characteristics')
            return('')
          }
          p2 = plot_monotoneFit(block2, dosedependent_auc)
        }
        if(check_box2 & !check_box0 & check_box1){
          block = extract_dose_block(df_list, drug, patient, treatment, samp)
          block2 = preprocess_data(block, mean_switch, outlier_switch, onehunda_switch )
          if(is.null(block)){
            showNotification('Please select other drug-patient characteristics')
            return('')
          }
          p2 = plot_sigmodiFit(block2, dosedependent_auc)
        }
        if(!check_box2 & check_box0 & check_box1){
          block = extract_dose_block(df_list, drug, patient, treatment, samp)
          block2 = preprocess_data(block, mean_switch, outlier_switch, onehunda_switch )
          if(is.null(block)){
            showNotification('Please select other drug-patient characteristics')
            return('')
          }
          p2 = plot_sigmodiFit(block2, dosedependent_auc)
        }
        if((!check_box2 & !check_box0 & check_box1 & check_box3) |
           (!check_box2 & check_box0 & !check_box1 & check_box3) |
           (check_box2 & !check_box0 & !check_box1 & check_box3) ) {
          block = extract_dose_block(df_list, drug, patient, treatment, samp)
          block2 = preprocess_data(block, mean_switch, outlier_switch, onehunda_switch , drop_values=F)
          if(is.null(block)){
            showNotification('Please select other drug-patient characteristics')
            return('')
          }
          p2 = plot_npbFit(block2, dosedependent_auc)
        }
        
        if(!is.null(p2)) p2
      })
      
        output$Plot6 <- renderPlot({
          drug = d1()
          patient = p1()
          treatment = t1()
          samp = s1()
          key_plot = kp()
          
          check_box0 = cb0_()
          check_box1 = cb1_()
          check_box2 = cb2_()
          check_box3 = cb3_()
          mean_switch = ms_()
          outlier_switch = os_()
          onehunda_switch = hs_()
          dosedependent_auc = dd_()
          
          p3 = NULL
          
          if(check_box0 & check_box1 & check_box2){
            block = extract_dose_block(df_list, drug, patient, treatment, samp)
            if(is.null(block)){
              showNotification('Please select other drug-patient characteristics')
              return('')
            }
            block2 = preprocess_data(block, mean_switch, outlier_switch, onehunda_switch )
            p3 = plot_sigmodiFit(block2, dosedependent_auc)
          }
          
          if( (!check_box2 & check_box0 & check_box1 & check_box3) |
              (check_box2 & !check_box0 & check_box1 & check_box3) |
              (check_box2 & check_box0 & !check_box1 & check_box3)
              ){
            block = extract_dose_block(df_list, drug, patient, treatment, samp)
            if(is.null(block)){
              showNotification('Please select other drug-patient characteristics')
              return('')
            }
            block2 = preprocess_data(block, mean_switch, outlier_switch, onehunda_switch, drop_values=F )
            p3 = plot_npbFit(block2, dosedependent_auc)
          }
          
          if(!is.null(p3)) p3
      })
      
        output$Plot6.1 <- renderPlot({
          
          drug = d1()
          patient = p1()
          treatment = t1()
          samp = s1()
          key_plot = kp()
          
          check_box0 = cb0_()
          check_box1 = cb1_()
          check_box2 = cb2_()
          check_box3 = cb3_()
          mean_switch = ms_()
          outlier_switch = os_()
          onehunda_switch = hs_()
          dosedependent_auc = dd_()
          
          p4=NULL
          if(check_box0 & check_box1 & check_box2 & check_box3){
            block = extract_dose_block(df_list, drug, patient, treatment, samp)
            if(is.null(block)){
              showNotification('Please select other drug-patient characteristics')
              return('')
            }
            block2 = preprocess_data(block, mean_switch, outlier_switch, onehunda_switch, drop_values=F )
            
            p4 = plot_npbFit(block2, dosedependent_auc)
            p4
          }
          showNotification("Plot generated", duration = 1, id = "message")
          if(!is.null(p4)) p4
          
        })  
        
        
      # For downloading once plot is generated
      output$download <- downloadHandler(
        
        filename =  function() {
          if (input$file_name == ""){
            file <- "ENDS"
          } else {
            file <- input$file_name
          }
          return(paste(file, input$file_type, sep="."))
        },
        
        content = function(file){
          showNotification("Downloading plot please wait...", duration = NULL, id = "message")
     
          
          # Correct block 2 for input signals when added
          block = extract_dose_block(df_list, input$Drug, input$Patient, input$Treatment, input$Sample)
          block2 = preprocess_data(block)
          
          p1 = PlotOverlay(block2, input$checkgroup1)
          
          p2 = plot_sigmodiFit(block2)
          p3 = plot_monotoneFit(block2)
          if(input$NPBPlot){
            block2 = preprocess_data(block, drop_values=F)
            p4 = plot_npbFit(block2)  
          }
          
          
          wid  = 6*4
          hei = 4*4
          
          # This can be inside auxiliary function
          # Only 1 plot
          if(input$SplinePlot & !input$SigmoidPlot & !input$MonotonePlot & !input$NPBPlot){
            p = gridExtra::grid.arrange(p1, ncol=1, widths = wid, heights = hei)
            ggsave(file,plot = p,device = input$file_type,dpi = input$resolution,width = wid, height = hei, unit = 'cm' )
          }
          if(!input$SplinePlot & input$SigmoidPlot & !input$MonotonePlot & !input$NPBPlot){
            p = gridExtra::grid.arrange(p2, ncol=1, widths = wid, heights = hei)
            ggsave(file,plot = p,device = input$file_type,dpi = input$resolution,width = wid, height = hei, unit = 'cm' )
          }
          if(!input$SplinePlot & !input$SigmoidPlot & input$MonotonePlot & !input$NPBPlot){
            p = gridExtra::grid.arrange(p3, ncol=1, widths = wid, heights = hei)
            ggsave(file,plot = p,device = input$file_type,dpi = input$resolution,width = wid, height = hei, unit = 'cm' )
          }
          if(!input$SplinePlot & !input$SigmoidPlot & !input$MonotonePlot & input$NPBPlot){
            p = gridExtra::grid.arrange(p4, ncol=1, widths = wid, heights = hei)
            ggsave(file,plot = p,device = input$file_type,dpi = input$resolution,width = wid, height = hei, unit = 'cm' )
          }
          
           # 2 plots
          if(input$SplinePlot & input$SigmoidPlot & !input$MonotonePlot & !input$NPBPlot){
              p = gridExtra::grid.arrange(p1,  p2, ncol=1, nrow=2,  widths = wid, heights = c(hei,hei))
            ggsave(file,plot = p,device = input$file_type,dpi = input$resolution,width = wid, height = hei*2, unit = 'cm' )
          }
            
          if(input$SplinePlot & !input$SigmoidPlot & input$MonotonePlot & !input$NPBPlot){
            p = gridExtra::grid.arrange( p1, p3, ncol=1, nrow=2,  widths = wid, heights = c(hei,hei))
            ggsave(file,plot = p,device = input$file_type,dpi = input$resolution,width = wid, height = hei*2, unit = 'cm' )
          }
          
          if(!input$SplinePlot & input$SigmoidPlot & input$MonotonePlot & !input$NPBPlot){
            p = gridExtra::grid.arrange(p2,  p3, ncol=1, nrow=2,  widths = wid, heights = c(hei,hei))
            ggsave(file,plot = p,device = input$file_type,dpi = input$resolution,width = wid, height = hei*2, unit = 'cm' )
          }
          if(input$SplinePlot & !input$SigmoidPlot & !input$MonotonePlot & input$NPBPlot){
            p = gridExtra::grid.arrange(p1,  p4, ncol=1, nrow=2,  widths = wid, heights = c(hei,hei))
            ggsave(file,plot = p,device = input$file_type,dpi = input$resolution,width = wid, height = hei*2, unit = 'cm' )
          }
          if(!input$SplinePlot & input$SigmoidPlot & !input$MonotonePlot & input$NPBPlot){
            p = gridExtra::grid.arrange(p2,  p4, ncol=1, nrow=2,  widths = wid, heights = c(hei,hei))
            ggsave(file,plot = p,device = input$file_type,dpi = input$resolution,width = wid, height = hei*2, unit = 'cm' )
          }
          if(!input$SplinePlot & !input$SigmoidPlot & input$MonotonePlot & input$NPBPlot){
            p = gridExtra::grid.arrange(p3,  p4, ncol=1, nrow=2,  widths = wid, heights = c(hei,hei))
            ggsave(file,plot = p,device = input$file_type,dpi = input$resolution,width = wid, height = hei*2, unit = 'cm' )
          }
          
          
          # 3 plots
          if(input$SplinePlot & input$SigmoidPlot & input$MonotonePlot & !input$NPBPlot){
            p = gridExtra::grid.arrange(p1, p2, p3, ncol=1, nrow=3, widths = wid, heights = c(hei,hei,hei) )
            ggsave(file,plot = p,device = input$file_type,dpi = input$resolution,width = wid, height = hei*3, unit = 'cm' )
          }
          if(!input$SplinePlot & input$SigmoidPlot & input$MonotonePlot & input$NPBPlot){
            p = gridExtra::grid.arrange(p2, p3, p4, ncol=1, nrow=3, widths = wid, heights = c(hei,hei,hei) )
            ggsave(file,plot = p,device = input$file_type,dpi = input$resolution,width = wid, height = hei*3, unit = 'cm' )
          }
          if(input$SplinePlot & !input$SigmoidPlot & input$MonotonePlot & input$NPBPlot){
            p = gridExtra::grid.arrange(p1, p3, p4, ncol=1, nrow=3, widths = wid, heights = c(hei,hei,hei) )
            ggsave(file,plot = p,device = input$file_type,dpi = input$resolution,width = wid, height = hei*3, unit = 'cm' )
          }
          if(input$SplinePlot & input$SigmoidPlot & !input$MonotonePlot & input$NPBPlot){
            p = gridExtra::grid.arrange(p1, p2, p4, ncol=1, nrow=3, widths = wid, heights = c(hei,hei,hei) )
            ggsave(file,plot = p,device = input$file_type,dpi = input$resolution,width = wid, height = hei*3, unit = 'cm' )
          }
          
          # 4 plots
          if(input$SplinePlot & input$SigmoidPlot & !input$MonotonePlot & input$NPBPlot){
            p = gridExtra::grid.arrange(p1, p2, p3, p4, ncol=1, nrow=4, widths = wid, heights = c(hei,hei,hei,hei) )
            ggsave(file,plot = p,device = input$file_type,dpi = input$resolution,width = wid, height = hei*4, unit = 'cm' )
          }
          
          removeNotification(id = "message")
        }
      )
      
        }
    )
    
    observeEvent(input$add_graph0, {
      
      output$Plot1 <- renderPlot({
        
        showNotification("Generating plot...", duration = NULL, id = "message")
        check_boxes  = cbs()
        
        check_box0 = cb0()
        check_box1 = cb1()
        check_box2 = cb2()
        check_box3 = cb3()
        
        mean_switch = ms()
        outlier_switch = os()
        onehunda_switch = hs()
        dosedependent_auc = dd()
        
        
        block = mydata()
        
        p1 = NULL
        if(is.null(block)){
          showNotification('Please upload a file')
          return('')
        }
        block2 = preprocess_data(block, mean_samples = mean_switch, keep_outliers = outlier_switch, over_viability = onehunda_switch)
        
        if( check_box2 ){
          p1 = PlotOverlay(block2, check_boxes, dosedependent_auc)
        }
        if( (!check_box2) & check_box0){
          p1 = plot_monotoneFit(block2, dosedependent_auc)
        }
        if( (!check_box2) & (!check_box0) & check_box1 ){
          p1 = plot_sigmodiFit(block2, dosedependent_auc)
        }
        if( (!check_box2) & (!check_box0) & (!check_box1) & check_box3 ){
          block2 = preprocess_data(block, mean_samples = mean_switch, keep_outliers = outlier_switch, over_viability = onehunda_switch, drop_values=F)
          p1 = plot_npbFit(block2, dosedependent_auc)
        }
        
        p1
      })
      
      output$Plot2 <- renderPlot({
        
        
        check_box0 = cb0()
        check_box1 = cb1()
        check_box2 = cb2()
        check_box3 = cb3()
        
        mean_switch = ms()
        outlier_switch = os()
        onehunda_switch = hs()
        dosedependent_auc = dd()
        
        block = mydata()
        
        p2 = NULL
        
        if(is.null(block)){
          showNotification('Please upload a file')
          return('')
        }
        block2 = preprocess_data(block, mean_samples = mean_switch, keep_outliers = outlier_switch, over_viability = onehunda_switch)
        
        
        if(check_box2 & check_box0 ){
          p2 = plot_monotoneFit(block2, dosedependent_auc)
        }
        if(check_box2 & !check_box0 & check_box1){
          p2 = plot_sigmodiFit(block2, dosedependent_auc)
        }
        if((!check_box2 & !check_box0 & check_box1 & check_box3) |
           (!check_box2 & check_box0 & !check_box1 & check_box3) |
           (check_box2 & !check_box0 & !check_box1 & check_box3)){
          block2 = preprocess_data(block, mean_samples = mean_switch, keep_outliers = outlier_switch, over_viability = onehunda_switch, drop_values=F)
          p2 = plot_npbFit(block2, dosedependent_auc)
        }
        
        if(!is.null(p2)) p2
        
      })
      
      output$Plot3 <- renderPlot({
        
        check_box0 = cb0()
        check_box1 = cb1()
        check_box2 = cb2()
        check_box3 = cb3()
        mean_switch = ms()
        outlier_switch = os()
        onehunda_switch = hs()
        dosedependent_auc = dd()
        
        block = mydata()
        
        p3 = NULL
        
        if(is.null(block)){
          showNotification('Please upload a file')
          return('')
        }
        block2 = preprocess_data(block, mean_samples = mean_switch, keep_outliers = outlier_switch, over_viability = onehunda_switch)
        
        if(check_box2 & check_box0 & check_box1){
          p3 = plot_sigmodiFit(block2, dosedependent_auc)
        }
        if( (!check_box2 & check_box0 & check_box1 & check_box3) |
            (check_box2 & !check_box0 & check_box1 & check_box3) |
            (check_box2 & check_box0 & !check_box1 & check_box3) ){
          block2 = preprocess_data(block, mean_samples = mean_switch, keep_outliers = outlier_switch, over_viability = onehunda_switch, drop_values=F)
          p3 = plot_npbFit(block2, dosedependent_auc)
        }
        p3
        
      })
      
      output$Plot3.1 <- renderPlot({
        
        check_box0 = cb0()
        check_box1 = cb1()
        check_box2 = cb2()
        check_box3 = cb3()
        mean_switch = ms()
        outlier_switch = os()
        onehunda_switch = hs()
        dosedependent_auc = dd()
        
        block = mydata()
        
        p4 = NULL
        
        if(is.null(block)){
          showNotification('Please upload a file')
          return('')
        }
        block2 = preprocess_data(block, mean_samples = mean_switch, keep_outliers = outlier_switch, over_viability = onehunda_switch, drop_values=F)
        
        if(check_box2 & check_box0 & check_box1 & check_box3){
          p4 = plot_npbFit(block2, dosedependent_auc)
        }
        showNotification("Plot generated", duration = 1, id = "message")
        p4
        
      })
      
      # For downloading once plot is generated
      output$download <- downloadHandler(
        
        filename =  function() {
          if (input$file_name == ""){
            file <- "ENDS"
          } else {
            file <- input$file_name
          }
          return(paste(file, input$file_type, sep="."))
        },
        
        content = function(file){
          showNotification("Downloading plot...", duration = NULL, id = "message")
          
          # Note that it does not exist for all parameter combinations
          block = mydata()
          block2 = preprocess_data(block, mean_samples = input$mean_switch, keep_outliers = input$outlier_switch, over_viability = input$onehunda_switch)
          
          
          p1 = PlotOverlay(block2, input$checkgroup1, input$dosedep_auc)
          p2 = plot_monotoneFit(block2, input$dosedep_auc)
          p3 = plot_sigmodiFit(block2, input$dosedep_auc)
          if(input$NPBPlot){
            block2 = preprocess_data(block, mean_samples = input$mean_switch, keep_outliers = input$outlier_switch, over_viability = input$onehunda_switch, drop_values=F)
            p4 = plot_npbFit(block2, input$dosedep_auc)
          }
          
          wid  = 6*4
          hei = 4*4
          
          # This can be inside auxiliary function
          # Only 1 plot
          if(input$SplinePlot & !input$SigmoidPlot & !input$MonotonePlot & !input$NPBPlot){
            p = gridExtra::grid.arrange(p1, ncol=1, widths = wid, heights = hei)
            ggsave(file,plot = p,device = input$file_type,dpi = input$resolution,width = wid, height = hei, unit = 'cm' )
          }
          if(!input$SplinePlot & input$SigmoidPlot & !input$MonotonePlot & !input$NPBPlot){
            p = gridExtra::grid.arrange(p2, ncol=1, widths = wid, heights = hei)
            ggsave(file,plot = p,device = input$file_type,dpi = input$resolution,width = wid, height = hei, unit = 'cm' )
          }
          if(!input$SplinePlot & !input$SigmoidPlot & input$MonotonePlot & !input$NPBPlot){
            p = gridExtra::grid.arrange(p3, ncol=1, widths = wid, heights = hei)
            ggsave(file,plot = p,device = input$file_type,dpi = input$resolution,width = wid, height = hei, unit = 'cm' )
          }
          if(!input$SplinePlot & !input$SigmoidPlot & !input$MonotonePlot & input$NPBPlot){
            p = gridExtra::grid.arrange(p4, ncol=1, widths = wid, heights = hei)
            ggsave(file,plot = p,device = input$file_type,dpi = input$resolution,width = wid, height = hei, unit = 'cm' )
          }
          
          # 2 plots
          if(input$SplinePlot & input$SigmoidPlot & !input$MonotonePlot & !input$NPBPlot){
            p = gridExtra::grid.arrange(p1,  p2, ncol=1, nrow=2,  widths = wid, heights = c(hei,hei))
            ggsave(file,plot = p,device = input$file_type,dpi = input$resolution,width = wid, height = hei*2, unit = 'cm' )
          }
          
          if(input$SplinePlot & !input$SigmoidPlot & input$MonotonePlot & !input$NPBPlot){
            p = gridExtra::grid.arrange( p1, p3, ncol=1, nrow=2,  widths = wid, heights = c(hei,hei))
            ggsave(file,plot = p,device = input$file_type,dpi = input$resolution,width = wid, height = hei*2, unit = 'cm' )
          }
          
          if(!input$SplinePlot & input$SigmoidPlot & input$MonotonePlot & !input$NPBPlot){
            p = gridExtra::grid.arrange(p2,  p3, ncol=1, nrow=2,  widths = wid, heights = c(hei,hei))
            ggsave(file,plot = p,device = input$file_type,dpi = input$resolution,width = wid, height = hei*2, unit = 'cm' )
          }
          if(input$SplinePlot & !input$SigmoidPlot & !input$MonotonePlot & input$NPBPlot){
            p = gridExtra::grid.arrange(p1,  p4, ncol=1, nrow=2,  widths = wid, heights = c(hei,hei))
            ggsave(file,plot = p,device = input$file_type,dpi = input$resolution,width = wid, height = hei*2, unit = 'cm' )
          }
          if(!input$SplinePlot & input$SigmoidPlot & !input$MonotonePlot & input$NPBPlot){
            p = gridExtra::grid.arrange(p2,  p4, ncol=1, nrow=2,  widths = wid, heights = c(hei,hei))
            ggsave(file,plot = p,device = input$file_type,dpi = input$resolution,width = wid, height = hei*2, unit = 'cm' )
          }
          if(!input$SplinePlot & !input$SigmoidPlot & input$MonotonePlot & input$NPBPlot){
            p = gridExtra::grid.arrange(p3,  p4, ncol=1, nrow=2,  widths = wid, heights = c(hei,hei))
            ggsave(file,plot = p,device = input$file_type,dpi = input$resolution,width = wid, height = hei*2, unit = 'cm' )
          }
          
          
          # 3 plots
          if(input$SplinePlot & input$SigmoidPlot & input$MonotonePlot & !input$NPBPlot){
            p = gridExtra::grid.arrange(p1, p2, p3, ncol=1, nrow=3, widths = wid, heights = c(hei,hei,hei) )
            ggsave(file,plot = p,device = input$file_type,dpi = input$resolution,width = wid, height = hei*3, unit = 'cm' )
          }
          if(!input$SplinePlot & input$SigmoidPlot & input$MonotonePlot & input$NPBPlot){
            p = gridExtra::grid.arrange(p2, p3, p4, ncol=1, nrow=3, widths = wid, heights = c(hei,hei,hei) )
            ggsave(file,plot = p,device = input$file_type,dpi = input$resolution,width = wid, height = hei*3, unit = 'cm' )
          }
          if(input$SplinePlot & !input$SigmoidPlot & input$MonotonePlot & input$NPBPlot){
            p = gridExtra::grid.arrange(p1, p3, p4, ncol=1, nrow=3, widths = wid, heights = c(hei,hei,hei) )
            ggsave(file,plot = p,device = input$file_type,dpi = input$resolution,width = wid, height = hei*3, unit = 'cm' )
          }
          if(input$SplinePlot & input$SigmoidPlot & !input$MonotonePlot & input$NPBPlot){
            p = gridExtra::grid.arrange(p1, p2, p4, ncol=1, nrow=3, widths = wid, heights = c(hei,hei,hei) )
            ggsave(file,plot = p,device = input$file_type,dpi = input$resolution,width = wid, height = hei*3, unit = 'cm' )
          }
          
          # 4 plots
          if(input$SplinePlot & input$SigmoidPlot & !input$MonotonePlot & input$NPBPlot){
            p = gridExtra::grid.arrange(p1, p2, p3, p4, ncol=1, nrow=4, widths = wid, heights = c(hei,hei,hei,hei) )
            ggsave(file,plot = p,device = input$file_type,dpi = input$resolution,width = wid, height = hei*4, unit = 'cm' )
          }
          
          
          removeNotification(id = "message")
        }
      )
      
    }
    )
    
    output$downloadData <- downloadHandler(
      filename = function() {
        'ExampleData.csv'
      },
      content = function(con) {
        write.csv(df_example, con, row.names=FALSE)
      }
    )
    output$downloadData2 <- downloadHandler(
      filename = function() {
        'ExampleData.csv'
      },
      content = function(con) {
        write.csv(df_example, con, row.names=FALSE)
      }
    )
    
    
}

# Run the application 
shinyApp(ui = ui, server = server)
