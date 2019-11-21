#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

# Load necessary libraries
library(Matrix)
library(readxl)
library(deSolve)
library(pracma)
library(neldermead)
library(profvis)
library(Rcpp)
library(ggplot2)
library(gridExtra)

source("fun/get_addresses.R")
source("fun/get_distribution_fns.R")
source("fun/get_objective.R")
source("fun/goveqs_basis.R")
source("fun/scale_up.R")
source("fun/allocate_parameters.R")
source("fun/make_model.R")
source("fun/Make_distr_fns.R")
source("fun/return_output.R")
source("fun/plot_targets.R")
source("fun/get_input.R")
source("fun/get_interventions.R")
source("fun/make_intervention.R")


sourceCpp("fun/compute_dx.cpp")
sourceCpp("fun/scale_matrix.cpp")


library(shiny)


# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("DTO TB model"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
     sidebarPanel(
       
       helpText(),
       selectInput('x', 'Location:',
                   choices = c("Bihar","Kerala")),
     
     
       
       
       
       
         
       radioButtons('format', 'Document format', c('PDF', 'HTML', 'Word'),
                    inline = TRUE),+
       downloadButton('downloadReport')
     ),
      
      # Show a plot of the generated distribution
      mainPanel(
         plotOutput("targetPlot")
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  locations<-'Bihar_u'
  location<-"Bihar_u"
  
  action  <- "plot_mle"
  lsq<-0     # Use least squares or llk
  nmcmc<-1000
  nint<-5    #Number of inyterventions
  
  
  # for (q in 1:length(locations)){
  #   
  # location=locations{q}
  
  modein<-get_input(location)
  prm<-modein$parm ; ref<-modein$ref; 
  sel<-modein$sel  ; agg<-modein$agg;
  gps<-modein$gps
  datapoints<-modein$datapoints
  #__________________________________________________________________________
  #  Create instances of model functions
  #__________________________________________________________________________
  
  obj       <-function(x)  get_objective(x,prm, ref, sel, agg, gps, lhd,FALSE)
  obj_spx   <-function(x) -get_objective(x,prm, ref, sel, agg, gps, lhd)
  obj_mcmc  <-function(x)  get_objective(x,prm, ref, sel, agg, gps, lhd)
  
  
  
  
  f<-paste("res/","output","_",location,"_","mle",sep="")
     runs<-readRDS(f)
   
  
  
   
   output$targetPlot <- renderPlot({
     plot_targets(runs,datapoints,location)
      # # generate bins based on input$bins from ui.R
      # x    <- faithful[, 2] 
      # bins <- seq(min(x), max(x), length.out = input$bins + 1)
      # 
      # # draw the histogram with the specified number of bins
      # hist(x, breaks = bins, col = 'darkgray', border = 'white')
   })
}

# Run the application 
shinyApp(ui = ui, server = server)

