#' ConcentrationGames UI
#' 
#' Shiny user interface for CocentrationGames
#'
#' @import shinydashboard
ui = dashboardPage(
  dashboardHeader(title = "Concentration Games"),
  dashboardSidebar(disable = TRUE),
  dashboardBody(
    box(title = "Cross Section Geometry", width = 12, collapsible = TRUE, solidHeader = TRUE,
      box(title = "Left Overbank", width = 4, status = "warning", solidHeader = TRUE,
        sliderInput("LOB.height", "Height", min = 1, max = 100, value = 5),
        sliderInput("LOB.width", "Width", min = 1, max = 1000, value = 50)
      ),
      box(title = "Channel", width = 4, status = "primary", solidHeader = TRUE,
        sliderInput("C.height", "Height", min = 0, max = 0, value = 0),
        sliderInput("C.width", "Width", min = 1, max = 100, value = 10)      
      ),
      box(title = "Right Overbank", width = 4, status = "danger", solidHeader = TRUE,
        sliderInput("ROB.height", "Height", min = 1, max = 100, value = 8),
        sliderInput("ROB.width", "Width", min = 1, max = 1000, value = 50)      
      )
    ),
    fluidRow(
      tabBox(width = 4, selected = "Channel Settings",
        tabPanel("Channel Settings",
          numericInput("H", "Water Level", min = 1, value = 15),
          numericInput("mass", "Mass Deposited", min = 0, value = 1e5),
          numericInput("decay", "Overbank Decay Rate", min = 0, value = 0.01),
          numericInput("truncate.dist", "Truncation Distance", min = 1, value = 1000)
        ),
        tabPanel("Sediment Settings",
          selectInput("cfun", "Gradient Function", choices = c("Constant", "Rouse", 
            "Lane and Kalinske", "Barenblatt", "Tanaka and Sugimoto")),
          numericInput("b", "Active Layer Height", value = 0.05),
          numericInput("w", "Settling Velocity", value = 0.03),
          numericInput("us", "Shear Velocity", value = 0.2),
          numericInput("k", "von Karmen Constant", value = 0.4)
        )
      ),
      box(title = "Deposition result", status = "primary", width = 8,
         plotOutput("xschange")
      )
    )
  )
)

#' ConcentrationGames Server
#' 
#' Shiny Server for CocentrationGames
#'
#' @import ggplot2
server <- function(input, output, session) {
  session$onSessionEnded(stopApp)
  observeEvent(input$close, stopApp() )
  output$xschange = renderPlot({
    conc.fun = switch(input$cfun,
      "Constant" = NULL,
      "Rouse" = rouse,
      "Lane and Kalinske" = lane,
      "Barenblatt" = barenblatt,
      "Tanaka and Sugimoto" = tanaka
    )
    xs = define_xs(input$C.height, input$C.width, input$ROB.height, 
      input$ROB.width, input$LOB.height, input$LOB.width)
    res = distribute_mass(input$mass, xs, input$H, input$b, conc.fun, us = input$us, k = input$k, w = input$w)
    massdep = deposit_mass(res, 0.1, input$decay, input$truncate.dist)
    bedchange = mass_to_bedchange(massdep)
    ggplot(bedchange) + aes(x = y, y = z + dz) + geom_line(linetype = "dashed") +
      geom_line(aes(y = z), color = "brown") + geom_hline(yintercept = input$H, color = "blue")
  })
}

#' ConcentrationGames
#'
#' Open the Concentration Games app.
#'
#' @import shiny
#' @export
ConcentrationGames = function() {
  shinyApp(ui, server)
}
