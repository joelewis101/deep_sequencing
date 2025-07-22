

library(tidyverse)
library(here)
library(janitor)
library(lubridate)

library(shiny)


# Define UI for dataset viewer app ----
ui <- fluidPage(

  # App title ----
  titlePanel("DASSIM stool DNA boxplans"),

  # Sidebar layout with a input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(

      # Input: Selector for choosing dataset ----

      # Input: Numeric entry for number of obs to view ----
      numericInput(inputId = "box",
                   label = "Box ID",
                   value = 1)
    ),

    # Main panel for displaying outputs ----
    mainPanel(
  plotOutput(outputId = "distPlot")
    )
  )
)

server <- function(input, output) {
  output$distPlot <- renderPlot({


source("load_and_clean_mainfest.R")
    df_out %>%
  left_join(filter(df_extraction, !grepl("NC$", sample_id))) %>%
  group_by(sample_id) %>%
  mutate(labz =
         case_when(grepl("NC", sample_id) ~ sample_id,
                   TRUE ~ paste(unique(extraction_id), collapse = ","))
         ) %>%
  filter(box == input$box) %>%
  ggplot(aes(column,
             fct_rev(row),
             fill = as.character(extraction_date),
             shape = !is.na(negative_control))) +
  geom_tile() +
  geom_point() +
  scale_x_continuous(breaks = 1:9) +
  labs(title = paste0("Box ", input$box),
       shape = "has linked neg control") +
  geom_text(aes(label = labz), nudge_y = 0.2, size = 3)

    })

}
shinyApp(ui,server)
