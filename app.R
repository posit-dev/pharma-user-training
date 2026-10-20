# app.R -- TEAE explorer for study CDISCPILOT01
# A reviewer picks a severity level and a system organ class; the table and the
# Kaplan-Meier figure update to show only the adverse events that match.
#
# Run: open this file in Positron and click "Run App" (or shiny::runApp()).

library(shiny)
library(dplyr)
library(tidyr)
library(ggsurvfit)

arms <- c("Placebo", "Xanomeline Low Dose", "Xanomeline High Dose")  # protocol order

# Both files are limited to the safety population; ADAE holds treatment-emergent events only
adsl <- read.csv("data/ADSL.csv") |>
  mutate(TRT01A = factor(TRT01A, levels = arms),
         TRTSDT = as.Date(TRTSDT), TRTEDT = as.Date(TRTEDT), LSTALVDT = as.Date(LSTALVDT))
adae <- read.csv("data/ADAE.csv") |>
  mutate(TRT01A = factor(TRT01A, levels = arms), ASTDT = as.Date(ASTDT))

arm_n <- count(adsl, TRT01A, name = "N")
soc_choices <- sort(unique(adae$AEBODSYS))

ui <- fluidPage(
  titlePanel("CDISCPILOT01 Safety Explorer"),
  sidebarLayout(
    sidebarPanel(
      checkboxGroupInput("severity", "Severity",
                         choices = c("MILD", "MODERATE", "SEVERE"),
                         selected = c("MILD", "MODERATE", "SEVERE")),
      selectInput("soc", "System organ class", choices = c("All", soc_choices)),
      helpText("Subjects with at least one matching TEAE, safety population.")
    ),
    mainPanel(
      textOutput("summary"),
      h4("Most frequent preferred terms, n (%)"),
      tableOutput("teae_table"),
      h4("Time to first matching TEAE"),
      plotOutput("km_plot")
    )
  )
)

server <- function(input, output, session) {

  # Adverse events that match the filters. Everything below builds on this.
  adae_f <- reactive({
    d <- filter(adae, AESEV %in% input$severity)
    if (input$soc != "All") d <- filter(d, AEBODSYS == input$soc)
    d
  })

  output$summary <- renderText({
    n <- n_distinct(adae_f()$USUBJID)
    sprintf("%d of %d subjects (%.1f%%) have a matching TEAE.", n, nrow(adsl), 100 * n / nrow(adsl))
  })

  output$teae_table <- renderTable({
    req(nrow(adae_f()) > 0)
    top10 <- adae_f() |>
      distinct(USUBJID, AEDECOD) |>
      count(AEDECOD, sort = TRUE) |>
      slice_head(n = 10) |>
      pull(AEDECOD)

    adae_f() |>
      filter(AEDECOD %in% top10) |>
      distinct(USUBJID, TRT01A, AEDECOD) |>
      count(AEDECOD, TRT01A, name = "n") |>
      complete(AEDECOD, TRT01A, fill = list(n = 0L)) |>
      left_join(arm_n, by = "TRT01A") |>
      mutate(cell = sprintf("%d (%.1f%%)", n, 100 * n / N)) |>
      select(AEDECOD, TRT01A, cell) |>
      pivot_wider(names_from = TRT01A, values_from = cell) |>
      mutate(AEDECOD = factor(AEDECOD, levels = top10)) |>
      arrange(AEDECOD) |>
      rename(`Preferred term` = AEDECOD)
  })

  output$km_plot <- renderPlot({
    req(nrow(adae_f()) > 0)
    first_teae <- adae_f() |>
      group_by(USUBJID) |>
      summarise(teae_date = min(ASTDT), .groups = "drop")

    tte <- adsl |>
      left_join(first_teae, by = "USUBJID") |>
      mutate(
        event = as.integer(!is.na(teae_date)),
        end_date = if_else(!is.na(TRTEDT), TRTEDT + 30, LSTALVDT),   # censor 30 days after last dose
        time = pmax(1, as.numeric(if_else(event == 1, teae_date, end_date) - TRTSDT) + 1)
      )

    survfit2(Surv(time, event) ~ TRT01A, data = tte) |>
      ggsurvfit() +
      add_confidence_interval() +
      labs(x = "Days since first dose", y = "Probability of no matching TEAE")
  })
}

shinyApp(ui, server)
