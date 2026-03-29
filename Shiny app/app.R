library(shiny)
library(tidyverse)
library(ordbetareg)
library(truncnorm)
library(shinydashboard)
library(statConfR)

# Source simulation functions
source("simulation_functions.R")

# Define UI
ui <- dashboardPage(
  dashboardHeader(title = "Heuristic Model Simulator"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Model Simulator", tabName = "simulator", icon = icon("sliders")),
      menuItem("About Model", tabName = "about", icon = icon("info"))
    )
  ),
  
  dashboardBody(
    tabItems(
      # Tab 1: Simulator
      tabItem(tabName = "simulator",
        fluidRow(
          # Parameter panel
          column(3,
            box(title = "Model Parameters", status = "primary", solidHeader = TRUE, width = 12,
              
              # Basic parameters
              h4("Perceptual Noise", class = "section-header"),
              sliderInput("sigma_e", 
                label = "Evidence Noise (sigma_e) [log scale]", 
                min = -3, max = 0, value = -2, step = 0.1),
              
              h4("Choice Noise", class = "section-header"),  
              sliderInput("sigma_k", 
                label = "Criterion Variability (sigma_k) [log scale]", 
                min = -3, max = 0, value = -2, step = 0.1),
              
              h4("Metacognitive Noise", class = "section-header"),
              sliderInput("sigma_m", 
                label = "Metacognitive Noise (sigma_m) [log scale]", 
                min = -3, max = 0, value = -1, step = 0.1),
              
              h4("Other Parameters", class = "section-header"),
              sliderInput("bias", 
                label = "Stimulus Bias (beta)", 
                min = -1, max = 1, value = 0, step = 0.05),
              
              h4("Other Parameters", class = "section-header"),
              sliderInput("metabias", 
                          label = "Confidence Bias (beta)", 
                          min = -3, max = 3, value = 0, step = 0.05),
              
              
              h4("Simulation Settings", class = "section-header"),
              sliderInput("n_trials",
                label = "Number of Trials",
                min = 100, max = 2000, value = 500, step = 100),
              
              sliderInput("n_coherence_levels",
                label = "Number of Stimulus Levels",
                min = 3, max = 15, value = 7, step = 1),
              
              sliderInput("seed",
                label = "Random Seed",
                min = 1, max = 10000, value = 123, step = 1),
              
              actionButton("simulate_btn", "Simulate Data", 
                class = "btn-success", width = "100%"),
              
              actionButton("simulate_seed_btn", "Simulate From Seed", 
                class = "btn-info", width = "100%"),
              
              # Reset button
              actionButton("reset_btn", "Reset to Defaults", 
                class = "btn-default", width = "100%")
            )
          ),
          
          # Main visualization area
          column(9,
            # Psychometric function
            box(title = "Psychometric Function (P(Choice = 1))", 
              status = "success", solidHeader = TRUE, width = 12,
              plotOutput("psychometric_plot", height = "400px"),
              htmlOutput("psychometric_stats")
            ),
            
            # Confidence ratings
            box(title = "Mean Confidence Ratings by Accuracy", 
              status = "success", solidHeader = TRUE, width = 12,
              plotOutput("confidence_plot", height = "400px"),
              htmlOutput("confidence_stats")
            )
          )
        )
      ),
      
      # Tab 2: About the model
      tabItem(tabName = "about",
        fluidRow(
          column(12,
            box(title = "About the Heuristic Model of Metacognition", 
              status = "info", solidHeader = TRUE, width = 12,
              includeHTML("model_description.html")
            )
          )
        )
      )
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  
  # Store simulated data
  sim_data <- reactiveVal(NULL)
  
  # Simulate data when button is clicked (random)
  observeEvent(input$simulate_btn, {
    withProgress(message = "Simulating...", value = 0, {
      data <- simulate_data_mcmc_x1_custom(
        n = input$n_trials,
        sigma_e = input$sigma_e,
        sigma_k = input$sigma_k,
        sigma_m = input$sigma_m,
        bias = input$bias,
        metabias = input$metabias,
        n_coherence_levels = input$n_coherence_levels,
        seed = NULL
      )
      sim_data(data)
      incProgress(1)
    })
  }, ignoreNULL = TRUE)
  
  # Simulate data from seed when button is clicked
  observeEvent(input$simulate_seed_btn, {
    withProgress(message = "Simulating from seed...", value = 0, {
      data <- simulate_data_mcmc_x1_custom(
        n = input$n_trials,
        sigma_e = input$sigma_e,
        sigma_k = input$sigma_k,
        sigma_m = input$sigma_m,
        bias = input$bias,
        metabias = input$metabias,
        n_coherence_levels = input$n_coherence_levels,
        seed = input$seed
      )
      sim_data(data)
      incProgress(1)
    })
  }, ignoreNULL = TRUE)
  
  # Reset parameters to defaults
  observeEvent(input$reset_btn, {
    updateSliderInput(session, "sigma_e", value = -2)
    updateSliderInput(session, "sigma_k", value = -1)
    updateSliderInput(session, "sigma_m", value = -1)
    updateSliderInput(session, "bias", value = 0)
    updateSliderInput(session, "metabias", value = 0)
    updateSliderInput(session, "conf_prec", value = 50)
    updateSliderInput(session, "n_trials", value = 500)
    updateSliderInput(session, "seed", value = 123)
  })
  
  # Render psychometric function plot
  output$psychometric_plot <- renderPlot({
    if (is.null(sim_data())) return(NULL)
    
    df <- sim_data()
    
    # Bin the data by signed stimulus (X * D)
    plot_data <- df %>%
      mutate(
        XD = X * D,
        XD_binned = round(XD * 20) / 20,
        choice = a_gauss
      ) %>%
      group_by(XD_binned) %>%
      summarize(
        mean_choice = mean(choice, na.rm = TRUE),
        n = n(),
        se = sqrt(mean_choice * (1 - mean_choice) / n),
        .groups = "drop"
      ) %>%
      filter(!is.na(XD_binned))
    
    # Create plot
    plot_data %>%
      ggplot(aes(x = XD_binned, y = mean_choice)) +
      geom_point(size = 3, alpha = 0.6, color = "#333333") +
      geom_errorbar(aes(ymin = mean_choice - 1.96*se, 
                        ymax = mean_choice + 1.96*se),
                    width = 0.02, alpha = 0.5, color = "#333333") +
      geom_smooth(se = FALSE, alpha = 0.2, color = "#555555", fill = "#cccccc", method = "glm", 
                  method.args = list(family = "binomial")) +
      labs(x = "Signed Stimulus (X × D)", y = "P(Choice = 1)",
           title = "Psychometric Function") +
      geom_hline(yintercept = 0.5, linetype = 2)+
      geom_vline(xintercept = 0, linetype = 2)+
      theme_classic(base_size = 12) +
      scale_y_continuous(limits = c(-0.05, 1.05)) +
      scale_x_continuous(limits = c(-1, 1)) +
      theme(legend.position = "top")
  })
  
  # Render confidence plot
  output$confidence_plot <- renderPlot({
    if (is.null(sim_data())) return(NULL)
    
    df <- sim_data()
    
    # Prepare data for confidence plot
    conf_data <- df %>%
      mutate(
        XD = X * D,
        XD_binned = round(XD * 20) / 20,
        ACC = as.factor(ifelse(ACC == 1, "Correct", "Error"))
      ) %>%
      group_by(XD_binned, ACC) %>%
      summarize(
        mean_conf = mean(c_mu, na.rm = TRUE),
        sd_conf = sd(c_mu, na.rm = TRUE),
        n = n(),
        se = sd_conf / sqrt(n),
        .groups = "drop"
      ) %>%
      filter(!is.na(XD_binned))
    
    # Create plot
    conf_data %>%
      ggplot(aes(x = XD_binned, y = mean_conf,
                 color = ACC, fill = ACC)) +
      geom_point(size = 3, alpha = 0.6) +
      geom_errorbar(aes(ymin = mean_conf - 1.96*se,
                        ymax = mean_conf + 1.96*se),
                    width = 0.02, alpha = 0.5) +
      geom_smooth(method = "loess", se = FALSE, alpha = 0.2) +
      geom_hline(yintercept = 0.5, linetype = 2)+
      geom_vline(xintercept = 0, linetype = 2)+
      scale_color_manual(
        name = "Accuracy",
        values = c("Error" = "#d62728", "Correct" = "#2ca02c")
      ) +
      scale_fill_manual(
        name = "Accuracy",
        values = c("Error" = "#d62728", "Correct" = "#2ca02c")
      ) +
      labs(x = "Signed Stimulus (X × D)", y = "Mean Confidence Rating",
           title = "Confidence Ratings by Accuracy") +
      theme_classic(base_size = 12) +
      scale_y_continuous(limits = c(0, 1)) +
      scale_x_continuous(limits = c(-1, 1)) +
      theme(legend.position = "top")
  })
  
  # Psychometric function statistics
  output$psychometric_stats <- renderUI({
    if (is.null(sim_data())) return(NULL)
    
    df <- sim_data()
    
    overall_acc <- mean(df$ACC, na.rm = TRUE)
    
    HTML(paste(
      "<p><strong>Overall Accuracy:</strong>", round(overall_acc, 3),
      "</p>"
    ))
  })
  
  # Confidence statistics
  output$confidence_stats <- renderUI({
    if (is.null(sim_data())) return(NULL)
    
    df <- sim_data()
    
    data = df %>% mutate(conf_bin = cut(c_mu,4),
                         stimulus  = as.factor(D),
                         participant = 1,
                         rating = as.numeric(as.factor(conf_bin)),
                         correct = ACC) %>% drop_na()
    
    
    mratio = fitMetaDprime(data)

    mean_conf_correct <- mean(df$c_mu[df$ACC == 1], na.rm = TRUE)
    mean_conf_error <- mean(df$c_mu[df$ACC == 0], na.rm = TRUE)
    
    HTML(paste(
      "<p><strong>Mean Confidence (Correct):</strong>", round(mean_conf_correct, 3),
      "</p><p><strong>Mean Confidence (Error):</strong>", round(mean_conf_error, 3),
      "</p><p><strong> Sensitivity (d'):</strong>", (round(mratio$dprime,3)),
      "</p><p><strong>Metacognitive Sensitivity (meta d'):</strong>", (round(mratio$metaD,3)),
      "</p><p><strong>Metacognitive efficiency (M-ratio):</strong>", (round(mratio$Ratio,3)),
      "</p>"
    ))
  })
  

}

# Run the application
shinyApp(ui = ui, server = server)
