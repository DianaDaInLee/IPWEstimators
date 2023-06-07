rm(list = ls())
set.seed(890324)
library(tidyverse)
library(magrittr)
library(ggplot2)
library(data.table)
library(ggrepel)
library(gridExtra)
library(ggridges)
require(shiny)
require(shinydashboard)
require(shinyWidgets)
require(htmlwidgets)
require(plotly)
library(DT)
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source('fn.R')

ui <- dashboardPage(
  dashboardHeader(title = "Multi-arm Bandit Simulation", titleWidth = 300),
  dashboardSidebar(
    width = 300,
    sidebarMenu(id="menu1",
                sliderInput(inputId  = "K", label = "No. of arms:", min = 2, max = 50, value = 9)),
    sidebarMenu(id="menu3",
                sliderInput(inputId  = "period", label = "No. of periods:", min = 2, max = 100, value = 10)),
    sidebarMenu(id="menu2",
                textInput(inputId  = "N", label = "Total sample size:", value = '1000')),
    sidebarMenu(id="menu4",
                textInput(inputId  = "true", label = "True parameters:", value = "0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1"),
                column(12, helpText('Insert true success rates for all arms 1 thru K'))),
                # column(12, helpText('from 1 through K.'))),
sidebarMenu(id="menu5",
                textInput(inputId  = "floor", label = "Floor rate :", value = "0.01")),
    sidebarMenu(id="menu6",
                selectInput(inputId  = "type", label = "Sampling method:", 
                            choices  = c('Static', 'TS', 'TS, BW'),
                            selected = c('Static', 'TS', 'TS, BW'),
                            multiple = T)),
    sidebarMenu(id="menu5",
                sliderInput(inputId  = "iter", label = "Iterations:", min = 1, max = 50, value = 1)),
    br(),
    column(12, align = 'left', offset = 0, submitButton("Run simulation", icon = icon("play")))
  ),
  dashboardBody(
    tabItem("ml",
            fluidRow(
              tabBox(
                width = 20,
                tabPanel(title = "Results",
                         fluidRow(
                           box(width = 6,  status = "info", solidHeader = FALSE, 
                               title = "Simulated Posterior Probabilities over Period", 
                               'If iteration > 1, smooth curves of posterior probabilities will be fitted via Generalized Additive Model (GAM).', br(),
                               #'TS, Uniform = Thompson Sampling; TS, BW = Thompson Sampling with Balancing Weight', 
                               plotOutput("posterior")),
                           box(width = 6,  status = "info", solidHeader = FALSE, 
                               title = "Cumulative Sample Size over Period", 
                               'TS = Thompson Sampling; TS, BW = Thompson Sampling with Balancing Weight', 
                               plotOutput("nsize")),
                           box(width = 12,  status = "info", solidHeader = FALSE, 
                               title = "Distribution of Bias", 
                               'Points represent average bias, and error bars are 0.025 and 0.975 quantiles across all iterations. Error bars not shown if iteration = 1.',
                               plotlyOutput("bias")),
                           box(width = 12,  status = "info", solidHeader = FALSE, 
                               title = "Distribution of CI Radius", 
                               'Points represent average CI radius, and error bars are 0.025 and 0.975 quantiles across all iterations. Error bars not shown if iteration = 1.',
                               plotlyOutput("ci")),
                           box(width = 12, status = "info", solidHeader = FALSE,  
                               title = "Summary Results", 
                               'RMSE is average root mean squared error of the estimates of the mean of each arm, and CI Radius is average CI radius across all iterations.',
                               dataTableOutput("result"))
                         )
                )
              )
              )
    )
  )
)

server <- function(input, output){
  d <- reactive({
    print('# SET PARAMETERS') 
    S_s     <- 10000
    iter_s  <- as.numeric(input$iter)
    floor_s <- as.numeric(input$floor)
    periods_s    <- as.numeric(input$period)
    K_s          <- as.numeric(input$K)
    N_s          <- as.numeric(input$N)
    true_theta_s <- c(as.numeric(strsplit(input$true, ",")[[1]]))
    type_s <- input$type
    pcaption <- paste0("Simulated result (iterations = ", iter_s, ") under static (left), Thompson sampling (center) and Thompson sampling with balanced weights (right).\n", 
                       "Based on a sampling of ", N_s/periods_s, " observations per period, assigning treatment to ", K_s,  " arms.\n",
                       "In the adaptive algorithms, uniform priors are assumed for all arms, and arms are sampled with equal probability in the first period.\n",
                       "Success probabilies are: ", input$true,". Floor rate of ", floor_s, " is used.")
    
    print('# SIMULATE')
    sim  <- list()
    for (i in type_s) sim[[i]] <- sim_iter_equalsize(periods = periods_s, N = N_s, K = K_s, true_theta = true_theta_s, 
                                                     iter = iter_s, S = S_s, FLR = floor_s, type = i)
    
    print('# POSTERIOR PROBABILITY')
    post <- plyr::ldply(names(sim), function(x){
      sim[[x]]$posterior_full %>% mutate(type = x) %>% pivot_longer(cols = starts_with('T', ignore.case = F))
    })
    post$type   <- factor(post$type, levels = type_s)
    post$p.cat  <- as.factor(post$name)
    post_t <- post %>%
      filter(period == max(post$period)) %>%
      group_by(type, period, name, p.cat) %>%
      summarize(value = mean(value))
    
    if (input$iter == 1){
      p_posterior <- post %>%
        ggplot(aes(x = period, y = value, group = name, color = p.cat, linetype = p.cat)) +
        geom_line()
    }
    if (input$iter > 1){
      p_posterior <- post %>%
        ggplot(aes(x = period, y = value, group = name, color = p.cat, linetype = p.cat)) +
        geom_smooth(method = 'gam', linewidth = 0.5, alpha = 0.2)
    }
    ymin   <- floor(min(post$value)*10)/10 - 0.1
    ymax   <- ceiling(max(post$value)*10)/10
    ptheme <- theme_bw() +
              theme(legend.position  = 'none',
                    panel.grid.minor = element_blank(),
                    plot.caption = element_text(hjust = 0))
    
    p_posterior <- p_posterior +
      facet_wrap(vars(type), nrow = 1) +
      coord_cartesian(xlim = c(0, (periods_s + 4)),  ylim = c(ymin, ymax), clip = 'off') + 
      geom_text_repel(data = post_t, aes(label = p.cat), nudge_x = 4, hjust = 1, segment.size = .2, seed = 343, direction = 'y', size = 3) +
      scale_x_continuous(breaks = c(seq(1, periods_s, 1))) +
      scale_y_continuous(breaks = c(seq(0, ymax, 0.1))) +
      scale_colour_manual(name = 'Position in last period',
                          breaks = levels(post$p.cat),
                          labels = levels(post$p.cat),
                          values = paste0('gray', seq(1, 90, floor((90-1)/K_s)))) +
      ylab('Posterior probability of being the best arm') + xlab('Number of period') + labs(caption = pcaption) + ptheme
    
    print('# CUMULATIVE SAMPLE SIZE')
    nsize <- plyr::ldply(names(sim), function(x){
      sim[[x]]$n_full %>% mutate(type = x) %>% pivot_longer(cols = starts_with('T', ignore.case = F))
    })
    nsize$type   <- factor(nsize$type, levels = type_s)
    nsize$p.cat  <- as.factor(nsize$name)
    nsize_t <- nsize[nsize$period == max(periods_s), ] %>%
      group_by(type, period, name, p.cat) %>%
      summarize(value = mean(value))
    
    if (input$iter == 1){
      p_nsize <- nsize %>%
        ggplot(aes(x = period, y = value, group = name, color = p.cat, linetype = p.cat)) +
        geom_line()
    }
    if (input$iter > 1){
      p_nsize <- nsize %>%
        ggplot(aes(x = period, y = value, group = name, color = p.cat, linetype = p.cat)) +
        geom_smooth(method = 'gam', linewidth = 0.5, alpha = 0.2)
    }
    
    p_nsize <- p_nsize + 
      facet_wrap(vars(type), nrow = 1) +
      theme_bw() +
      coord_cartesian(xlim = c(0, (periods_s + 4)),  ylim = c(0, max(ceiling(nsize$value/10))*10), clip = 'off') + 
      scale_x_continuous(breaks = c(seq(1, periods_s, 1))) +
      scale_colour_manual(name = 'Position in last period',
                          breaks = levels(post$p.cat),
                          labels = levels(post$p.cat),
                          values = paste0('gray', seq(1, 90, floor((90-1)/K_s)))) +
      ylab('Cumulative sample size') + xlab('Number of period') + labs(caption = pcaption) +
      geom_text_repel(data = nsize_t, aes(label = p.cat), nudge_x = 4, hjust = 1, segment.size = .2,
                      seed = 343, direction = 'y', size = 3) +
      ptheme
    
    print('# MEAN VALUE ESTIMATE & CONFIDENCE RANGE')
    mean <- plyr::ldply(names(sim), function(x){
      m <- v <- NULL
      for (i in c('sample', 'ipw', 'haj', 'aipw', 'awaipw')){
        m <- bind_rows(m, data.frame(sim[[x]][[paste0(i, '_mean')]]) %>% mutate(e = i, type = x, iter = 1:iter_s) %>% 
                         pivot_longer(cols = starts_with('X'), names_to = 'arm', values_to = 'est'))
        v <- bind_rows(v, data.frame(sim[[x]][[paste0(i, '_var')]])  %>% mutate(e = i, type = x, iter = 1:iter_s) %>% 
                         pivot_longer(cols = starts_with('X'), names_to = 'arm', values_to = 'var'))
      }
      d <- merge(m, v)
    })
    mean$arm  <- factor(gsub('X', 'T', mean$arm), levels = paste0('T', 1:K_s))
    mean$e    <- ifelse(mean$e == 'sample', 'Naive', ifelse(mean$e == 'haj', 'Hajek', ifelse(mean$e == 'awaipw', 'AW-AIPW', toupper(mean$e))))
    mean$e    <- factor(mean$e, level = c('Naive', 'IPW', 'Hajek', 'AIPW', 'AW-AIPW'))
    mean$se   <- sqrt(mean$var)
    mean$ci_radius <- qnorm(1 - 0.05 / 2) * mean$se
    mean$true <- 0
    for (i in 1:K_s) mean$true[mean$arm == paste0('T', i)] <- true_theta_s[i]
    mean$bias <- mean$true - mean$est 
    
    p_bias <- plotit(mean, 'bias') + geom_vline(xintercept = 0, color = 'red', linetype = 'dotted') + theme(panel.spacing = unit(0.01, "lines"))
    p_ci   <- plotit(mean, 'ci_radius') + theme(panel.spacing = unit(0.01, "lines"))
    
    # TABLE FORMAT
    stat <- mean %>%
      group_by(e, type, arm) %>%
      summarize(true = min(true),
                est  = mean(est),
                se   = mean(se),
                ci_radius = mean(ci_radius),
                rmse = mean((est - true)^2)) %>%
      mutate(rmse = sqrt(rmse)) %>%
      pivot_wider(id_cols = c(type, arm, true), names_from = e, values_from = c(rmse, ci_radius)) 
    stat_sum <- stat[,1:8] %>%
      rename_all(~gsub('rmse_', '', .x)) %>%
      mutate(stat = 'RMSE') %>%
      bind_rows(stat[,c(1:3,9:13)] %>% rename_all(~gsub('ci_radius_', '', .x)) %>% mutate(stat = 'CI Radius')) %>%
      arrange(stat, arm, type) %>%
      dplyr::select(stat, type, arm, everything()) %>%
      mutate_if(is.numeric, ~ round(.x, 3))
    
    res = list(p_posterior = p_posterior, p_nsize = p_nsize, p_bias = p_bias, p_ci = p_ci, stat = stat_sum)
    return(res)
  })

  output$posterior <- renderPlot({d()$p_posterior})  
  output$nsize     <- renderPlot({d()$p_nsize})
  output$bias      <- renderPlotly({d()$p_bias})
  output$ci        <- renderPlotly({d()$p_ci})
  output$result    <- renderDataTable({d()$stat}, filter="top", options = list(searchable = TRUE, scrollX = TRUE, pageLength = 20))
}

shinyApp(ui, server, options = list(launch.browser = TRUE))
