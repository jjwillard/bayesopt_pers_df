### Example plotting functions used in "Bayesian optimization for personalized
### dose-finding trials with combination therapies"
### 2023-10-26

pacman::p_load(dplyr, tidyr, purrr, ggplot2, cowplot, foreach, doParallel, latex2exp)

##################### Data #############################

regular <- readRDS("PATH/bimodal_sim_data_standard.RDS")
personalized <- readRDS("PATHbimodal_sim_data_personalized.RDS")

### Example objective function

objective_function <- function(x_1, x_2, x_3){
  data_x <- as.matrix(tibble(x_1 = x_1, x_2 = x_2))
  
  # same for scenario 3
  if(x_3 == 0){
    f_x <-  -mvtnorm::dmvnorm(data_x,
                              mean = c(0.3, 0.7),
                              sigma = matrix(c(0.2, 0.05, 0.05, 0.1), ncol = 2))
  } else {
    f_x <-  -mvtnorm::dmvnorm(data_x,
                              mean = c(0.7, 0.3),
                              sigma = matrix(c(0.2, 0.05, 0.05, 0.1), ncol = 2))
  }
  tibble(f_x = f_x)
  
}

optimal_value_x_3_0 <- pull(objective_function(0.3, 0.7, 0))
optimal_value_x_3_1 <- pull(objective_function(0.7, 0.3, 1))
optimal_dose_combo_x_3_0 <- c(0.3, 0.7)
optimal_dose_combo_x_3_1 <- c(0.7, 0.3)


full_data <- regular %>% mutate(algorithm = "standard") %>%
  bind_rows(personalized %>% mutate(algorithm = "personalized"))

## labels by iteration/sample size
tick_labels <- c("1\n24", "2\n28", "3\n32", "4\n36", "5\n40", "6\n44", "7\n48", 
                 "8\n52", "9\n56", "10\n60", "11\n64", "12\n68", "13\n72", "14\n76",
                 "15\n80")

########################### Summarize data ###############################
summarized_data <- full_data  %>%
  filter(iteration >= 1) %>%
  select(simulation, algorithm, iteration, mu_rec, euc_dist_dose, rmse, x_3) %>%
  distinct() %>%
  group_by(x_3, iteration, algorithm) %>%
  summarise(mean_rmse = mean(rmse),
         median_rmse = median(rmse),
         mean_mu_rec = mean(mu_rec),
         median_mu_rec = median(mu_rec),
         mean_euc_dist_dose = mean(euc_dist_dose),
         median_euc_dist_dose = median(euc_dist_dose),.groups = "drop") %>%
  mutate(optimal_value = if_else(x_3 == 0, optimal_value_x_3_0, optimal_value_x_3_1),
         x_3 = factor(if_else(x_3 == 0,
                              "z = 0",
                              "z = 1")))
levels(summarized_data$x_3) <- c(`z = 0` = TeX("$z_1 = 0$"), `z = 1` = TeX("$z_1 = 1$"))

 
############################### Plotting function ###############################

plot_summary <- function(data, summary_metric, ylab, facet_scales = "free_y",
                         parse_ylab = TRUE){
  summary_metric <- enquo(summary_metric)
  p <- ggplot(data, 
              aes(x = iteration, y = !!summary_metric, color = algorithm, shape = algorithm)) +
    facet_grid(x_3 ~ ., scales = facet_scales, labeller = label_parsed) +
    geom_point() +
    geom_line(linetype = "dotted") +
    scale_x_continuous(breaks = unique(data$iteration),
                       labels = tick_labels) +
    scale_color_manual(values = c("#26495c", "#c66b3d")) +
    theme_bw() +
    xlab("iteration \n total number of participants") +
    theme(panel.grid.minor = element_blank(),
          legend.position = "top",
          legend.title = element_blank(),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 7))
  
  if(parse_ylab == TRUE){
    p + ylab(parse(text = ylab)) 
  } else {
    p + ylab(ylab)
  }
}

## Euclidean distance plot
euclidean_distance_ylab <- TeX("expected $||\\widehat{\\textbf{d}}_{opt} - \\textbf{d}_{opt}||$")
euc_dist_dose_plot <- plot_summary(summarized_data,
                                   summary_metric = mean_euc_dist_dose,
                                   ylab = euclidean_distance_ylab)


## RPSEL plot
rmse_ylab <- "expected rpsel"
rmse_plot <- plot_summary(summarized_data,
                          summary_metric = mean_rmse,
                          ylab = rmse_ylab,
                          parse_ylab = FALSE)


## Posterior mean plot
post_mean_ylab <- TeX("expected $E\\lbrack f | D, \\widehat{\\textbf{d}}\\rbrack$")
post_mean_plot <- plot_summary(summarized_data,
                               summary_metric = mean_mu_rec,
                               ylab = post_mean_ylab) +
  geom_hline(aes(yintercept = optimal_value), linetype = "dashed") +
  coord_cartesian(ylim = c(-1.3, NA))


#################################### Objective function plot ####################################

ex_data <- expand_grid(x_1 = seq(0, 1, 0.01),
                       x_2 = seq(0, 1, 0.01),
                       x_3 = c(0,1))


doParallel::registerDoParallel(cores = parallel::detectCores())
f_x <- foreach(i = 1:nrow(ex_data), .combine = 'rbind') %dopar% {
  objective_function(ex_data$x_1[i], 
                     ex_data$x_2[i],
                     ex_data$x_3[i])
} 
doParallel::stopImplicitCluster()

plot_data <- bind_cols(ex_data, f_x) %>%
  mutate(x_3 = factor(if_else(x_3 == 0, "z = 0", "z = 1")))
levels(plot_data$x_3) <- c(`z = 0` = TeX("$z_1 = 0$"), `z = 1` = TeX("$z_1 = 1$"))

point_data <- tibble(x_3 = factor(c("z = 0", "z = 1")),
                     x_1 = c(0.3, 0.7),
                     x_2 = c(0.7, 0.3))
levels(point_data$x_3) <- c(`z = 0` = TeX("$z_1 = 0$"), `z = 1` = TeX("$z_1 = 1$"))

true_dgm_plot <- ggplot(plot_data, aes(x = x_1, y = x_2, z = f_x)) +
  facet_wrap(x_3~., labeller = label_parsed) +
  geom_contour_filled(breaks = c(-2, -1.33, -1.07, -0.8, -0.53, -0.27, 2)) + 
  geom_point(data = point_data,
             aes(x = x_1, y = x_2), shape = 8, color = "white", inherit.aes = FALSE) +
  coord_cartesian(xlim = c(-0.1, 1.1), ylim = c(-0.1, 1.1)) +
  guides(fill = guide_legend(title = "", nrow = 2, byrow = FALSE)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) + #match scale of other plots
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_fill_viridis_d(option = "D") +
  ylab(parse(text = TeX('$d_2$')))+
  xlab(parse(text = TeX('$d_1$')))+
  theme(panel.background = element_rect(fill='darkgrey'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 8),
        legend.position = "none",
        aspect.ratio = 1)



############################## Plots for single run of algorithms ############################## 

plot_mu <- function(fits, evaluated_points = NULL, nrows_legend = 2,
                    legend_position = "bottom",
                    labeller_parsed = TRUE){
  
  if(labeller_parsed == TRUE){
    p <- ggplot(fits, aes(x = x_1, y = x_2, z = mu)) +
      facet_wrap(x_3~., labeller = label_parsed) +
      geom_contour_filled(breaks = c(-2, -1.33, -1.07, -0.8, -0.53, -0.27, 2))
  } else {
    p <- ggplot(fits, aes(x = x_1, y = x_2, z = mu)) +
      facet_wrap(x_3~.) +
      geom_contour_filled(breaks = c(-2, -1.33, -1.07, -0.8, -0.53, -0.27, 2))
  }
  
  
  
  if(!is.null(evaluated_points)){
    p <- p + 
      geom_point(data = evaluated_points %>%
                   filter(iteration < 1),
                 aes(x = x_1, y = x_2), 
                 inherit.aes = FALSE) + 
      geom_point(data = evaluated_points %>%
                   filter(iteration >= 1),
                 aes(x = x_1, y = x_2), color = "white",
                 inherit.aes = FALSE) + 
      ggrepel::geom_text_repel(data = evaluated_points %>% filter(iteration >= 1), aes(x = x_1, y = x_2, label = iteration),
                               inherit.aes = FALSE, 
                               color = "white", size = 2.5, max.overlaps = 100) +
      coord_cartesian(xlim = c(-0.1, 1.1), ylim = c(-0.1, 1.1)) +
      guides(fill = guide_legend("", nrow = nrows_legend, byrow = FALSE)) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill='darkgrey'))
  }
  
  p +
    scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) + #match scale of other plots
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) + #match scale of other plots
    scale_fill_viridis_d(option = "D", drop = FALSE) +
    ylab(parse(text = TeX('$d_2$')))+
    xlab(parse(text = TeX('$d_1$')))+
    theme(legend.text = element_text(size = 8),
          legend.position = legend_position,
          aspect.ratio = 1) 
}

single_run_regular <- readRDS("PATH/single_run_bimodal_resp_het_regular.RDS") 
single_run_regular_fits <- single_run_regular$yhat %>%
  mutate(x_3 = "")
single_run_regular_evals <- single_run_regular$evals %>% 
  select(x_1, x_2, iteration) %>% 
  distinct() %>%
  mutate(x_3 = "")

single_run_personalized <- readRDS("PATH/single_run_bimodal_resp_het_personalized.RDS")
single_run_personalized_fits <- single_run_personalized$yhat %>%
  mutate(x_3 = factor(if_else(x_3 == 0, "z = 0", "z = 1")))
levels(single_run_personalized_fits$x_3) <- c(`z = 0` = TeX("$z_1 = 0$"), `z = 1` = TeX("$z_1 = 1$"))
single_run_personalized_evals <- single_run_personalized$evals %>% 
  select(x_1, x_2, x_3, iteration) %>% 
  distinct() %>%
  mutate(x_3 = factor(if_else(x_3 == 0, "z = 0", "z = 1")))
levels(single_run_personalized_evals$x_3) <- c(`z = 0` = TeX("$z_1 = 0$"), `z = 1` = TeX("$z_1 = 1$"))

mu_plot_regular <- plot_mu(single_run_regular_fits, 
                           single_run_regular_evals,
                           legend_position = "none",
                           labeller_parsed = FALSE) 

mu_plot_personalized <- plot_mu(single_run_personalized_fits,
                                single_run_personalized_evals,
                                legend_position = "none")

############### Final plot used in manuscript under Scenario 2 (bimodal response heterogeneity) #######

plot_legend <- get_legend(plot_mu(single_run_personalized$yhat, 
                                  single_run_personalized$evals %>% select(x_1, x_2, x_3, iteration) %>% distinct(),
                                  nrows_legend = 1,
                                  legend_position = "bottom"))

prow <- plot_grid(true_dgm_plot,
                  mu_plot_regular,
                  mu_plot_personalized,
                  labels = LETTERS[1:3],
                  align = "h",
                  axis = "bt",
                  hjust = -1,
                  ncol = 3,
                  rel_widths = c(1.75, 1, 1.75))

prow_legend <- plot_grid(prow, plot_legend, ncol = 1, rel_heights = c(1, .1))


rate_prop <- plot_grid(euc_dist_dose_plot,
                       rmse_plot,
                       post_mean_plot, 
                       nrow = 1,
                       labels = c("D", "E", "F"))

final_plot <- plot_grid(prow_legend,
          rate_prop,
          nrow = 2,
          scale = 0.95)

ggsave2(filename = "PATH/bimodal_resp_het.eps",
        plot = final_plot,
        width = 11,
        height = 6.5,
        units = "in",
        dpi = 800)
