#nolint start
rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(scales)


source("analysis/equations/coalescence/deme_drift_bottleneck_equations.R")

m_values <- seq(-4.5, -0.01, 0.01)
m_values <- 10^m_values

N <- 1e3
u <- 1e-6
d <- 100


simulated_data <- read.csv("analysis/pi_data_plotting.csv")
k_values <- unique(simulated_data$k)

data <- expand.grid(m = m_values, N = N, k = k_values, u = u, d = d)


data <- data %>%
    rowwise() %>%
    mutate(
        result = list(deme_drift_bottleneck_coalescence(m, d, N, k, u))
    ) %>%
    unnest_wider(result)

data <- data %>%
    mutate(
        piT = 2 * u * (
            1/N * T0 + (1 - 1/N) * T1
        )
    )



for (k in k_values) {
# Filter simulated data to match the k value we're plotting
simulated_data <- simulated_data %>% filter(k == k)

    # piS_plot <- ggplot(data, aes(x = m, y = piS, color = "#E69F00")) +
    #     geom_line(linewidth = 0.75) +
    #     # geom_line(aes(x = m, y = piS_WF, color = "#009E73"), linewidth = 0.75) +
    #     geom_point(data = simulated_data, aes(x = m, y = piS), color = "#0072B2", size = 2) +
    #     geom_errorbar(data = simulated_data, aes(x = m, ymin = piS - piS_se, ymax = piS + piS_se), color = "#0072B2", width = 0.1, linewidth = 0.5, inherit.aes = FALSE) +
    #     facet_wrap(~k, scales = "free_y", nrow = 1) +
    #     scale_x_log10() +
    #     scale_y_continuous(labels = scientific) +
    #     labs(x = "Migration Rate (m)", y = "πS") +
    #     theme_minimal(base_size = 12) +
    #     theme(legend.position = "none")

    piT_plot <- ggplot(data, aes(x = m, y = piT, color = "#E69F00")) +
        geom_line(linewidth = 0.75) +
        geom_hline(yintercept = 2*d*u, linetype = "dashed", color = "grey") +
        # geom_line(aes(x = m, y = piT_WF, color = "#009E73"), linewidth = 0.75) +
        geom_point(data = simulated_data, aes(x = m, y = pi), color = "#0072B2", size = 2) +
        geom_errorbar(data = simulated_data, aes(x = m, ymin = pi - se_pi, ymax = pi + se_pi), color = "#0072B2", width = 0.1, linewidth = 0.5, inherit.aes = FALSE) +
        facet_wrap(~k, scales = "free_y", nrow = 1) +
        scale_x_log10() +
        scale_y_continuous(labels = scientific) +
        labs(x = "Migration Rate (m)", y = "πT") +
        theme_minimal(base_size = 12) +
        theme(legend.position = "none")

    # FST_plot <- ggplot(data, aes(x = m, y = FST, color = "#E69F00")) +
    #     geom_line(linewidth = 0.75) +
    #     # geom_line(aes(x = m, y = FST_WF, color = "#009E73"), linewidth = 0.75) +
    #     geom_point(data = simulated_data, aes(x = m, y = FST), color = "#0072B2", size = 2) +
    #     # geom_errorbar(data = simulated_data, aes(x = m, ymin = FST - FST_se, ymax = FST + FST_se), color = "#0072B2", width = 0.1, linewidth = 0.5, inherit.aes = FALSE) +
    #     facet_wrap(~k, scales = "free_y", nrow = 1) +
    #     scale_x_log10() +
    #     scale_y_continuous(limits = c(0, 1)) +
    #     labs(x = "Migration Rate (m)", y = "FST") +
    #     theme_minimal(base_size = 12) +
    #     theme(legend.position = "none")

    combined_plot <- (piT_plot) +
        plot_layout(ncol = 1, guides = "collect")



}
scale <- 1.5
ggsave("analysis/coalescent_plot.png", combined_plot, width = 10*scale, height = 2.5*scale, dpi = 300)
#nolint end