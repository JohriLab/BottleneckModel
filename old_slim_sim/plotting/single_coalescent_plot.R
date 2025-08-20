#nolint start
rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)


source("deme_drift_bottleneck_equations.R")
source("WF_equations.R")

m_values <- seq(-4, -0.01, 0.01)
m_values <- 10^m_values

k_values <- c(10)

N <- 1e3
u <- 1e-7
d <- 100

data <- expand.grid(m = m_values, N = N, k = k_values, u = u, d = d)
simulated_data <- read.csv("summary_data.csv")

# Filter simulated data to match the k value we're plotting
simulated_data <- simulated_data %>% filter(k == 10)

# Debug: print the structure
cat("Simulated data structure:\n")
str(simulated_data)
cat("\nSimulated data columns:\n")
print(colnames(simulated_data))
cat("\nSimulated data content:\n")
print(simulated_data)

data <- data %>%
    rowwise() %>%
    mutate(
        result = list(deme_drift_bottleneck_coalescence(m, d, N, k, u)),
        WF_result = list(WF_coalescent_k(m, d, k, u))
    ) %>%
    unnest_wider(result) %>%
    unnest_wider(WF_result)

data <- data %>%
    mutate(
        piS = 2 * u * T0,
        piT = 2 * u * (
            1/N * T0 + (1 - 1/N) * T1
        ),
        FST = 1 - piS/piT,
        piS_WF = 2 * u * T0_WF,
        piT_WF = 2 * u * (
            1/N * T0_WF + (1 - 1/N) * T1_WF
        ),
        FST_WF = 1 - piS_WF/piT_WF
    )

# Create the plots
piS_plot <- ggplot(data, aes(x = m, y = piS, color = "#E69F00")) +
    geom_line(linewidth = 0.75) +
    # geom_line(aes(x = m, y = piS_WF, color = "#009E73"), linewidth = 0.75) +
    geom_point(data = simulated_data, aes(x = m, y = mean_piS), color = "#0072B2", size = 2) +
    geom_errorbar(data = simulated_data, aes(x = m, ymin = mean_piS - se_piS, ymax = mean_piS + se_piS), color = "#0072B2", width = 0.1, linewidth = 0.5, inherit.aes = FALSE) +
    facet_wrap(~k, scales = "free_y", nrow = 1) +
    scale_x_log10() +
    labs(x = "Migration Rate (m)", y = "πS") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")

piT_plot <- ggplot(data, aes(x = m, y = piT, color = "#E69F00")) +
    geom_line(linewidth = 0.75) +
    geom_hline(yintercept = 2*d*u, linetype = "dashed", color = "grey") +
    # geom_line(aes(x = m, y = piT_WF, color = "#009E73"), linewidth = 0.75) +
    geom_point(data = simulated_data, aes(x = m, y = mean_piT), color = "#0072B2", size = 2) +
    geom_errorbar(data = simulated_data, aes(x = m, ymin = mean_piT - se_piT, ymax = mean_piT + se_piT), color = "#0072B2", width = 0.1, linewidth = 0.5, inherit.aes = FALSE) +
    facet_wrap(~k, scales = "free_y", nrow = 1) +
    scale_x_log10() +
    labs(x = "Migration Rate (m)", y = "πT") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")

FST_plot <- ggplot(data, aes(x = m, y = FST, color = "#E69F00")) +
    geom_line(linewidth = 0.75) +
    # geom_line(aes(x = m, y = FST_WF, color = "#009E73"), linewidth = 0.75) +
    geom_point(data = simulated_data, aes(x = m, y = FST), color = "#0072B2", size = 2) +
    facet_wrap(~k, scales = "free_y", nrow = 1) +
    scale_x_log10() +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = "Migration Rate (m)", y = "FST") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")

# Combine plots using patchwork
combined_plot <- (FST_plot + piS_plot + piT_plot) +
    plot_layout(ncol = 1, guides = "collect")


# Save the plot
ggsave("single_coalescent_plot.png", combined_plot, width = 7, height = 12, dpi = 300)



#nolint end