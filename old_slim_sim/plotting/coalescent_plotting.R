#nolint start
rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)


source("equations/deme_drift_bottleneck_equations.R")
source("equations/WF_equations.R")
source("equations/deme_drift_no_bottleneck_equations.R")
source("equations/no_deme_drift_bottleneck_equations.R")

m_values <- seq(-4, -0.01, 0.01)
m_values <- 10^m_values

k_values <- c(1, 2, 5, 10)

N <- 1e3
u <- 1e-7
d <- 100



data <- expand.grid(m = m_values, N = N, k = k_values, u = u, d = d)

data <- data %>%
    rowwise() %>%
    mutate(
        result = list(deme_drift_bottleneck_coalescence(m, d, N, k, u)),
        WF_result = list(WF_coalescent_k(m, d, k, u)),
        deme_drift_result = list(deme_drift_coalescence(m, d, k, u)),
        no_deme_drift_bottleneck_result = list(no_deme_drift_bottleneck_coalescent(m, d, N, k, u))
    ) %>%
    unnest_wider(result) %>%
    unnest_wider(WF_result) %>%
    unnest_wider(deme_drift_result) %>%
    unnest_wider(no_deme_drift_bottleneck_result)

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
        FST_WF = 1 - piS_WF/piT_WF,
        piS_deme_drift = 2 * u * deme_drift_T0,
        piT_deme_drift = 2 * u * (
            1/N * deme_drift_T0 + (1 - 1/N) * deme_drift_T1
        ),
        FST_deme_drift = 1 - piS_deme_drift/piT_deme_drift,
        piS_no_deme_drift_bottleneck = 2 * u * no_deme_drift_bottleneck_T0,
        piT_no_deme_drift_bottleneck = 2 * u * (
            1/N * no_deme_drift_bottleneck_T0 + (1 - 1/N) * no_deme_drift_bottleneck_T1
        ),
        FST_no_deme_drift_bottleneck = 1 - piS_no_deme_drift_bottleneck/piT_no_deme_drift_bottleneck
    )

# Define consistent color mapping for all models
model_colors <- c(
  "Bottleneck Model" = "#E69F00",         # orange
  "Deme Drift" = "#0072B2",               # blue
  "No Deme Drift Bottleneck" = "#D55E00", # dark orange/red
  "Wright-Fisher" = "#009E73"             # green
)

# Create the plots
piS_plot <- ggplot(data, aes(x = m, y = piS, color = "Bottleneck Model")) +
    geom_line(linewidth = 0.75) +
    # geom_line(aes(x = m, y = piS_WF, color = "Wright-Fisher"), linewidth = 0.75) +
    # geom_line(aes(x = m, y = piS_deme_drift, color = "Deme Drift"), linewidth = 0.75) +
    # geom_line(aes(x = m, y = piS_no_deme_drift_bottleneck, color = "No Deme Drift Bottleneck"), linewidth = 0.75) +
    facet_wrap(~k, scales = "free_y", nrow = 1) +
    scale_x_log10() +
    scale_color_manual(values = model_colors) +
    labs(x = "Migration Rate (m)", y = "πS", color = "Model") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")

piT_plot <- ggplot(data, aes(x = m, y = piT, color = "Bottleneck Model")) +
    geom_line(linewidth = 0.75, alpha = 0.6) +
    geom_hline(yintercept = 2*d*u, linetype = "dashed", color = "grey") +
    # geom_line(aes(x = m, y = piT_WF, color = "Wright-Fisher"), linewidth = 0.75, alpha = 0.6) +
    # geom_line(aes(x = m, y = piT_deme_drift, color = "Deme Drift"), linewidth = 0.75, alpha = 0.6) +
    # geom_line(aes(x = m, y = piT_no_deme_drift_bottleneck, color = "No Deme Drift Bottleneck"), linewidth = 0.75, alpha = 0.6) +
    facet_wrap(~k, scales = "free_y", nrow = 1) +
    scale_x_log10() +
    scale_color_manual(values = model_colors) +
    labs(x = "Migration Rate (m)", y = "πT", color = "Model") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")

FST_plot <- ggplot(data, aes(x = m, y = FST, color = "Bottleneck Model")) +
    geom_line(linewidth = 0.75) +
    # geom_line(aes(x = m, y = FST_WF, color = "Wright-Fisher"), linewidth = 0.75) +
    # geom_line(aes(x = m, y = FST_deme_drift, color = "Deme Drift"), linewidth = 0.75) +
    # geom_line(aes(x = m, y = FST_no_deme
    facet_wrap(~k, scales = "free_y", nrow = 1) +
    scale_x_log10() +
    scale_y_continuous(limits = c(0, 1)) +
    scale_color_manual(values = model_colors) +
    labs(x = "Migration Rate (m)", y = "FST", color = "Model") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")

# Dummy plot for legend only
legend_data <- data.frame(x = 1:4, y = 1:4, model = names(model_colors))
legend_plot <- ggplot(legend_data, aes(x = x, y = y, color = model)) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = model_colors) +
    labs(color = "Model") +
    theme_minimal() +
    theme(
        legend.position = "bottom",
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        plot.background = element_blank()
    )

# Combine plots using patchwork
combined_plot <- (FST_plot + piS_plot + piT_plot + legend_plot) +
    patchwork::plot_layout(ncol = 1, heights = c(1, 1, 1, 0.15))


# Save the plot
ggsave("coalescent_plots.png", combined_plot, width = 19, height = 12, dpi = 300)



#nolint end