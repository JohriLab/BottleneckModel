#nolint start
rm(list = ls())

suppressMessages({
  library(tidyverse)
  library(ggplot2)
})

  df <- read_csv("bottleneck_simulation_results.csv", col_names = c("piS", "piT", "allele_count", "k", "m"), show_col_types = FALSE) %>%
    filter(k == 10, m == 0.4)



  mean_piS <- mean(df$piS)
  mean_piT <- mean(df$piT)

  sd_piS <- sd(df$piS)
  sd_piT <- sd(df$piT)

  se_piS <- sd_piS / sqrt(length(df$piS))
  se_piT <- sd_piT / sqrt(length(df$piT))

  FST <- 1 - (mean_piS / mean_piT)

# Create a data frame for the summary statistics
summary_stats <- data.frame(
  Statistic = c("piS", "piT", "FST"),
  Mean = c(mean_piS, mean_piT, FST),
  Standard_Error = c(se_piS, se_piT, NA)  # FST doesn't have a standard error in this context
)

  # Print the table
  print("Summary Statistics Table:")
  print(summary_stats)

  wattersons_theta <- 2*100*1e-7
  print(paste("Watterson's Theta (N = d):", wattersons_theta))


  # SFS
  n_obs <- length(df$allele_count)
  sfs_df <- data.frame(allele_count = 1:1e5)

  sfs_table <- table(df$allele_count)

  sfs_df <- sfs_df %>%
      mutate(
          bottleneck_frequency = as.numeric(sfs_table[as.character(allele_count)]) / n_obs,
          WF_frequency = mean_piT / allele_count
      )
  sfs_df$bottleneck_frequency[is.na(sfs_df$bottleneck_frequency)] <- 0
  sfs_df <- sfs_df %>%
      filter(allele_count > 1) %>%
      mutate(
          bottleneck_frequency = bottleneck_frequency / sum(bottleneck_frequency),
          WF_frequency = WF_frequency / sum(WF_frequency)
      )


  sfs_df_plot <- sfs_df %>%
      pivot_longer(
          cols = c(bottleneck_frequency, WF_frequency),
          names_to = "type",
          values_to = "frequency"
      )

  plot <- sfs_df_plot %>%
      filter(allele_count <=20, allele_count >= 2) %>%
      ggplot(aes(x = allele_count, y = frequency, fill = type)) +
      geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
      scale_fill_manual(values = c("bottleneck_frequency" = "#E69F00", "WF_frequency" = "#009E73"),
                      labels = c("Bottleneck", "Wright-Fisher")) +
    labs(x = "Allele Count", y = "Frequency", title = "Site Frequency Spectrum Comparison") +
      theme_minimal() +
      theme(legend.title =  element_blank(),
            panel.background = element_rect(fill = "white"),
            plot.background = element_rect(fill = "white"))

ggsave("sfs_comparison_plot.png", plot, width = 10, height = 6, dpi = 300)



# Group by m and k, then calculate summary statistics
summary_data <- df %>%
  group_by(m, k) %>%
  summarise(
    mean_piS = mean(piS),
    se_piS = sd(piS) / sqrt(n()),
    mean_piT = mean(piT),
    se_piT = sd(piT) / sqrt(n()),
    FST = FST,
    .groups = 'drop'
  )
# Write summary statistics table to CSV
write.csv(summary_data, "plotting/summary_data.csv", row.names = FALSE)



    
#nolint end