#nolint start
library(tidyverse)
library(ggplot2)
library(boot)

sample_size <- 40

sfs_data <- read.csv("analysis/sfs.csv") %>%
    select(-replicate)

k_values <- unique(sfs_data$k)
m_values <- unique(sfs_data$m)

# Corrected replicate assignment logic
sfs_data$replicate <- 0
rep_id <- 0
cumulative_count <- 0

# Assign correct replicate numbers to each row

for (i in 1:nrow(sfs_data)) {
    cumulative_count <- cumulative_count + sfs_data$count[i]
    sfs_data$replicate[i] <- rep_id
    
    # When we reach exactly 100,000, increment the replicate ID
    if (cumulative_count == 100000) {
        rep_id <- rep_id + 1
        cumulative_count <- 0
    }
    
    if (i %% 100000 == 0) {
        print(paste0("Processed ", i, " rows, current replicate: ", rep_id))
    }
}

# Verify that each replicate sums to exactly 100,000
incorrect_replicates <- sfs_data %>%
    group_by(k, m, replicate) %>%
    summarise(total_count = sum(count), .groups = 'drop') %>%
    filter(total_count != 100000)

if (nrow(incorrect_replicates) > 0) {
    print("Warning: Found replicates that don't sum to 100,000:")
    print(incorrect_replicates)
} else {
    print("All replicates sum to exactly 100,000 ✓")
}

print(paste("Total number of replicates:", max(sfs_data$replicate) + 1))

# Pivot the data so that each row is a replicate, and each count is a column
sfs_data <- sfs_data %>%
    group_by(k, m, replicate) %>%
    mutate(count_id = row_number()) %>%
    pivot_wider(
        id_cols = c(k, m, replicate),
        names_from = count_id,
        values_from = count,
        names_prefix = "count_",
        values_fill = list(count = 0)
    )


write.csv(sfs_data, "analysis/sfs_replicate_corrected.csv", row.names = FALSE)

sfs_data <- read.csv("analysis/sfs_replicate_corrected.csv")


# For each row (replicate), sample sample_size individuals and store the resulting mutation frequencies as a vector

# Identify the count columns
count_cols <- grep("^count_", colnames(sfs_data), value = TRUE)

# Function to sample SFS for a single replicate row
sampled_alleles_vector <- function(row_data, sample_size = sample_size) {
    # row_data is a data frame row, so we can access columns by name
    counts <- as.numeric(row_data[count_cols])
    mutation_classes <- rep(seq_along(counts), counts)
    # replace_flag <- length(mutation_classes) < sample_size
    sampled <- sample(mutation_classes, size = sample_size, replace = FALSE)
    sampled_sfs <- tabulate(sampled, nbins = length(counts))
    return(sampled_sfs)
}

sampled_alleles_vector_simplified <- function(row_data, sample_size = sample_size) {
    
}



sfs_data$sampled_alleles <- lapply(1:nrow(sfs_data), function(row_idx) {
    sampled_alleles_vector(sfs_data[row_idx, count_cols], sample_size = sample_size)
})


# Save the full data with list column as RDS (preserves all data types)
saveRDS(sfs_data, "analysis/sfs_sampled.rds")


# Calculate the SFS from the sampled alleles for each unique m and k combination
sfs_sampled_data <- readRDS("analysis/sfs_sampled.rds")

# Function to calculate folded SFS for a given dataset
calculate_folded_sfs <- function(allele_data, sample_size) {
    # Convert list of vectors to matrix
    allele_matrix <- do.call(rbind, allele_data)
    
    # Calculate folded SFS
    folded_sfs_df <- data.frame(i = 1:ceiling((sample_size - 1)/2), frequency = 0)
    
    for (i in 1:ceiling((sample_size - 1)/2)) {
        folded_sfs_df[folded_sfs_df$i == i, "frequency"] <- sum(allele_matrix == i) + sum(allele_matrix == (sample_size - i))
        if (i == sample_size - i) {
            folded_sfs_df[folded_sfs_df$i == i, "frequency"] <- folded_sfs_df[folded_sfs_df$i == i, "frequency"] / 2
        }
    }

    # Make frequency SFS instead of count SFS
    # Calculate proportion of loci (not sites)
    folded_sfs_df$frequency <- folded_sfs_df$frequency / length(allele_data)

    # Calculate π from SFS
    theta <- sum(folded_sfs_df$i * (sample_size - folded_sfs_df$i) * folded_sfs_df$frequency) / (sample_size * (sample_size - 1) / 2)
    
    # Calculate expected SFS (theta/i) scaled to match total sum
    expected_freq <- theta / folded_sfs_df$i
    
    # Add expected frequencies to dataframe
    folded_sfs_df$expected_frequency <- expected_freq
    
    return(folded_sfs_df)
}

# Function to calculate folded SFS for bootstrap resampling
bootstrap_folded_sfs <- function(allele_data, indices, sample_size) {
    # Resample the allele data
    resampled_data <- allele_data[indices]
    
    # Calculate folded SFS for resampled data
    allele_matrix <- do.call(rbind, resampled_data)
    
    # Calculate folded SFS
    folded_sfs <- rep(0, ceiling((sample_size - 1)/2))
    
    for (i in 1:ceiling((sample_size - 1)/2)) {
        folded_sfs[i] <- sum(allele_matrix == i) + sum(allele_matrix == (sample_size - i))
        if (i == sample_size - i) {
            folded_sfs[i] <- folded_sfs[i] / 2
        }
    }
    
    # Convert to proportions
    folded_sfs <- folded_sfs / length(resampled_data)
    
    return(folded_sfs)
}

# Function to calculate folded SFS with error bars for a given dataset
calculate_folded_sfs_with_errors <- function(allele_data, sample_size, n_bootstrap = 1000) {
    # Calculate the main SFS
    main_sfs <- calculate_folded_sfs(allele_data, sample_size)
    
    # Perform bootstrap resampling
    bootstrap_results <- boot(allele_data, bootstrap_folded_sfs, R = n_bootstrap, sample_size = sample_size)
    
    # Calculate confidence intervals for each frequency bin
    ci_lower <- rep(0, nrow(main_sfs))
    ci_upper <- rep(0, nrow(main_sfs))
    
    for (i in 1:nrow(main_sfs)) {
        # Get bootstrap samples for this frequency bin
        bootstrap_samples <- bootstrap_results$t[, i]
        
        # Calculate 95% confidence intervals
        ci_lower[i] <- quantile(bootstrap_samples, 0.025)
        ci_upper[i] <- quantile(bootstrap_samples, 0.975)
    }
    
    # Add confidence intervals to the dataframe
    main_sfs$ci_lower <- ci_lower
    main_sfs$ci_upper <- ci_upper
    
    return(main_sfs)
}

# Calculate folded SFS with error bars for each unique m and k combination
all_folded_sfs <- list()

for (k_val in unique(sfs_sampled_data$k)) {
    for (m_val in unique(sfs_sampled_data$m)) {
        # Filter data for this specific k and m combination
        subset_data <- sfs_sampled_data %>%
            filter(k == k_val, m == m_val)
        
        if (nrow(subset_data) > 0) {
            # Calculate folded SFS with error bars for this subset
            folded_sfs <- calculate_folded_sfs_with_errors(subset_data$sampled_alleles, sample_size)
            
            # Add metadata columns
            folded_sfs$k <- k_val
            folded_sfs$m <- m_val
            
            # Store in list with descriptive name
            list_name <- paste0("k_", k_val, "_m_", m_val)
            all_folded_sfs[[list_name]] <- folded_sfs
        }
        print(paste("Processed:", k_val, m_val))
    }
}

saveRDS(all_folded_sfs, file = "analysis/all_folded_sfs.rds")


# Combine all folded SFS data into one dataframe
combined_folded_sfs <- do.call(rbind, all_folded_sfs)

# Reset row names
rownames(combined_folded_sfs) <- NULL

# Create faceted plot with error bars
faceted_plot_with_errors <- combined_folded_sfs %>%
    ggplot() +
    # Observed SFS bars with error bars
    geom_bar(aes(x = i, y = frequency), stat = "identity", fill = "#2E86AB", alpha = 0.9, color = "#2E86AB") +
    geom_errorbar(aes(x = i, ymin = ci_lower, ymax = ci_upper), 
                  width = 0.2, color = "#34495E", linewidth = 0.5) +
    # Expected SFS (theta/i) as red points connected by lines
    geom_point(aes(x = i, y = expected_frequency), color = "#E74C3C", linewidth = 1.2) +
    # geom_point(aes(x = i, y = expected_frequency), color = "#E74C3C", size = 2.5, fill = "white", shape = 21) +
    facet_grid(m ~ k, labeller = label_both, scales = "free_y") +
    labs(
        x = "Derived Allele Count", 
        y = "Frequency"
    ) +
    theme_light() +
    theme(
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12, face = "bold", color = "#2C3E50"),
        strip.background = element_rect(fill = "#ECF0F1", color = "#BDC3C7"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold", color = "#2C3E50", margin = margin(b = 10)),
        plot.subtitle = element_text(hjust = 0.5, size = 12, color = "#7F8C8D", margin = margin(b = 20)),
        plot.caption = element_text(hjust = 0.5, size = 12, color = "#7F8C8D", margin = margin(t = 10)),
        axis.title = element_text(size = 12, face = "bold", color = "#2C3E50"),
        axis.text = element_text(size = 12, color = "#34495E"),
        panel.grid.major = element_line(color = "#E8E8E8", linewidth = 0.5),
        panel.grid.minor = element_line(color = "#F5F5F5", linewidth = 0.3),
        panel.border = element_rect(color = "#BDC3C7", linewidth = 0.5)
    )

# Save the plot
ggsave("analysis/sfs_faceted_by_k_and_m_with_errors_40.png", faceted_plot_with_errors, width = 14, height = 12, dpi = 300)

# # Also save the combined data with error bars
# write.csv(combined_folded_sfs, "analysis/folded_sfs_by_k_and_m_with_errors.csv", row.names = FALSE)

# Print summary statistics
print("Summary of folded SFS with error bars:")
print(paste("Total combinations:", length(all_folded_sfs)))
print(paste("Unique k values:", length(unique(combined_folded_sfs$k))))
print(paste("Unique m values:", length(unique(combined_folded_sfs$m))))
print(paste("Bootstrap iterations per combination:", 1000))


print("Number of variant sites at each point:")


length(which(sapply(sfs_sampled_data$sampled_alleles, function(x) x[1] != 40)))

for (k_val in unique(sfs_sampled_data$k)) {
    for (m_val in unique(sfs_sampled_data$m)) {
        subset_data <- sfs_sampled_data %>%
            filter(k == k_val, m == m_val)
        cat("k =", k_val, "m =", m_val, "count =", length(which(sapply(subset_data$sampled_alleles, function(x) x[1] != 40))), "\n")
    }
}


#nolint end
