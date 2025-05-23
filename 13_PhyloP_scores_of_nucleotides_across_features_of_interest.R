## 13th (bonus) ##

```{r}
library("tidyverse")
library("patchwork")
library("ggpubr")
library("rtracklayer")
library("GenomicRanges")
```

# import .bw (conservation score) file:
```{r}
bw <- import.bw("ce11.phyloP26way.bw")
```

# Extract start and end coordinates of each IR event:
```{r}
extracted_coords <- lapply(intron_sums_heatmap$AS_event_ID, function(id) {
  parts <- unlist(strsplit(id, "_"))
  chr <- as.character(parts[1])
  start <- as.numeric(parts[3])
  end <- as.numeric(parts[4])
  return(c(chr, start, end))
})

extracted_coords <- do.call(rbind, extracted_coords)
rownames(extracted_coords) <- intron_sums_heatmap$AS_event_ID
#write.csv(extracted_coords, file = "extracted_coords_IR_events.csv")
```

# calculate mean phylop scores between extracted coordinates ranges
```{r}
# Initialize output data frame
phyloP_scores_combined <- data.frame(AS_event_ID = character(), Chr = character(), Start = integer(), 
                        End = integer(), Average_Score = numeric(), stringsAsFactors = FALSE)

# Process each row of extracted_coords
for (i in seq_len(nrow(extracted_coords))) {
  AS_event_ID <- rownames(extracted_coords)[i]
  chr <- paste0("chr", extracted_coords[i, 1])
  start <- as.integer(extracted_coords[i, 2])
  end <- as.integer(extracted_coords[i, 3])
  
  # Check if start and end are valid
  if (!is.na(start) && !is.na(end) && start < end) {
    # Extract phyloP scores from the BigWig file
    scores <- import.bw("ce11.phyloP26way.bw", which = GRanges(chr, IRanges(start, end)))$score
    
    # Compute average score
    average <- if (length(scores) > 0) mean(scores, na.rm = TRUE) else NA
    
  # Debugging output
  cat(sprintf("Processing: AS_event_ID=%s, Chr=%s, Start=%d, End=%d, average=%f\n", AS_event_ID, chr, start, end, average))
  
    # Append result to output data frame
    phyloP_scores_combined <- rbind(phyloP_scores_combined, data.frame(AS_event_ID, chr, start, end, average, stringsAsFactors = FALSE))
  } else {
    cat(sprintf("Skipping invalid coordinates: AS_event_ID=%s (%s, %d, %d)\n", AS_event_ID, chr, start, end))
  }
}
```

# Modify and format phyloP_scores_combined into intron_phylop:
```{r}
intron_phylop <- as.data.frame(phyloP_scores_combined)
col_name <- "PhyloP score"
colnames(intron_phylop)[5] <- col_name
```

# Correlation of uniqueness index values with PhyloP score:
```{r}
intron_sums_stats <- intron_sums_heatmap %>%
  rowwise() %>%
  summarize(AS_event_ID = AS_event_ID, Min_Value = min(c_across(-c(AS_event_ID, shift_or_preserve)), na.rm = TRUE), Max_Value = max(c_across(-c(AS_event_ID, shift_or_preserve)), na.rm = TRUE), shift_or_preserve = first(shift_or_preserve))

# "bin" min and max values into 5 separate percentiles:
intron_sums_stats <- intron_sums_stats %>%
  mutate(across(c(Min_Value, Max_Value), ~ cut(.x, quantile(.x, probs = seq(0, 1, 0.2), na.rm = TRUE), 
                                               labels = c("0-20", "21-40", "41-60", "61-80", "81-100"), include.lowest = TRUE), 
                .names = "{.col}_Percentile"))

combined_data <- merge(intron_phylop, intron_sums_stats, by = "AS_event_ID")
combined_data <- combined_data[,!names(combined_data) %in% c("shift_or_preserve")]
combined_data <- combined_data[!is.na(combined_data$`PhyloP score`),]

lm_model <- lm((`Max_Value`) ~ `PhyloP score`, data = combined_data)
model_summary <- summary(lm_model)

slope <- coef(lm_model)[2]
intercept <- coef(lm_model)[1]
r_squared <- model_summary$adj.r.squared
p_value <- coef(model_summary)[2,4]

custom_label <- paste0(
  "y = ", round(intercept, 3), " + ", round(slope, 3), "x",
  "\nAdj. R² = ", round(r_squared, 3),
  "\nP = ", format.pval(p_value, digits = 3))

ggplot(combined_data, aes(x = Max_Value, y = `PhyloP score`)) +
  geom_point(color = "dodgerblue", size = 2) +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  annotate("text", x = max(combined_data$Max_Value, na.rm = TRUE) * 0.7, y = max(combined_data$`PhyloP score`, na.rm = TRUE) * 0.9, label = custom_label, size = 5, color = "black", hjust = 0) +
  labs(x = "Max Uniqueness Value", y = "Mean PhyloP Score", title = "Linear Regression: Max Uniqueness Value vs Mean PhyloP Score for Each IR Event") +
  theme_minimal()

lm_model <- lm((`Min_Value`) ~ `PhyloP score`, data = combined_data)
model_summary <- summary(lm_model)

slope <- coef(lm_model)[2]
intercept <- coef(lm_model)[1]
r_squared <- model_summary$adj.r.squared
p_value <- coef(model_summary)[2,4]

custom_label <- paste0(
  "y = ", round(intercept, 3), " + ", round(slope, 3), "x",
  "\nAdj. R² = ", round(r_squared, 3),
  "\nP = ", format.pval(p_value, digits = 3))

ggplot(combined_data, aes(x = Min_Value, y = `PhyloP score`)) +
  geom_point(color = "dodgerblue", size = 2) +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  annotate("text", x = min(combined_data$Min_Value, na.rm = TRUE) * 0.7, y = max(combined_data$`PhyloP score`, na.rm = TRUE) * 0.9, label = custom_label, size = 5, color = "black", hjust = 0) +
  labs(x = "Min Uniqueness Value", y = "Mean PhyloP Score", title = "Linear Regression: Min Uniqueness Value vs Mean PhyloP Score for Each IR Event")+
  theme_minimal()
```

# Do the same linear regression, but alter it for unique neuron types:
```{r}
regression_results <- data.frame(Column = character(),Intercept = numeric(),Slope = numeric(),R_Squared = numeric(),P_Value = numeric(),stringsAsFactors = FALSE)

for (i in singletypes) {
  uniqueness_values <- intron_sums_heatmap
  valid_indices <- row_number(intersect(phyloP_scores_combined$AS_event_ID, uniqueness_values$AS_event_ID))
  uniqueness_values <- intron_sums_heatmap[valid_indices,i]
  regression_data <- data.frame(uniqueness_values, PhyloP_score = intron_phylop$`PhyloP score`)
  regression_data <- regression_data[complete.cases(regression_data),]
  regression_data <- regression_data[regression_data$uniqueness_values != 0,]
  
  lm_model <- lm(PhyloP_score ~ uniqueness_values, data = regression_data)
  model_summary <- summary(lm_model)
  
  intercept <- coef(lm_model)[1]
  slope <- coef(lm_model)[2]
  r_squared <- model_summary$r.squared
  p_value <- coef(model_summary)[2,4]
  
  regression_results <- rbind(regression_results, data.frame(
    Column = names(intron_sums_heatmap)[i],
    Intercept = intercept,
    Slope = slope,
    R_Squared = r_squared,
    P_Value = p_value))

custom_label <- paste0(
  "y = ", round(intercept, 3), " + ", round(slope, 3), "x",
  "\nR² = ", round(r_squared, 3),
  "\nP = ", format.pval(p_value, digits = 3))

p <- ggplot(regression_data, aes(x = uniqueness_values, y = PhyloP_score)) +
  geom_point(color = "dodgerblue", size = 2) +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  annotate("text", x = max(regression_data$uniqueness_values, na.rm = TRUE) * 0.6, y = max(regression_data$PhyloP_score, na.rm = TRUE) * 0.9, label = custom_label, size = 4, color = "black", hjust = 0) +
  labs(x = paste0("Uniqueness Values (", i, ") for Retained Intron"), y = "Mean PhyloP Score of Retained Intron", title = paste0("Linear Regression: ", i, " Uniqueness Values vs PhyloP Score")) +
  theme_minimal()

print(p)
}
```

# Now do a logistic regression with high, low, and mid uniqueness values:
```{r}
combined_data <- merge(intron_phylop, intron_sums_stats, by = "AS_event_ID") %>%
  select(-shift_or_preserve) %>%
  filter(!is.na(`PhyloP score`)) %>%
  mutate(Max_Bin = factor(Max_Value_Percentile, levels = c("0-20", "21-40", "41-60", "61-80", "81-100")),
         Min_Bin = factor(Min_Value_Percentile, levels = c("0-20", "21-40", "41-60", "61-80", "81-100")))

logit_max <- glm(Max_Bin ~ `PhyloP score`, data = combined_data, family = binomial)
summary_max <- summary(logit_max)
logit_min <- glm(Min_Bin ~ `PhyloP score`, data = combined_data, family = binomial)
summary_min <- summary(logit_min)

p_value_max <- coef(summary_max)[2,4]    # p-value for PhyloP score in logit_max
p_value_min <- coef(summary_min)[2,4]

plot_logit <- function(model, data, xvar, yvar, title) {
  ggplot(data, aes(x = .data[[yvar]], y = .data[[xvar]], color = .data[[yvar]])) +
    geom_jitter(width = 0.35, height = 0.1, size = 1.8) +
    scale_x_discrete(labels = levels(data[[yvar]])) +  # Set x-axis labels to bin labels
    scale_color_manual(values = c("0-20" = "red", "21-40" = "red2", "41-60" = "red3", "61-80" = "maroon", "81-100" = "red4")) +
    labs(x = "Percentile Bin (higher percentile = higher uniqueness value)", y = xvar, title = title, color = "Percentile") +
    theme_minimal()
}

plot_max <- plot_logit(logit_max, combined_data, "PhyloP score", "Max_Bin", "Logistic Regression of IR events: PhyloP Score vs Max Uniqueness Value")
plot_min <- plot_logit(logit_min, combined_data, "PhyloP score", "Min_Bin", "Logistic Regression of IR events: PhyloP Score vs Min Uniqueness Value")

plot_max
plot_min
```

# Now, instead of plotting uniqueness values, plot the individual % usage values of each IR event and correlate them with phylop scores:
```{r}
regression_results <- data.frame(Column = character(),Intercept = numeric(),Slope = numeric(),R_Squared = numeric(),P_Value = numeric(),stringsAsFactors = FALSE)
temp_intron <- mega_intron[match(intron_phylop$AS_event_ID, mega_intron$AS_event_ID),]

for (i in singletypes) {
  usage_values <- temp_intron[[i]]
  usage_values <- as.numeric(gsub("%", "", usage_values))
  regression_data <- data.frame(usage_values, PhyloP_score = intron_phylop$`PhyloP score`)
  regression_data <- regression_data[complete.cases(regression_data),]
  
  lm_model <- lm(PhyloP_score ~ usage_values, data = regression_data)
  model_summary <- summary(lm_model)
  
  intercept <- coef(lm_model)[1]
  slope <- coef(lm_model)[2]
  r_squared <- model_summary$r.squared
  p_value <- coef(model_summary)[2,4]
  
  regression_results <- rbind(regression_results, data.frame(
    Column = names(temp_intron)[i],
    Intercept = intercept,
    Slope = slope,
    R_Squared = r_squared,
    P_Value = p_value))

custom_label <- paste0(
  "y = ", round(intercept, 3), " + ", round(slope, 3), "x",
  "\nR² = ", round(r_squared, 3),
  "\nP = ", format.pval(p_value, digits = 3))

p <- ggplot(regression_data, aes(x = as.numeric(usage_values), y = PhyloP_score)) +
  geom_point(color = "darkblue", size = 2) +
  geom_smooth(method = "lm", color = "red", se = FALSE, na.rm = TRUE) +
  annotate("text", x = 60, y = 5, label = custom_label, size = 4, color = "black", hjust = 0) +
  labs(x = paste0("Percent Usage Values (", i, ") for Retained Intron"), 
       y = "Mean PhyloP Score of Retained Intron", 
       title = paste0("Linear Regression: ", i, " Percent Usage Values vs PhyloP Score")) +
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100)) +
  theme_minimal()

print(p)
}
```

# Sort intron_phylop by highest to lowest phylop score:
```{r}
intron_phylop <- intron_phylop[order(-intron_phylop$`PhyloP score`),]
intron_phylop$Gene <- mega_intron$Gene[match(intron_phylop$AS_event_ID, mega_intron$AS_event_ID)]
intron_phylop <- intron_phylop[,c("AS_event_ID","chr","start","end","Gene","PhyloP score")]
write.csv(intron_phylop, file = "intron_retention_phylop_scores.csv", row.names = F)
```

### Cassette exons instead of retained introns: ###
# Extract start and end coordinates of each cassette exon event:
```{r}
extracted_coords <- lapply(rownames(cassette_sums_heatmap), function(id) {
  event_id <- sub("^[^ ]+ ", "", id)
  parts <- unlist(strsplit(event_id, "_"))
  chr <- as.character(parts[1])
  region_start <- as.numeric(parts[3])
  cassette_start <- as.numeric(parts[4])
  cassette_end <- as.numeric(parts[5])
  region_end <- as.numeric(parts[6])
  return(c(chr, region_start, cassette_start, cassette_end, region_end))
})

extracted_coords <- do.call(rbind, extracted_coords)
rownames(extracted_coords) <- sub("^[^ ]+ ", "", rownames(cassette_sums_heatmap))
#write.csv(extracted_coords, file = "extracted_coords_cassette_events.csv")
```

# calculate mean phylop scores between extracted coordinates ranges
```{r}
# Initialize output data frame
phyloP_scores_combined <- data.frame(AS_event_ID = character(), Chr = character(), region_start = integer(), cassette_start = integer(), cassette_end = integer(), region_end = integer(), Average_Score = numeric(), stringsAsFactors = FALSE)

# Process each row of extracted_coords
for (i in seq_len(nrow(extracted_coords))) {
  AS_event_ID <- rownames(extracted_coords)[i]
  chr <- paste0("chr", extracted_coords[i, 1])
  region_start <- as.integer(extracted_coords[i, 2])
  cassette_start <- as.integer(extracted_coords[i, 3])
  cassette_end <- as.integer(extracted_coords[i, 4])
  region_end <- as.integer(extracted_coords[i, 5])
  
  # Check if start and end are valid
  if (!is.na(region_start) && !is.na(region_end) && region_start < region_end &&
      !is.na(cassette_start) && !is.na(cassette_end) && cassette_start < cassette_end) {
    # Extract phyloP scores from the BigWig file
    Total_Region_Score <- import.bw("ce11.phyloP26way.bw", which = GRanges(chr, IRanges(region_start, region_end)))$score
    Cassette_Exon_Score <- import.bw("ce11.phyloP26way.bw", which = GRanges(chr, IRanges(cassette_start, cassette_end)))$score
    Flanking_Regions_Score_1 <- import.bw("ce11.phyloP26way.bw", which = GRanges(chr, IRanges(region_start, cassette_start)))$score
    Flanking_Regions_Score_2 <- import.bw("ce11.phyloP26way.bw", which = GRanges(chr, IRanges(cassette_end, region_end)))$score
    Flanking_Regions_Score <- c(Flanking_Regions_Score_1, Flanking_Regions_Score_2)
    
    # Compute average scores
    Total_Region_Score <- if (length(Total_Region_Score) > 0) mean(Total_Region_Score, na.rm = TRUE) else NA
    Cassette_Exon_Score <- if (length(Cassette_Exon_Score) > 0) mean(Cassette_Exon_Score, na.rm = TRUE) else NA
    Flanking_Regions_Score <- if (length(Flanking_Regions_Score) > 0) mean(Flanking_Regions_Score, na.rm = TRUE) else NA
    
  # Debugging output
  cat(sprintf("Processing: AS_event_ID=%s, Chr=%s, Start=%d, End=%d, Total Region Score=%f\n", AS_event_ID, chr, region_start, region_end, Total_Region_Score))
  
    # Append result to output data frame
    phyloP_scores_combined <- rbind(phyloP_scores_combined, data.frame(AS_event_ID, chr, region_start, cassette_start, cassette_end, region_end, Total_Region_Score, Cassette_Exon_Score, Flanking_Regions_Score, stringsAsFactors = FALSE))
  } else {
    cat(sprintf("Skipping invalid coordinates: AS_event_ID=%s (%s, %d, %d)\n", AS_event_ID, chr, region_start, region_end))
  }
}
```

# Modify and format phyloP_scores_combined into cassette_phylop:
```{r}
cassette_phylop <- as.data.frame(phyloP_scores_combined)
cassette_phylop$Total_Region_Score <- as.numeric(cassette_phylop$Total_Region_Score)   # convert the score columns to numeric
cassette_phylop$Cassette_Exon_Score <- as.numeric(cassette_phylop$Cassette_Exon_Score)   # convert the score columns to numeric
cassette_phylop$Flanking_Regions_Score <- as.numeric(cassette_phylop$Flanking_Regions_Score)   # convert the score columns to numeric
```

# Correlation of uniqueness index values with PhyloP score:
```{r}
cassette_sums_stats <- cassette_sums_heatmap %>%
  rownames_to_column("AS_event_ID") %>%
  mutate(
    Min_Value = pmin(!!!syms(singletypes), na.rm = TRUE),
    Max_Value = pmax(!!!syms(singletypes), na.rm = TRUE)) %>%
  mutate(
    across(c(Min_Value, Max_Value), 
           ~ cut(.x, 
                 breaks = quantile(.x, probs = seq(0, 1, 0.2), na.rm = TRUE), 
                 labels = c("0-20", "21-40", "41-60", "61-80", "81-100"), 
                 include.lowest = TRUE), 
           .names = "{.col}_Percentile"))

cassette_sums_stats <- cassette_sums_stats %>%
  mutate(AS_event_ID = sub("^[^ ]+ ", "", AS_event_ID))
combined_data <- merge(cassette_phylop, cassette_sums_stats, by.y = "AS_event_ID")
combined_data <- combined_data[!is.na(combined_data$Cassette_Exon_Score),]

create_regression_plot <- function(data, x_var, y_var = "Cassette_Exon_Score") {
  plot_data <- data %>%
    filter(!is.na(.data[[x_var]]), !is.na(.data[[y_var]]))
  
  lm_model <- lm(as.formula(paste(y_var, "~", x_var)), data = plot_data)
  model_summary <- summary(lm_model)
  
  slope <- coef(lm_model)[2]
  intercept <- coef(lm_model)[1]
  r_squared <- model_summary$adj.r.squared
  p_value <- coef(model_summary)[2,4]
  
  custom_label <- paste0(
    "y = ", round(intercept, 3), " + ", round(slope, 3), "x",
    "\nAdj. R² = ", round(r_squared, 3),
    "\nP = ", format.pval(p_value, digits = 3))
  
  ggplot(plot_data, aes(x = .data[[x_var]], y = .data[[y_var]])) +
    geom_point(color = "dodgerblue", size = 2) +
    geom_smooth(method = "lm", color = "red", se = TRUE) +
    annotate("text", 
             x = quantile(plot_data[[x_var]], 0.7, na.rm = TRUE), 
             y = quantile(plot_data[[y_var]], 0.9, na.rm = TRUE), 
             label = custom_label, size = 5, color = "black", hjust = 0) +
    labs(
      x = paste(gsub("_", " ", x_var), "Uniqueness Value"), 
      y = "Cassette Exon PhyloP Score", 
      title = paste("Linear Regression:", gsub("_", " ", x_var), "Uniqueness Value vs Cassette Exon PhyloP Score")
    ) +
    theme_minimal()
}

max_value_plot <- create_regression_plot(combined_data, "Max_Value")
print(max_value_plot)
min_value_plot <- create_regression_plot(combined_data, "Min_Value")
print(min_value_plot)

# Optional: Combine plots if needed
max_value_plot + min_value_plot

## Absolute value of uniqueness index:
cassette_sums_stats <- cassette_sums_heatmap %>%
  rownames_to_column("AS_event_ID") %>%
  mutate(
    Min_Value_Abs = abs(pmin(!!!syms(singletypes), na.rm = TRUE)),
    Max_Value_Abs = abs(pmax(!!!syms(singletypes), na.rm = TRUE)),
    Highest_Abs = pmax(Min_Value_Abs, Max_Value_Abs)) %>%
  mutate(
    across(c(Highest_Abs), 
           ~ cut(.x, 
                 breaks = quantile(.x, probs = seq(0, 1, 0.2), na.rm = TRUE), 
                 labels = c("0-20", "21-40", "41-60", "61-80", "81-100"), 
                 include.lowest = TRUE), 
           .names = "{.col}_Percentile"))

cassette_sums_stats <- cassette_sums_stats %>%
  mutate(AS_event_ID = sub("^[^ ]+ ", "", AS_event_ID))
combined_data <- merge(cassette_phylop, cassette_sums_stats, by.y = "AS_event_ID")
combined_data <- combined_data[!is.na(combined_data$Cassette_Exon_Score),]

create_regression_plot <- function(data, x_var, y_var = "Cassette_Exon_Score") {
  plot_data <- data %>%
    filter(!is.na(.data[[x_var]]), 
           !is.na(.data[[y_var]]),
           is.finite(.data[[x_var]]),
           is.finite(.data[[y_var]]))
  
  if(nrow(plot_data) < 10) {
    warning("Not enough valid data points for regression")
    return(NULL)
  }
  
  lm_model <- lm(as.formula(paste(y_var, "~", x_var)), data = plot_data)
  model_summary <- summary(lm_model)
  
  slope <- coef(lm_model)[2]
  intercept <- coef(lm_model)[1]
  r_squared <- model_summary$adj.r.squared
  p_value <- coef(model_summary)[2,4]
  
  custom_label <- paste0(
    "y = ", round(intercept, 3), " + ", round(slope, 3), "x",
    "\nAdj. R² = ", round(r_squared, 3),
    "\nP = ", format.pval(p_value, digits = 3))
  
  ggplot(plot_data, aes(x = .data[[x_var]], y = .data[[y_var]])) +
    geom_point(color = "dodgerblue", size = 2) +
    geom_smooth(method = "lm", color = "red", se = TRUE) +
    annotate("text", 
             x = quantile(plot_data[[x_var]], 0.7, na.rm = TRUE), 
             y = quantile(plot_data[[y_var]], 0.9, na.rm = TRUE), 
             label = custom_label, size = 5, color = "black", hjust = 0) +
    labs(
      x = paste(gsub("_", " ", x_var), "Uniqueness Value"), 
      y = "Cassette Exon PhyloP Score", 
      title = paste("Linear Regression:", gsub("_", " ", x_var), "Uniqueness Value vs Cassette Exon PhyloP Score")) +
    theme_minimal()
}

# Add additional data cleaning before regression
combined_data <- combined_data %>%
  mutate(
    Highest_Abs = as.numeric(as.character(Highest_Abs)),
    Cassette_Exon_Score = as.numeric(as.character(Cassette_Exon_Score))) %>%
  filter(!is.na(Highest_Abs), 
         !is.na(Cassette_Exon_Score),
         is.finite(Highest_Abs),
         is.finite(Cassette_Exon_Score))

highest_uniqueness_plot <- create_regression_plot(combined_data, "Highest_Abs")
print(highest_uniqueness_plot)
```

# In the last regression, which data points represent cassette exons with high abs. uniqueness values AND high conservation score?
```{r}
# Define thresholds:
uniqueness_threshold <- quantile(combined_data$Highest_Abs, 0.9, na.rm = TRUE)   ### change threshold here
phylop_threshold <- quantile(combined_data$Cassette_Exon_Score, 0.9, na.rm = TRUE)   ### change threshold here

# Filter and arrange
high_uniqueness_high_phylop <- combined_data %>%
  filter(
    Highest_Abs >= uniqueness_threshold,
    Cassette_Exon_Score >= phylop_threshold) %>%
  arrange(desc(Highest_Abs), desc(Cassette_Exon_Score))

high_uniqueness_high_phylop <- high_uniqueness_high_phylop %>%
  dplyr::select(AS_event_ID, Highest_Abs, Cassette_Exon_Score, everything())

# Optional - add gene name for each AS_event_ID:
high_uniqueness_high_phylop <- high_uniqueness_high_phylop %>%
  rowwise() %>%
  mutate(
    Gene = mega_cassette$Gene[mega_cassette$AS_event_ID == .data$AS_event_ID][1]) %>%
  ungroup()

# Optional - reorder columns to insert Gene between Chr and Region_Start
high_uniqueness_high_phylop <- high_uniqueness_high_phylop[, c(
  "AS_event_ID", 
  "Highest_Abs", 
  "Cassette_Exon_Score", 
  "chr", 
  "Gene",
  "region_start", 
  setdiff(names(high_uniqueness_high_phylop), c("AS_event_ID", "Highest_Abs", "Cassette_Exon_Score", "chr", "region_start", "Gene")))]

# Return the number of events meeting the criteria
print(paste("Number of cassette exons above high uniqueness, high phyloP threshold:", nrow(high_uniqueness_high_phylop)))

# Optional - view the top cassette events:
head(high_uniqueness_high_phylop, 20)
# optional - write.csv of results:
write.csv(head(high_uniqueness_high_phylop, 20), file = "cassette_high_uniqueness_high_phylop.csv")

# Optional - create a scatter plot highlighting these events:
highlight_plot <- ggplot(combined_data, aes(x = Highest_Abs, y = Cassette_Exon_Score)) +
  geom_point(color = "gray", alpha = 0.5) +
  geom_point(data = high_uniqueness_high_phylop, color = "red2", size = 2.8) +
  geom_hline(yintercept = phylop_threshold, linetype = "dashed", color = "blue") +
  geom_vline(xintercept = uniqueness_threshold, linetype = "dashed", color = "blue") +
  labs(
    title = "Cassette Exons with High (Abs.) Uniqueness and High PhyloP Score",
    x = "Highest Absolute Uniqueness Value (any cell type)",
    y = "Cassette Exon PhyloP Score") +
  theme_minimal()
print(highlight_plot)
```

# Do the same linear regression, but alter it for unique neuron types:
```{r}
regression_results <- data.frame(Column = character(),Intercept = numeric(),Slope = numeric(),R_Squared = numeric(),P_Value = numeric(),stringsAsFactors = FALSE)

# Create a directory to store plots and initialize empty list to store them in R:
dir.create("uniqueness_regression_plots", showWarnings = FALSE)
plot_list <- list()

for (i in singletypes) {
  uniqueness_values <- as.numeric(cassette_sums_heatmap[[i]])
  regression_data <- data.frame(
    AS_event_ID = rownames(cassette_sums_heatmap),
    uniqueness_values = uniqueness_values,
    stringsAsFactors = FALSE)
  
  regression_data$AS_event_ID <- sub("^[^ ]+ ", "", regression_data$AS_event_ID)
  
  merged_data <- merge(regression_data, cassette_phylop, by = "AS_event_ID", all.x = TRUE)
  regression_data <- merged_data[!is.na(merged_data$Cassette_Exon_Score),]
  
  # eliminate zero-uniqueness values from graph/regression??
  # doesn't really help that much...
  regression_data <- regression_data[regression_data$uniqueness_values != 0,]
  regression_data <- regression_data[regression_data$Cassette_Exon_Score != "No Data",]
  regression_data$Total_Region_Score <- as.numeric(regression_data$Total_Region_Score)
  regression_data$Cassette_Exon_Score <- as.numeric(regression_data$Cassette_Exon_Score)
  regression_data$Flanking_Regions_Score <- as.numeric(regression_data$Flanking_Regions_Score)
  
  lm_model <- lm(Cassette_Exon_Score ~ uniqueness_values, data = regression_data)   ### change phylop score category here
  model_summary <- summary(lm_model)
  
  intercept <- coef(lm_model)[1]
  slope <- coef(lm_model)[2]
  r_squared <- model_summary$r.squared
  p_value <- coef(model_summary)[2,4]
  
  regression_results <- rbind(regression_results, data.frame(
    Column = names(cassette_sums_heatmap)[i],
    Intercept = intercept,
    Slope = slope,
    R_Squared = r_squared,
    P_Value = p_value))

custom_label <- paste0(
  "y = ", round(intercept, 3), " + ", round(slope, 3), "x",
  "\nR² = ", round(r_squared, 3),
  "\nP = ", format.pval(p_value, digits = 3))

p <- ggplot(regression_data, aes(x = uniqueness_values, y = Cassette_Exon_Score)) +
  geom_point(color = "dodgerblue", size = 2) +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  annotate("text", x = max(regression_data$uniqueness_values, na.rm = TRUE) * 0.6, y = max(regression_data$Cassette_Exon_Score, na.rm = TRUE) * 0.9, label = custom_label, size = 4, color = "black", hjust = 0) +
  labs(x = paste0("Uniqueness Values (", i, ") for Cassette Exon"), y = "Mean PhyloP Score of Cassette Exon", title = paste0("Linear Regression: ", i, " Uniqueness Values vs PhyloP Score")) +
  theme_minimal()

print(p)
ggsave(filename = paste0("uniqueness_regression_plots/", i, "_regression_plot.png"), plot = p, width = 10, height = 6)
}

# Create a zip file of the plots:
zip(zipfile = "uniqueness_regression_plots.zip", files = "uniqueness_regression_plots")

# Optional: Remove the directory after zipping
unlink("uniqueness_regression_plots", recursive = TRUE)
```

# Now, instead of plotting uniqueness values, plot the individual % usage values of each cassette exon and correlate them with phylop scores:
```{r}
regression_results <- data.frame(Column = character(),Intercept = numeric(),Slope = numeric(),R_Squared = numeric(),P_Value = numeric(),stringsAsFactors = FALSE)
temp_cassette <- mega_cassette[match(cassette_phylop$AS_event_ID, mega_cassette$AS_event_ID),]

# Create a directory to store plots and initialize empty list to store them in R:
dir.create("usage_regression_plots", showWarnings = FALSE)
plot_list <- list()

for (i in singletypes) {
  usage_values <- temp_cassette[[i]]
  usage_values <- as.numeric(gsub("%", "", usage_values))
  regression_data <- data.frame(usage_values, PhyloP_score = cassette_phylop$Cassette_Exon_Score)    ### change phylop score category here
  regression_data <- regression_data[complete.cases(regression_data),]
  regression_data <- regression_data[regression_data$PhyloP_score != "No Data",]
  regression_data$PhyloP_score <- as.numeric(regression_data$PhyloP_score)
  
  lm_model <- lm(PhyloP_score ~ usage_values, data = regression_data)
  model_summary <- summary(lm_model)
  
  intercept <- coef(lm_model)[1]
  slope <- coef(lm_model)[2]
  r_squared <- model_summary$r.squared
  p_value <- coef(model_summary)[2,4]
  
  regression_results <- rbind(regression_results, data.frame(
    Column = names(temp_cassette)[i],
    Intercept = intercept,
    Slope = slope,
    R_Squared = r_squared,
    P_Value = p_value))

custom_label <- paste0(
  "y = ", round(intercept, 3), " + ", round(slope, 3), "x",
  "\nR² = ", round(r_squared, 3),
  "\nP = ", format.pval(p_value, digits = 3))

p <- ggplot(regression_data, aes(x = as.numeric(usage_values), y = PhyloP_score)) +
  geom_point(color = "darkblue", size = 2) +
  geom_smooth(method = "lm", color = "red", se = FALSE, na.rm = TRUE) +
  annotate("text", x = 60, y = 5, label = custom_label, size = 4, color = "black", hjust = 0) +
  labs(x = paste0("Percent Usage Values (", i, ") for Cassette Exon"), 
       y = "Mean PhyloP Score of Cassette Exon", 
       title = paste0("Linear Regression: ", i, " Percent Usage Values vs PhyloP Score")) +
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100)) +
  theme_minimal()

print(p)
ggsave(filename = paste0("usage_regression_plots/", i, "_regression_plot.png"), plot = p, width = 10, height = 6)
}

# Create a zip file of the plots:
zip(zipfile = "usage_regression_plots.zip", files = "usage_regression_plots")

# Optional: Remove the directory after zipping
unlink("usage_regression_plots", recursive = TRUE)
```

### alternatively, if you already know the coordinates of interest you want to examine (through use of the more generically-named script average_phylop_score_for_regions_of_interest.txt), you can enter them manually in the chunk(s) below:

## implement exact commands from UNIX but in R language ## 3-11:
```{r}
# Define coordinate ranges for exons and introns
regions <- list(
  exon1 = c(8176481, 8176538),
  intron1 = c(8176539, 8177108),
  exon2 = c(8177109, 8177211),
  intron2 = c(8177212, 8177761),
  exon3 = c(8177762, 8177960),
  intron3 = c(8177961, 8178577),
  exon4 = c(8178578, 8178681)
)

# BigWig file path
bigwig_file <- "ce11.phyloP26way.bw"

# Function to calculate average phyloP score for a given range
calculate_average <- function(start, end, bigwig_file) {
  scores <- import.bw(bigwig_file, which = GRanges("chrIV", IRanges(start, end)))$score
  
  if (length(scores) > 0) {
    avg_score <- mean(scores, na.rm = TRUE)
  } else {
    avg_score <- NA
  }
  
  return(avg_score)
}

# Calculate and print the average scores for each region
cat("Average phyloP scores for the defined regions:\n")

for (region in names(regions)) {
  start <- regions[[region]][1]
  end <- regions[[region]][2]
  avg_score <- calculate_average(start, end, bigwig_file)
  
  cat(sprintf("%s (%d-%d): %s\n", region, start, end, ifelse(is.na(avg_score), "No data", round(avg_score, 5))))
}
```

# Load your PhyloP score dataset (pct-1):
```{r}
phylop <- read_tsv("chrIV_8176481_8178681_phyloP_formatted.txt")
phylop <- t(as.data.frame(phylop))

col_name <- "PhyloP score"
colnames(phylop)[1] <- col_name
```

# Declare which scores belong to which features (i.e., exons and introns):
```{r}
first_exon <- data.frame(phylop[as.character(8176481:8176538),])
first_intron <- data.frame(phylop[as.character(8176539:8177108),])
second_exon <- data.frame(phylop[as.character(8177109:8177211),])
second_intron <- data.frame(phylop[as.character(8177212:8177761),])
third_exon <- data.frame(phylop[as.character(8177762:8177960),])
third_intron <- data.frame(phylop[as.character(8177961:8178577),])
fourth_exon <- data.frame(phylop[as.character(8178578:8178681),])

colnames(first_exon) <- col_name
colnames(first_intron) <- col_name
colnames(second_exon) <- col_name
colnames(second_intron) <- col_name
colnames(third_exon) <- col_name
colnames(third_intron) <- col_name
colnames(fourth_exon) <- col_name
```

# List of datasets and corresponding titles:
```{r}
regions <- list(
  "1st Exon" = first_exon,
  "1st Intron" = first_intron,
  "2nd Exon" = second_exon,
  "2nd Intron" = second_intron,
  "3rd Exon" = third_exon,
  "3rd Intron" = third_intron,
  "4th Exon" = fourth_exon)
```

# Violin plot:
```{r}
for (region_name in names(regions)) {
  region_data <- data.frame(regions[[region_name]])
  colnames(region_data) <- "PhyloP score"
  
  sorted_data <- data.frame(region_data[order(-region_data[,1]),])
  colnames(sorted_data) <- "PhyloP score"
  
p <- ggplot(sorted_data, aes(x = region_name, y = `PhyloP score`)) +
  geom_violin(fill = "lightblue", color = "black", trim = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkorange4") +
  geom_hline(yintercept = mean(sorted_data$`PhyloP score`), linetype = "twodash", color = "black") +
  geom_text(aes(x = region_name, y = mean(sorted_data$`PhyloP score`), label = sprintf("Mean: %.2f", mean(sorted_data$`PhyloP score`))), vjust = -0.8, hjust = 0.5, color = "black", size = 3, fontface = "plain") +     ###    mean line label    ###
  labs(y = "PhyloP Scores", x = "Number of Nucleotides", title = paste("PhyloP Scores for Each Nucleotide in the", region_name, "of pct-1")) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.line.x = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold", size = 12),
    axis.title.x = element_text(size = 9, face = "bold"),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"))
  
  print(p)
  ggsave(filename = paste0(tolower(gsub(" ", "_", region_name)), "_violin_pct1.svg"), plot = p, width = 6, height = 4)
}

combined_data <- do.call(rbind, lapply(names(regions), function(region_name) {
  region_data <- data.frame(Region = region_name, PhyloP_score = regions[[region_name]][,1])
  return(region_data)
}))

p <- ggplot(combined_data, aes(x = Region, y = PhyloP_score, fill = Region)) +
  geom_violin(trim = FALSE, color = "black") +
  # Add individual mean lines using stat_summary
  stat_summary(fun = mean, geom = "errorbar", aes(ymin = ..y.., ymax = ..y..),
               color = "black", size = 1.5, width = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkorange4") +
  labs(y = "PhyloP Scores", x = "Regions", title = "PhyloP Scores for Each Nucleotide in pct-1") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.line.x = element_line(color = "black"),
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold", size = 12),
    axis.title.x = element_text(size = 9, face = "bold"),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"))

print(p)
ggsave("all_regions_combined_violin_pct1.svg", plot = p, width = 10, height = 6)
```
