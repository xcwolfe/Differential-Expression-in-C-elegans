## 10th ##

## The difference between this script and the script titled "7_cell_type_vs_cell_type_heatmap_script_for_JUM_analysis.R" is that this script focuses on the alternative splicing counts of cell types independently - we are no longer interested on cell type XXX vs YYY, but now on the abosolute counts of cell type XXX.

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```


```{r}
library("pheatmap")
library("RColorBrewer")
library("tidyverse")
library("data.table")
library("purrr")
library("foreach")
library("parallel")
```

# List your cell types/conditions:
```{r}
singletypes <- c("ADL","AFD","AIM","AIN","AIY","ASEL","ASER","ASG","ASI","ASK","AVA","AVE","AVG","AVH","AVK","AVL","AVM","AWA","AWB","AWC","BAG","CAN","DA","DD","DVC","I5","IL1","IL2","LUA","NSM","OLL","OLQ","PHA","PVC","PVD","PVM","RIA","RIC","RIM","RIS","RMD","SMB","SMD","VB","VC","VD")

singletypes_withreps <- c("ADL1", "ADL2", "ADL3", "ADL4", "AFD1", "AFD2", "AFD3", "AFD4", "AFD5", "AIM1", "AIM2", "AIM3", "AIM4", "AIN1", "AIN2", "AIN3", "AIN4", "AIN5", "AIN6", "AIY1", "AIY2", "AIY3", "ASEL1", "ASEL2", "ASEL3", "ASER1", "ASER2", "ASER3", "ASER4", "ASG1", "ASG2", "ASG3", "ASG4", "ASI1", "ASI2", "ASI3", "ASI4", "ASK1", "ASK2", "ASK3", "ASK4", "AVA1", "AVA2", "AVA3", "AVA4", "AVA5", "AVA6", "AVE1", "AVE2", "AVE3", "AVG1", "AVG2", "AVG3", "AVH1", "AVH2", "AVH3", "AVH4", "AVK1", "AVK2", "AVK3", "AVK4", "AVL1", "AVL2", "AVL3", "AVM1", "AVM2", "AVM3", "AWA1", "AWA2", "AWA3", "AWA4", "AWB1", "AWB2", "AWB3", "AWB4", "AWB5", "AWC1", "AWC2", "AWC3", "AWC4", "BAG1", "BAG2", "BAG3", "BAG4", "CAN1", "CAN2", "CAN3", "DA1", "DA2", "DA3", "DA4", "DD1", "DD2", "DD3", "DD4", "DVC1", "DVC2", "DVC3", "DVC4", "I51", "I52", "I53", "I54", "IL11", "IL12", "IL13", "IL21", "IL22", "IL23", "IL24", "LUA1", "LUA2", "LUA3", "LUA4", "NSM1", "NSM2", "NSM3", "OLL1", "OLL2", "OLQ1", "OLQ2", "OLQ3", "PHA1", "PHA2", "PHA3", "PHA4", "PVC1", "PVC2", "PVC3", "PVC4", "PVC5", "PVD1", "PVD2", "PVM1", "PVM2", "RIA1", "RIA2", "RIA3", "RIA4", "RIA5", "RIA6", "RIC1", "RIC2", "RIC3", "RIC4", "RIM1", "RIM2", "RIM3", "RIM4", "RIS1", "RIS2", "RIS3", "RMD1", "RMD2", "RMD3", "RMD4", "RMD5", "SMB1", "SMB2", "SMB3", "SMB4", "SMB5", "SMD1", "SMD2", "SMD3", "SMD4", "VB1", "VB2", "VB3", "VB4", "VC1", "VC2", "VC3", "VC4", "VC5", "VC6", "VD1", "VD2", "VD3", "VD4")
```

# Use customized functions in order to run future chunk(s):
```{r}
# create a function called "remove_mean_suffix" that takes a column name as input and cuts off the "_mean" part
remove_mean_suffix <- function(column_name) {
  # Find the index of "_mean" in the column name
  mean_index <- regexpr("_mean", column_name)
  
  # If "_mean" is found, extract the substring before it and concatenate with the part after it
  if (mean_index > 0) {
    prefix <- substring(column_name, 1, mean_index - 1)
    suffix <- substring(column_name, mean_index + 6)
    return(paste0(prefix, ".", suffix))
  } else {
    # If "_mean" is not found, return the original column name
    return(column_name)
  }
}

replace_non_finite_with_na <- function(x) {
  x[!is.finite(x)] <- NA
  return(x)
}

replace_na_percent <- function(x) {
  ifelse(x == "NaN%" | x == "NA%", NA, x)
}

replace_inf_percent <- function(x) {
  ifelse(x == "Inf%" | x == "INF%", NA, x)
}

remove_percent <- function(x) {
  as.numeric(gsub("%", "", x))
}

convert_to_numeric <- function(x) {
  as.numeric(sub("%", "", x))
}

convert_percent <- function(x) {
  ifelse(is.na(x), NA, as.numeric(sub("%", "", x)) / 100)
}

get_rows_with_all_na <- function(df) {
  # Extract column names starting with "percentage_usage_mean."
  percentage_cols <- grep("^percentage_usage_mean\\.", colnames(df), value = TRUE)
  
  # Find rows where all percentage_cols contain NA
  all_na_rows <- which(rowSums(is.na(df[, ..percentage_cols])) == length(percentage_cols))
  
  return(all_na_rows)
}
```

# Make a loop of individual cell type-cassette exon results:
```{r}
setwd("E:/Zach Wolfe's JUM analysis/JUM_results")

cassette_indiv_cell_types_list <- vector("list", length(singletypes))
names(cassette_indiv_cell_types_list) <- singletypes

for (i in singletypes) {
  for (j in singletypes) {
    if (i != j) {
      cassette_indiv_cell_types_list[[i]][j] <- vector("list", length(singletypes) - 1)
    }
  }
}

for (i in singletypes){
  for (j in singletypes){
    if (i != j){
    if (j > i){
cassette <- fread(paste0("FINAL_JUM_OUTPUT_pvalue_1", i, "vs", j, "/AS_differential_JUM_output_cassette_exon_events_pvalue_1_final_detailed.txt"))
    }
    if (i > j) next

# Remove the deltaPSI column (for now):
cassette_deltaPSI_col <- cassette$deltaPSI   # Save the column for re-addition later
cassette$deltaPSI <- NULL

# Find the columns with names starting with "percentage_usage."
percentage_usage_cols <- grep("^percentage_usage\\.", colnames(cassette), value = TRUE)

# Create indices for the middle 2 coordinates every 4 rows
first_row_indices <- seq(2, nrow(cassette), by = 4)  # Start from row 2 and increment by 4
last_row_indices <- seq(3, nrow(cassette), by = 4)  # Start from row 3 and increment by 4

for (m in percentage_usage_cols){
  
# Define the new column name
new_column_name <- paste0("percentage_usage_mean.", gsub("percentage_usage\\.", "", m))
  
# Create the new column filled with NAs
cassette[, (new_column_name) := NA_real_]
  
  for (n in 1:length(first_row_indices)){
    
  first_percentage_string <- as.character(cassette[first_row_indices[n],..m])
  last_percentage_string <- as.character(cassette[last_row_indices[n],..m])

  # Remove the percentage symbol using the sub() function and convert to numeric
  first_percentage_numeric <- as.numeric(sub("%", "", first_percentage_string))
  last_percentage_numeric <- as.numeric(sub("%", "", last_percentage_string))

  # Average the result
  percentage_mean <- mean(c(first_percentage_numeric, last_percentage_numeric))

  # Assign the calculated value to the new column in the current row
  cassette[first_row_indices[n], (new_column_name) := percentage_mean]
  cassette[last_row_indices[n], (new_column_name) := percentage_mean]
 }

# Define the column names to reorder
col_to_move <- new_column_name
target_col <- remove_mean_suffix(col_to_move)

# Get the current column order
current_cols <- colnames(cassette)

# Find the index of the target column
target_index <- which(current_cols == target_col)

# Reorder the columns
new_col_order <- c(current_cols[1:target_index], col_to_move, current_cols[(target_index + 2):length(current_cols) -1])
cassette <- cassette[, ..new_col_order]
}

# The very last new_column_name and target_column will be duplicated; remove these columns for simplicity:
cassette <- subset(cassette, select = 1:(ncol(cassette)-2))

# We can re-add our deltaPSI column:
cassette$deltaPSI <- cassette_deltaPSI_col

# And finally we can add columns that represent the means of the percentage_usage_means for each cell type:
i_cols <- grep(paste0("^percentage_usage_mean.", i), colnames(cassette), value = TRUE)
i_vals <- subset(cassette, select = i_cols)
i_vals[] <- lapply(i_vals, replace_non_finite_with_na)
i_vals$mean <- as.numeric()
j_cols <- grep(paste0("^percentage_usage_mean.", j), colnames(cassette), value = TRUE)
j_vals <- subset(cassette, select = j_cols)
j_vals[] <- lapply(j_vals, replace_non_finite_with_na)
j_vals$mean <- as.numeric()
for (k in first_row_indices){
i_vals[k,"mean"] <- mean(as.numeric(i_vals[k,-"mean"]), na.rm = T)
j_vals[k,"mean"] <- mean(as.numeric(j_vals[k,-"mean"]), na.rm = T)
}

# Define the column names for i and j
col_name_i <- paste0("average_", i)
col_name_j <- paste0("average_", j)

setDT(cassette)

# Modify mean columns
  cassette[, (col_name_i)] <- i_vals$mean
  cassette[, (col_name_j)] <- j_vals$mean

for (p in i_cols){
  cassette[[p]] <- paste0(cassette[[p]], "%")
}

for (q in j_cols){
  cassette[[q]] <- paste0(cassette[[q]], "%")
}

average_cols <- grep("^average_", colnames(cassette), value = TRUE)

for (p in average_cols){
  cassette[[p]] <- paste0(cassette[[p]], "%")
}

cassette[] <- lapply(cassette, replace_na_percent)
cassette[] <- lapply(cassette, replace_inf_percent)

# Save as a .csv:
write.csv(cassette, file = paste0(i, "_", j, "_detailed_cassette_R.csv"), row.names = F, na = "")
write.csv(cassette, file = paste0(j, "_", i, "_detailed_cassette_R.csv"), row.names = F, na = "")
    }
    if (i == j) next
  cassette_indiv_cell_types_list[[i]][[j]] <- cassette
  cassette_indiv_cell_types_list[[j]][[i]] <- cassette
  }
}
```

# Now we need to combine all cassette exon events from each cell type so that we get an event ID, gene name, and average % spliced in for each cell type:
```{r}
setwd("E:/Zach Wolfe's JUM analysis/JUM_results")

cassette_indiv_usage_combined <- vector("list", length(singletypes))
names(cassette_indiv_usage_combined) <- singletypes

cassette_indiv_usage_by_rep <- vector("list", length(singletypes_withreps))
names(cassette_indiv_usage_by_rep) <- singletypes_withreps

# rbind all files for cell type XXX:
for (i in singletypes) {
  cassette_indiv_usage_combined[[i]] <- data.table()
  cassette_indiv_usage_by_rep[[i]] <- data.table()

  matches <- grep(paste0("^", i), singletypes_withreps, value = TRUE)
  
  for (j in singletypes) {
    if (i != j) {
      cassette <- cassette_indiv_cell_types_list[[i]][[j]]
      cassette <- cassette[first_row_indices, ]
      cassette_by_rep <- cassette[, c("Gene", "AS_event_ID", paste0("percentage_usage_mean.", matches)), with = FALSE]
      cassette <- cassette[, c("Gene", "AS_event_ID", paste0("average_", i)), with = FALSE]
      
      # Use rbindlist to combine data tables with different numbers of columns
      cassette_indiv_usage_combined[[i]] <- rbindlist(list(cassette_indiv_usage_combined[[i]], cassette), fill = TRUE)
      cassette_indiv_usage_combined[[i]] <- unique(cassette_indiv_usage_combined[[i]])
      
      cassette_indiv_usage_by_rep[[i]] <- rbindlist(list(cassette_indiv_usage_by_rep[[i]], cassette_by_rep), fill = TRUE)
      cassette_indiv_usage_by_rep[[i]] <- unique(cassette_indiv_usage_by_rep[[i]])
    }
    if (i == j) next
  }
  write.csv(cassette_indiv_usage_combined[[i]], file = paste0("combined_cassette_exon_usage_", i, ".csv"), row.names = F)
  write.csv(cassette_indiv_usage_by_rep[[i]], file = paste0("combined_cassette_exon_usage_by_rep_", i, ".csv"), row.names = F)
}

# Lastly, remove the sublists which are labeled by individual replicate (while keeping the sublists representing whole cell types)
cassette_indiv_usage_by_rep <- Filter(function(x) !is.null(x), cassette_indiv_usage_by_rep)
```

# Make a heatmap for all A.S. events represented by % usage for each cell type:
```{r}
mega_cassette <- data.frame(Gene = character(), AS_event_ID = character())
mega_cassette_by_rep <- data.frame(Gene = character(), AS_event_ID = character())

for (i in rev(singletypes)) {
cassette <- cassette_indiv_usage_combined[[i]]
cassette <- as.data.frame(cassette)
mega_cassette <- merge(cassette, mega_cassette, by = c("Gene", "AS_event_ID"), all = T)
unique_rows <- !duplicated(mega_cassette$AS_event_ID, fromLast = TRUE) & !duplicated(mega_cassette$AS_event_ID)
mega_cassette <- mega_cassette[unique_rows,]
rownames(mega_cassette) <- paste(mega_cassette$Gene, mega_cassette$AS_event_ID) # You can label both the gene name and the AS_event_ID as rownames, or edit this line to label rownames as just one or the other
}

for (i in rev(singletypes)) {
cassette <- cassette_indiv_usage_by_rep[[i]]
cassette <- as.data.frame(cassette)
mega_cassette_by_rep <- merge(cassette, mega_cassette_by_rep, by = c("Gene", "AS_event_ID"), all.x = T, all.y = T)
unique_rows <- !duplicated(mega_cassette_by_rep$AS_event_ID, fromLast = TRUE) & !duplicated(mega_cassette_by_rep$AS_event_ID)
mega_cassette_by_rep <- mega_cassette_by_rep[unique_rows,]
rownames(mega_cassette_by_rep) <- paste(mega_cassette_by_rep$Gene, mega_cassette_by_rep$AS_event_ID) # You can label both the gene name and the AS_event_ID as rownames, or edit this line to label rownames as just one or the other
}

mega_cassette <- mega_cassette[,c(3:ncol(mega_cassette))]
colnames(mega_cassette) <- gsub("^average_", "", colnames(mega_cassette))

mega_cassette_by_rep <- mega_cassette_by_rep[,c(3:ncol(mega_cassette_by_rep))]
colnames(mega_cassette_by_rep) <- gsub("^percentage_usage_mean.", "", colnames(mega_cassette_by_rep))

mega_cassette_matrix <- apply(mega_cassette, 2, convert_to_numeric)
rownames(mega_cassette_matrix) <- rownames(mega_cassette)

mega_cassette_matrix_by_rep <- apply(mega_cassette_by_rep, 2, convert_to_numeric)
rownames(mega_cassette_matrix_by_rep) <- rownames(mega_cassette_by_rep)

colors <- colorRampPalette(brewer.pal(6, "Greens"))(250)
pheatmap(mega_cassette_matrix,
         na_col = "#FFFFFF",
         border_color = "black",
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = F,
         number_format = "%.0f",
         number_color = "black",
         fontsize = 7,
         fontsize_number = 6.5,
         show_rownames = T,
         show_colnames = T,
         angle_col = 90,
         col = colors,
         main = "Percent Usage of Cassette Exons by Cell Type")

pheatmap(mega_cassette_matrix_by_rep,
         na_col = "#FFFFFF",
         border_color = "black",
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = F,
         number_format = "%.0f",
         number_color = "black",
         fontsize = 7,
         fontsize_number = 6.5,
         show_rownames = T,
         show_colnames = T,
         angle_col = 90,
         col = colors,
         main = "Percent Usage of Cassette Exons by Replicate")

write.csv(mega_cassette_matrix, file = "mega_cassette_matrix.csv")
write.csv(mega_cassette_matrix_by_rep, file = "mega_cassette_matrix_by_rep.csv")
```


# Individual cell type-intron retention results:
```{r}
setwd("E:/Zach Wolfe's JUM analysis/JUM_results")

intron_indiv_cell_types_list <- vector("list", length(singletypes))
names(intron_indiv_cell_types_list) <- singletypes

for (i in singletypes) {
  for (j in singletypes) {
    if (i != j) {
      intron_indiv_cell_types_list[[i]][j] <- vector("list", length(singletypes) - 1)
    }
  }
}

for (i in singletypes){
  for (j in singletypes){
    if (i != j){
    if (j > i){
intron <- fread(paste0("FINAL_JUM_OUTPUT_pvalue_1", i, "vs", j, "/AS_differential_JUM_output_intron_retention_pvalue_1_final_detailed.txt"))
    }
    if (i > j) next

# Remove the deltaPSI column (for now):
intron_deltaPSI_col <- intron$deltaPSI   # Save the column for re-addition later
intron$deltaPSI <- NULL

# Find the columns with names starting with "percentage_usage."
percentage_usage_cols <- grep("^percentage_usage\\.", colnames(intron), value = TRUE)

# Create indices for the middle 2 coordinates every 4 rows
first_row_indices <- seq(2, nrow(intron), by = 4)  # Start from row 2 and increment by 4
last_row_indices <- seq(4, nrow(intron), by = 4)  # Start from row 4 and increment by 4

for (m in percentage_usage_cols){
  
# Define the new column name
new_column_name <- paste0("percentage_usage_mean.", gsub("percentage_usage\\.", "", m))
  
# Create the new column filled with NAs
intron[, (new_column_name) := NA_real_]
  
  for (n in 1:length(first_row_indices)){
    
  first_percentage_string <- as.character(intron[first_row_indices[n],..m])
  last_percentage_string <- as.character(intron[last_row_indices[n],..m])

  # Remove the percentage symbol using the sub() function and convert to numeric
  first_percentage_numeric <- as.numeric(sub("%", "", first_percentage_string))
  last_percentage_numeric <- as.numeric(sub("%", "", last_percentage_string))

  # Average the result
  percentage_mean <- mean(c(first_percentage_numeric, last_percentage_numeric))

  # Assign the calculated value to the new column in the current row
  intron[first_row_indices[n], (new_column_name) := percentage_mean]
  intron[last_row_indices[n], (new_column_name) := percentage_mean]
 }

# Define the column names to reorder
col_to_move <- new_column_name
target_col <- remove_mean_suffix(col_to_move)

# Get the current column order
current_cols <- colnames(intron)

# Find the index of the target column
target_index <- which(current_cols == target_col)

# Reorder the columns
new_col_order <- c(current_cols[1:target_index], col_to_move, current_cols[(target_index + 2):length(current_cols) -1])
intron <- intron[, ..new_col_order]
}

# The very last new_column_name and target_column will be duplicated; remove these columns for simplicity:
intron <- subset(intron, select = 1:(ncol(intron)-2))

# We can re-add our deltaPSI column:
intron$deltaPSI <- intron_deltaPSI_col

# And finally we can add columns that represent the means of the percentage_usage_means for each cell type:
i_cols <- grep(paste0("^percentage_usage_mean.", i), colnames(intron), value = TRUE)
i_vals <- subset(intron, select = i_cols)
i_vals[] <- lapply(i_vals, replace_non_finite_with_na)
i_vals$mean <- as.numeric()
j_cols <- grep(paste0("^percentage_usage_mean.", j), colnames(intron), value = TRUE)
j_vals <- subset(intron, select = j_cols)
j_vals[] <- lapply(j_vals, replace_non_finite_with_na)
j_vals$mean <- as.numeric()
for (k in first_row_indices){
i_vals[k,"mean"] <- mean(as.numeric(i_vals[k,-"mean"]), na.rm = T)
j_vals[k,"mean"] <- mean(as.numeric(j_vals[k,-"mean"]), na.rm = T)
}

# Define the column names for i and j
col_name_i <- paste0("average_", i)
col_name_j <- paste0("average_", j)

# Modify mean columns
  intron[, (col_name_i)] <- i_vals$mean
  intron[, (col_name_j)] <- j_vals$mean

for (p in i_cols){
  intron[[p]] <- paste0(intron[[p]], "%")
}

for (q in j_cols){
  intron[[q]] <- paste0(intron[[q]], "%")
}

average_cols <- grep("^average_", colnames(intron), value = TRUE)

for (p in average_cols){
  intron[[p]] <- paste0(intron[[p]], "%")
}

intron[] <- lapply(intron, replace_na_percent)
intron[] <- lapply(intron, replace_inf_percent)

# Save as a .csv:
write.csv(intron, file = paste0(i, "_", j, "_detailed_intron_R.csv"), row.names = F, na = "")
write.csv(intron, file = paste0(j, "_", i, "_detailed_intron_R.csv"), row.names = F, na = "")
    }
    if (i == j) next
  intron_indiv_cell_types_list[[i]][[j]] <- intron
  intron_indiv_cell_types_list[[j]][[i]] <- intron
  }
}
```

# Now we need to combine all intron retention events from each cell type so that we get an event ID, gene name, and average % spliced in for each cell type:
```{r}
setwd("E:/Zach Wolfe's JUM analysis/JUM_results")

intron_indiv_usage_combined <- vector("list", length(singletypes))
names(intron_indiv_usage_combined) <- singletypes

intron_indiv_usage_by_rep <- vector("list", length(singletypes_withreps))
names(intron_indiv_usage_by_rep) <- singletypes_withreps

# rbind all files for cell type XXX:
for (i in singletypes) {
  intron_indiv_usage_combined[[i]] <- data.table()
  intron_indiv_usage_by_rep[[i]] <- data.table()

  matches <- grep(paste0("^", i), singletypes_withreps, value = TRUE)
  
  for (j in singletypes) {
    if (i != j) {
      intron <- intron_indiv_cell_types_list[[i]][[j]]
      intron <- intron[first_row_indices, ]
      intron_by_rep <- intron[, c("Gene", "AS_event_ID", paste0("percentage_usage_mean.", matches)), with = FALSE]
      intron <- intron[, c("Gene", "AS_event_ID", paste0("average_", i)), with = FALSE]
      
      # Use rbindlist to combine data tables with different numbers of columns
      intron_indiv_usage_combined[[i]] <- rbindlist(list(intron_indiv_usage_combined[[i]], intron), fill = TRUE)
      intron_indiv_usage_combined[[i]] <- unique(intron_indiv_usage_combined[[i]])
      
      intron_indiv_usage_by_rep[[i]] <- rbindlist(list(intron_indiv_usage_by_rep[[i]], intron_by_rep), fill = TRUE)
      intron_indiv_usage_by_rep[[i]] <- unique(intron_indiv_usage_by_rep[[i]])
    }
    if (i == j) next
  }
  
  # Calculate the row means for columns starting with "average_"
  #avg_cols <- grep("^average_", colnames(intron_indiv_usage_by_rep[[i]]), value = TRUE)
  intron_indiv_usage_combined[[i]] <- intron_indiv_usage_combined[[i]] %>%
    mutate(across(all_of(avg_cols), ~ remove_percent(.) %>% mean(na.rm = TRUE), .names = "avg_{.col}"),
           across(starts_with("avg_"), ~ paste0(., "%")))
  
  write.csv(intron_indiv_usage_combined[[i]], file = paste0("combined_intron_retention_usage_", i, ".csv"), row.names = F)
  write.csv(intron_indiv_usage_by_rep[[i]], file = paste0("combined_intron_retention_usage_by_rep_", i, ".csv"), row.names = F)
}

# Lastly, remove the sublists which are labeled by individual replicate (while keeping the sublists representing whole cell types)
intron_indiv_usage_by_rep <- Filter(function(x) !is.null(x), intron_indiv_usage_by_rep)
```

# Make a heatmap for all A.S. events represented by % usage for each cell type:
```{r}
mega_intron <- data.frame(Gene = character(), AS_event_ID = character())
mega_intron_by_rep <- data.frame(Gene = character(), AS_event_ID = character())

for (i in rev(singletypes)) {
intron <- intron_indiv_usage_combined[[i]]
intron <- as.data.frame(intron)
intron[,3] <- as.numeric(sub("%", "", intron[,3])) / 100
intron_mean <- aggregate(. ~ Gene + AS_event_ID, data = intron, FUN = function(x) mean(as.numeric(x)))
intron_mean[,3] <- paste0(intron_mean[,3] * 100, "%")
mega_intron <- merge(intron_mean, mega_intron, by = c("Gene", "AS_event_ID"), all = TRUE)
rownames(mega_intron) <- paste(mega_intron$Gene, mega_intron$AS_event_ID) # You can label both the gene name and the AS_event_ID as rownames, or edit this line to label rownames as just one or the other
}

for (i in singletypes) {
intron <- intron_indiv_usage_by_rep[[i]]
intron <- as.data.frame(intron)
 numeric_cols <- 3:ncol(intron)
 intron[, numeric_cols] <- lapply(intron[, numeric_cols], function(x) {
   ifelse(grepl("%", x, fixed = TRUE), sub("%", "", x), x)
 })
intron_mean <- aggregate(. ~ Gene + AS_event_ID, data = intron, FUN = function(x) mean(as.numeric(x)), na.action = NULL) ### 8-10 you HAVE to add an NA ignore parameter here
mega_intron_by_rep <- full_join(mega_intron_by_rep, intron_mean, by = c("Gene", "AS_event_ID"))
mega_intron_by_rep <- mega_intron_by_rep %>% arrange(Gene)
rownames(mega_intron_by_rep) <- paste(mega_intron_by_rep$Gene, mega_intron_by_rep$AS_event_ID) # You can label both the gene name and the AS_event_ID as rownames, or edit this line to label rownames as just one or the other
}

mega_intron <- mega_intron[,c(3:ncol(mega_intron))]
colnames(mega_intron) <- gsub("^average_", "", colnames(mega_intron))

mega_intron_by_rep <- mega_intron_by_rep[,c(3:ncol(mega_intron_by_rep))]
colnames(mega_intron_by_rep) <- gsub("^percentage_usage_mean.", "", colnames(mega_intron_by_rep))

mega_intron_matrix <- apply(mega_intron, 2, convert_to_numeric)
rownames(mega_intron_matrix) <- rownames(mega_intron)

mega_intron_matrix_by_rep <- apply(mega_intron_by_rep, 2, convert_to_numeric)
rownames(mega_intron_matrix_by_rep) <- rownames(mega_intron_by_rep)

colors <- colorRampPalette(brewer.pal(6, "YlOrBr"))(250)
pheatmap(mega_intron_matrix,
         na_col = "#FFFFFF",
         border_color = "black",
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = F,
         number_format = "%.0f",
         number_color = "black",
         fontsize = 7,
         fontsize_number = 6.5,
         show_rownames = T,
         show_colnames = T,
         angle_col = 90,
         col = colors,
         main = "Percent Usage of Intron Retention Events by Cell Type")

pheatmap(mega_intron_matrix_by_rep,
         na_col = "#FFFFFF",
         border_color = "black",
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = F,
         number_format = "%.0f",
         number_color = "black",
         fontsize = 7,
         fontsize_number = 6.5,
         show_rownames = T,
         show_colnames = T,
         angle_col = 90,
         col = colors,
         main = "Percent Usage of Intron Retention Events by Replicate")

write.csv(mega_intron_matrix, file = "mega_intron_matrix.csv")
write.csv(mega_intron_matrix_by_rep, file = "mega_intron_matrix_by_rep.csv")
```

# Individual cell type-MXE results:
```{r}
setwd("E:/Zach Wolfe's JUM analysis/JUM_results")

MXE_indiv_cell_types_list <- vector("list", length(singletypes))
names(MXE_indiv_cell_types_list) <- singletypes

for (i in singletypes) {
  for (j in singletypes) {
    if (i != j) {
      MXE_indiv_cell_types_list[[i]][j] <- vector("list", length(singletypes) - 1)
    }
  }
}

for (i in singletypes){
  for (j in singletypes){
    if (i != j){
    if (j > i){
MXE <- fread(paste0("FINAL_JUM_OUTPUT_pvalue_1", i, "vs", j, "/AS_differential_JUM_output_MXE_events_pvalue_1_final_detailed.txt"))
    }
    if (i > j) next

if (nrow(MXE) != 0){

# Find the columns with names starting with "percentage_usage."
percentage_usage_cols <- grep("^percentage_usage\\.", colnames(MXE), value = TRUE)

MXE_plus <- MXE[MXE$sub_junction_strand == "+"]
MXE_minus <- MXE[MXE$sub_junction_strand == "-"]

# Remove the deltaPSI columns (for now):
MXE_deltaPSI_col_plus <- MXE_plus$deltaPSI   # Save the column for re-addition later
MXE_plus$deltaPSI <- NULL

# Remove the deltaPSI columns (for now):
MXE_deltaPSI_col_minus <- MXE_minus$deltaPSI   # Save the column for re-addition later
MXE_minus$deltaPSI <- NULL

if (nrow(MXE_plus) > 0){
# Create indices for the 1st and 3rd coordinates every 4 rows
first_row_indices <- seq(1, nrow(MXE_plus), by = 4)  # Start from row 2 and increment by 4
last_row_indices <- seq(3, nrow(MXE_plus), by = 4)  # Start from row 4 and increment by 4

for (m in percentage_usage_cols){
  
# Define the new column name
new_column_name <- paste0("percentage_usage_mean.", gsub("percentage_usage\\.", "", m))
  
# Create the new column filled with NAs
MXE_plus[, (new_column_name) := NA_real_]
  
  for (n in 1:length(first_row_indices)){
    
  first_percentage_string <- as.character(MXE_plus[first_row_indices[n],..m])
  last_percentage_string <- as.character(MXE_plus[last_row_indices[n],..m])

  # Remove the percentage symbol using the sub() function and convert to numeric
  first_percentage_numeric <- as.numeric(sub("%", "", first_percentage_string))
  last_percentage_numeric <- as.numeric(sub("%", "", last_percentage_string))

  # Average the result
  percentage_mean <- mean(c(first_percentage_numeric, last_percentage_numeric))

  # Assign the calculated value to the new column in the current row
  MXE_plus[first_row_indices[n], (new_column_name) := percentage_mean]
  MXE_plus[last_row_indices[n], (new_column_name) := percentage_mean]
 }

# Define the column names to reorder
col_to_move <- new_column_name
target_col <- remove_mean_suffix(col_to_move)

# Get the current column order
current_cols <- colnames(MXE_plus)

# Find the index of the target column
target_index <- which(current_cols == target_col)

# Reorder the columns
new_col_order <- c(current_cols[1:target_index], col_to_move, current_cols[(target_index + 2):length(current_cols) -1])
MXE_plus <- MXE_plus[, ..new_col_order]
}

# The very last new_column_name and target_column will be duplicated; remove these columns for simplicity:
MXE_plus <- subset(MXE_plus, select = 1:(ncol(MXE_plus)-2))

# We can re-add our deltaPSI column:
MXE_plus$deltaPSI <- MXE_deltaPSI_col_plus

# And finally we can add columns that represent the means of the percentage_usage_means for each cell type (so long as there is 1 or more MXE event):
i_cols <- grep(paste0("^percentage_usage_mean.", i), colnames(MXE_plus), value = TRUE)
i_vals <- subset(MXE_plus, select = i_cols)
i_vals[] <- lapply(i_vals, replace_non_finite_with_na)
i_vals$mean <- as.numeric()
j_cols <- grep(paste0("^percentage_usage_mean.", j), colnames(MXE_plus), value = TRUE)
j_vals <- subset(MXE_plus, select = j_cols)
j_vals[] <- lapply(j_vals, replace_non_finite_with_na)
j_vals$mean <- as.numeric()
for (k in first_row_indices){
i_vals[k,"mean"] <- mean(as.numeric(i_vals[k,-"mean"]), na.rm = T)
j_vals[k,"mean"] <- mean(as.numeric(j_vals[k,-"mean"]), na.rm = T)
}

# Define the column names for i and j
col_name_i <- paste0("average_", i)
col_name_j <- paste0("average_", j)

# Add or modify columns
MXE_plus <- MXE_plus %>%
  mutate(!!col_name_i := i_vals$mean,
         !!col_name_j := j_vals$mean)

}
if (nrow(MXE_plus) == 0) next
    
if (nrow(MXE_minus) > 0){
# Create indices for the 2nd and 4th coordinates every 4 rows
first_row_indices <- seq(2, nrow(MXE_minus), by = 4)  # Start from row 1 and increment by 4
last_row_indices <- seq(4, nrow(MXE_minus), by = 4)  # Start from row 3 and increment by 4

for (m in percentage_usage_cols){
  
# Define the new column name
new_column_name <- paste0("percentage_usage_mean.", gsub("percentage_usage\\.", "", m))
  
# Create the new column filled with NAs
MXE_minus[, (new_column_name) := NA_real_]
  
  for (n in 1:length(first_row_indices)){
    
  first_percentage_string <- as.character(MXE_minus[first_row_indices[n],..m])
  last_percentage_string <- as.character(MXE_minus[last_row_indices[n],..m])

  # Remove the percentage symbol using the sub() function and convert to numeric
  first_percentage_numeric <- as.numeric(sub("%", "", first_percentage_string))
  last_percentage_numeric <- as.numeric(sub("%", "", last_percentage_string))

  # Average the result
  percentage_mean <- mean(c(first_percentage_numeric, last_percentage_numeric))

  # Assign the calculated value to the new column in the current row
  MXE_minus[first_row_indices[n], (new_column_name) := percentage_mean]
  MXE_minus[last_row_indices[n], (new_column_name) := percentage_mean]
 }

# Define the column names to reorder
col_to_move <- new_column_name
target_col <- remove_mean_suffix(col_to_move)

# Get the current column order
current_cols <- colnames(MXE_minus)

# Find the index of the target column
target_index <- which(current_cols == target_col)

# Reorder the columns
new_col_order <- c(current_cols[1:target_index], col_to_move, current_cols[(target_index + 2):length(current_cols) -1])
MXE_minus <- MXE_minus[, ..new_col_order]
 }

MXE_minus <- subset(MXE_minus, select = 1:(ncol(MXE_minus)-2))

# We can re-add our deltaPSI column:
MXE_minus$deltaPSI <- MXE_deltaPSI_col_minus

# And finally we can add columns that represent the means of the percentage_usage_means for each cell type (so long as there is 1 or more MXE event):
i_cols <- grep(paste0("^percentage_usage_mean.", i), colnames(MXE_minus), value = TRUE)
i_vals <- subset(MXE_minus, select = i_cols)
i_vals[] <- lapply(i_vals, replace_non_finite_with_na)
i_vals$mean <- as.numeric()
j_cols <- grep(paste0("^percentage_usage_mean.", j), colnames(MXE_minus), value = TRUE)
j_vals <- subset(MXE_minus, select = j_cols)
j_vals[] <- lapply(j_vals, replace_non_finite_with_na)
j_vals$mean <- as.numeric()
for (k in first_row_indices){
i_vals[k,"mean"] <- mean(as.numeric(i_vals[k,-"mean"]), na.rm = T)
j_vals[k,"mean"] <- mean(as.numeric(j_vals[k,-"mean"]), na.rm = T)
}

# Define the column names for i and j
col_name_i <- paste0("average_", i)
col_name_j <- paste0("average_", j)

# Add or modify columns
MXE_minus <- MXE_minus %>%
  mutate(!!col_name_i := i_vals$mean,
         !!col_name_j := j_vals$mean)

}
if (nrow(MXE_minus) == 0) next

MXE <- rbind(MXE_plus, MXE_minus)
  
for (p in i_cols){
  MXE[[p]] <- paste0(MXE[[p]], "%")
}

for (q in j_cols){
  MXE[[q]] <- paste0(MXE[[q]], "%")
}

average_cols <- grep("^average_", colnames(MXE), value = TRUE)

for (p in average_cols){
  MXE[[p]] <- paste0(MXE[[p]], "%")
}

MXE[] <- lapply(MXE, replace_na_percent)
MXE[] <- lapply(MXE, replace_inf_percent)
    
}
  if (nrow(MXE) == 0) next
# Save as a .csv:
write.csv(MXE, file = paste0(i, "_", j, "_detailed_MXE_R.csv"), row.names = F, na = "")
write.csv(MXE, file = paste0(j, "_", i, "_detailed_MXE_R.csv"), row.names = F, na = "")
    }
    if (i == j) next
  MXE_indiv_cell_types_list[[i]][[j]] <- MXE
  MXE_indiv_cell_types_list[[j]][[i]] <- MXE
  }
}
```

# Now we need to combine all MXE events from each cell type so that we get an event ID, gene name, and average % spliced in for each cell type:
```{r}
setwd("E:/Zach Wolfe's JUM analysis/JUM_results")

MXE_indiv_usage_combined <- vector("list", length(singletypes))
names(MXE_indiv_usage_combined) <- singletypes

MXE_indiv_usage_by_rep <- vector("list", length(singletypes_withreps))
names(MXE_indiv_usage_by_rep) <- singletypes_withreps

# rbind all files for cell type XXX:
for (i in singletypes) {
  MXE_indiv_usage_combined[[i]] <- data.table()
  MXE_indiv_usage_by_rep[[i]] <- data.table()

  matches <- grep(paste0("^", i), singletypes_withreps, value = TRUE)
  
  for (j in singletypes) {
    if (i != j) {
      MXE <- MXE_indiv_cell_types_list[[i]][[j]]
      MXE <- MXE[!is.na(MXE[[paste0("average_", i)]]),]
      MXE_by_rep <- MXE[, c("Gene", "AS_event_ID", paste0("percentage_usage_mean.", matches)), with = FALSE]
      MXE <- MXE[, c("Gene", "AS_event_ID", paste0("average_", i)), with = FALSE]
      
      # Use rbindlist to combine data tables with different numbers of columns
      MXE_indiv_usage_combined[[i]] <- rbindlist(list(MXE_indiv_usage_combined[[i]], MXE), fill = TRUE)
      MXE_indiv_usage_combined[[i]] <- unique(MXE_indiv_usage_combined[[i]])
      
      MXE_indiv_usage_by_rep[[i]] <- rbindlist(list(MXE_indiv_usage_by_rep[[i]], MXE_by_rep), fill = TRUE)
      MXE_indiv_usage_by_rep[[i]] <- unique(MXE_indiv_usage_by_rep[[i]])
    }
    if (i == j) next
  }
  
  # Calculate the row means for columns starting with "average_"
  #avg_cols <- grep("^average_", colnames(MXE_indiv_usage_by_rep[[i]]), value = TRUE)
  MXE_indiv_usage_combined[[i]] <- MXE_indiv_usage_combined[[i]] %>%
    mutate(across(all_of(avg_cols), ~ remove_percent(.) %>% mean(na.rm = TRUE), .names = "avg_{.col}"),
           across(starts_with("avg_"), ~ paste0(., "%")))
  
  write.csv(MXE_indiv_usage_combined[[i]], file = paste0("combined_MXE_usage_", i, ".csv"), row.names = F)
  write.csv(MXE_indiv_usage_by_rep[[i]], file = paste0("combined_MXE_usage_by_rep_", i, ".csv"), row.names = F)
}

# Lastly, remove the sublists which are labeled by individual replicate (while keeping the sublists representing whole cell types)
MXE_indiv_usage_by_rep <- Filter(function(x) !is.null(x), MXE_indiv_usage_by_rep)
```

# Make a heatmap for all A.S. events represented by % usage for each cell type:
```{r}
mega_MXE <- data.frame(Gene = character(), AS_event_ID = character())
mega_MXE_by_rep <- data.frame(Gene = character(), AS_event_ID = character())

for (i in rev(singletypes)) {
MXE <- MXE_indiv_usage_combined[[i]]
MXE <- as.data.frame(MXE)
MXE[,3] <- as.numeric(sub("%", "", MXE[,3])) / 100
MXE_mean <- aggregate(. ~ Gene + AS_event_ID, data = MXE, FUN = function(x) mean(as.numeric(x)), na.action = NULL)
MXE_mean[,3] <- paste0(MXE_mean[,3] * 100, "%")
mega_MXE <- merge(MXE_mean, mega_MXE, by = c("Gene", "AS_event_ID"), all = TRUE)
rownames(mega_MXE) <- paste(mega_MXE$Gene, mega_MXE$AS_event_ID) # You can label both the gene name and the AS_event_ID as rownames, or edit this line to label rownames as just one or the other
}

for (i in singletypes) {
MXE <- MXE_indiv_usage_by_rep[[i]]
MXE <- as.data.frame(MXE)
 numeric_cols <- 3:ncol(MXE)
 MXE[, numeric_cols] <- lapply(MXE[, numeric_cols], function(x) {
   ifelse(grepl("%", x, fixed = TRUE), sub("%", "", x), x)
 })
MXE_mean <- aggregate(. ~ Gene + AS_event_ID, data = MXE, FUN = function(x) mean(as.numeric(x)), na.action = NULL)
mega_MXE_by_rep <- full_join(mega_MXE_by_rep, MXE_mean, by = c("Gene", "AS_event_ID"))
mega_MXE_by_rep <- mega_MXE_by_rep %>% arrange(tolower(Gene))
rownames(mega_MXE_by_rep) <- paste(mega_MXE_by_rep$Gene, mega_MXE_by_rep$AS_event_ID) # You can label both the gene name and the AS_event_ID as rownames, or edit this line to label rownames as just one or the other
}

mega_MXE <- mega_MXE[,c(3:ncol(mega_MXE))]
colnames(mega_MXE) <- gsub("^average_", "", colnames(mega_MXE))

mega_MXE_by_rep <- mega_MXE_by_rep[,c(3:ncol(mega_MXE_by_rep))]
colnames(mega_MXE_by_rep) <- gsub("^percentage_usage_mean.", "", colnames(mega_MXE_by_rep))

mega_MXE_matrix <- apply(mega_MXE, 2, convert_to_numeric)
rownames(mega_MXE_matrix) <- rownames(mega_MXE)

mega_MXE_matrix_by_rep <- apply(mega_MXE_by_rep, 2, convert_to_numeric)
rownames(mega_MXE_matrix_by_rep) <- rownames(mega_MXE_by_rep)

colors <- colorRampPalette(brewer.pal(6, "PuRd"))(250)
pheatmap(mega_MXE_matrix,
         na_col = "#FFFFFF",
         border_color = "black",
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = F,
         number_format = "%.0f",
         number_color = "black",
         fontsize = 7,
         fontsize_number = 6.5,
         show_rownames = T,
         show_colnames = T,
         angle_col = 90,
         col = colors,
         main = "Percent Usage of MXE Events by Cell Type")

pheatmap(mega_MXE_matrix_by_rep,
         na_col = "#FFFFFF",
         border_color = "black",
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = F,
         number_format = "%.0f",
         number_color = "black",
         fontsize = 7,
         fontsize_number = 6.5,
         show_rownames = T,
         show_colnames = T,
         angle_col = 90,
         col = colors,
         main = "Percent Usage of MXE Events by Replicate")

write.csv(mega_MXE_matrix, file = "mega_MXE_matrix.csv")
write.csv(mega_MXE_matrix_by_rep, file = "mega_MXE_matrix_by_rep.csv")
```

# Individual cell type-A3S results:
```{r}
setwd("E:/Zach Wolfe's JUM analysis/JUM_results")

A3S_indiv_cell_types_list <- vector("list", length(singletypes))
names(A3S_indiv_cell_types_list) <- singletypes

for (i in singletypes) {
  for (j in singletypes) {
    if (i != j) {
      A3S_indiv_cell_types_list[[i]][j] <- vector("list", length(singletypes) - 1)
    }
  }
}

for (i in singletypes){
  for (j in singletypes){
    if (i != j){
    if (j > i){
A3S <- fread(paste0("FINAL_JUM_OUTPUT_pvalue_1", i, "vs", j, "/AS_differential_JUM_output_A3SS_events_pvalue_1_final_detailed.txt"))
    }
    if (i > j) next

if (nrow(A3S) != 0){

# Find the columns with names starting with "percentage_usage."
percentage_usage_cols <- grep("^percentage_usage\\.", colnames(A3S), value = TRUE)

A3S_plus <- A3S[A3S$sub_junction_strand == "+"]
A3S_minus <- A3S[A3S$sub_junction_strand == "-"]

# add a function to eliminate all A3S values that include more than two 3' splice site options:
A3S_plus <- A3S_plus %>%
  group_by(AS_event_ID) %>%
  filter(n() < 3) %>%
  ungroup()

A3S_minus <- A3S_minus %>%
  group_by(AS_event_ID) %>%
  filter(n() < 3) %>%
  ungroup()

# Remove the deltaPSI columns (for now):
A3S_deltaPSI_col_plus <- A3S_plus$deltaPSI   # Save the column for re-addition later
A3S_plus$deltaPSI <- NULL

# Remove the deltaPSI columns (for now):
A3S_deltaPSI_col_minus <- A3S_minus$deltaPSI   # Save the column for re-addition later
A3S_minus$deltaPSI <- NULL

if (nrow(A3S_plus) > 0){
# Create indices for the 1st coordinate every 2 rows
row_indices <- seq(1, nrow(A3S_plus), by = 2)  # Start from row 1 and increment by 2

for (m in percentage_usage_cols){
  
# Define the new column name
new_column_name <- paste0("percentage_usage_mean.", gsub("percentage_usage\\.", "", m))
  
# Create the new column filled with NAs
A3S_plus[, new_column_name] <- rep(NA_real_, nrow(A3S_plus))
  
  for (n in 1:length(row_indices)){
    
  percentage_string <- as.character(A3S_plus[row_indices[n],m])

  # Remove the percentage symbol using the sub() function and convert to numeric
  percentage_numeric <- as.numeric(sub("%", "", percentage_string))

  # Assign the calculated value to the new column in the current row
  A3S_plus[row_indices[n], new_column_name] <- percentage_numeric
 }

# Define the column names to reorder
col_to_move <- new_column_name
target_col <- remove_mean_suffix(col_to_move)

# Get the current column order
current_cols <- colnames(A3S_plus)

# Find the index of the target column
target_index <- which(current_cols == target_col)

# Reorder the columns
new_col_order <- c(current_cols[1:target_index], col_to_move, current_cols[(target_index + 2):length(current_cols) -1])
A3S_plus <- A3S_plus[, new_col_order]
}

# The very last new_column_name and target_column will be duplicated; remove these columns for simplicity:
A3S_plus <- subset(A3S_plus, select = 1:(ncol(A3S_plus)-2))

# We can re-add our deltaPSI column:
A3S_plus$deltaPSI <- A3S_deltaPSI_col_plus

# And finally we can add columns that represent the means of the percentage_usage_means for each cell type (so long as there is 1 or more A3S event):
i_cols <- grep(paste0("^percentage_usage_mean.", i), colnames(A3S_plus), value = TRUE)
i_vals <- subset(A3S_plus, select = i_cols)
i_vals[] <- lapply(i_vals, replace_non_finite_with_na)
i_vals <- as_tibble(i_vals)
i_vals <- i_vals %>%
  mutate(mean = numeric(n()))
j_cols <- grep(paste0("^percentage_usage_mean.", j), colnames(A3S_plus), value = TRUE)
j_vals <- subset(A3S_plus, select = j_cols)
j_vals[] <- lapply(j_vals, replace_non_finite_with_na)
j_vals <- as_tibble(j_vals)
j_vals <- j_vals %>%
  mutate(mean = numeric(n()))
for (k in 1:nrow(i_vals)) {
  row_values <- as.numeric(i_vals[k, -which(names(i_vals) == "mean")])
  if (all(is.na(row_values)) == TRUE) {
    i_vals[k, "mean"] <- NA_real_
  } else {
    i_vals[k, "mean"] <- mean(row_values, na.rm = TRUE)
  }
  row_values <- as.numeric(j_vals[k, -which(names(j_vals) == "mean")])
  if (all(is.na(row_values))== TRUE) {
    j_vals[k, "mean"] <- NA_real_
  } else {
    j_vals[k, "mean"] <- mean(row_values, na.rm = TRUE)
  }
}

# Define the column names for i and j
col_name_i <- paste0("average_", i)
col_name_j <- paste0("average_", j)

# Add or modify columns
A3S_plus <- A3S_plus %>%
  mutate(!!col_name_i := i_vals$mean,
         !!col_name_j := j_vals$mean)

}
if (nrow(A3S_plus) == 0) next
    
if (nrow(A3S_minus) > 0){
# Create indices for the 2nd coordinate every 2 rows
row_indices <- seq(2, nrow(A3S_minus), by = 2)  # Start from row 2 and increment by 2

for (m in percentage_usage_cols){
  
# Define the new column name
new_column_name <- paste0("percentage_usage_mean.", gsub("percentage_usage\\.", "", m))
  
# Create the new column filled with NAs
A3S_minus[, new_column_name] <- rep(NA_real_, nrow(A3S_minus))

  for (n in 1:length(row_indices)){
    
  percentage_string <- as.character(A3S_minus[row_indices[n],m])

  # Remove the percentage symbol using the sub() function and convert to numeric
  percentage_numeric <- as.numeric(sub("%", "", percentage_string))

  # Assign the calculated value to the new column in the current row
  A3S_minus[row_indices[n], new_column_name] <- percentage_numeric
 }

# Define the column names to reorder
col_to_move <- new_column_name
target_col <- remove_mean_suffix(col_to_move)

# Get the current column order
current_cols <- colnames(A3S_minus)

# Find the index of the target column
target_index <- which(current_cols == target_col)

# Reorder the columns
new_col_order <- c(current_cols[1:target_index], col_to_move, current_cols[(target_index + 2):length(current_cols) -1])
A3S_minus <- A3S_minus[, new_col_order]
 }

A3S_minus <- subset(A3S_minus, select = 1:(ncol(A3S_minus)-2))

# We can re-add our deltaPSI column:
A3S_minus$deltaPSI <- A3S_deltaPSI_col_minus

# And finally we can add columns that represent the percentage_usage_means for each cell type (so long as there is 1 or more A3S event):
i_cols <- grep(paste0("^percentage_usage_mean.", i), colnames(A3S_minus), value = TRUE)
i_vals <- subset(A3S_minus, select = i_cols)
i_vals[] <- lapply(i_vals, replace_non_finite_with_na)
i_vals <- as_tibble(i_vals)
i_vals <- i_vals %>%
  mutate(mean = numeric(n()))
j_cols <- grep(paste0("^percentage_usage_mean.", j), colnames(A3S_minus), value = TRUE)
j_vals <- subset(A3S_minus, select = j_cols)
j_vals[] <- lapply(j_vals, replace_non_finite_with_na)
j_vals <- as_tibble(j_vals)
j_vals <- j_vals %>%
  mutate(mean = numeric(n()))
for (k in 1:nrow(i_vals)) {
  row_values <- as.numeric(i_vals[k, -which(names(i_vals) == "mean")])
  if (all(is.na(row_values)) == TRUE) {
    i_vals[k, "mean"] <- NA_real_
  } else {
    i_vals[k, "mean"] <- mean(row_values, na.rm = TRUE)
  }
  row_values <- as.numeric(j_vals[k, -which(names(j_vals) == "mean")])
  if (all(is.na(row_values))== TRUE) {
    j_vals[k, "mean"] <- NA_real_
  } else {
    j_vals[k, "mean"] <- mean(row_values, na.rm = TRUE)
  }
}

# Define the column names for i and j
col_name_i <- paste0("average_", i)
col_name_j <- paste0("average_", j)

# Add or modify columns
A3S_minus <- A3S_minus %>%
  mutate(!!col_name_i := i_vals$mean,
         !!col_name_j := j_vals$mean)

}
if (nrow(A3S_minus) == 0) next

A3S <- rbind(A3S_plus, A3S_minus)
  
for (p in i_cols){
  A3S[[p]] <- paste0(A3S[[p]], "%")
}

for (q in j_cols){
  A3S[[q]] <- paste0(A3S[[q]], "%")
}

average_cols <- grep("^average_", colnames(A3S), value = TRUE)

for (p in average_cols){
  A3S[[p]] <- paste0(A3S[[p]], "%")
}

A3S[] <- lapply(A3S, replace_na_percent)
A3S[] <- lapply(A3S, replace_inf_percent)
    
}
  if (nrow(A3S) == 0) next
# Save as a .csv:
write.csv(A3S, file = paste0(i, "_", j, "_detailed_A3S_R.csv"), row.names = F, na = "")
write.csv(A3S, file = paste0(j, "_", i, "_detailed_A3S_R.csv"), row.names = F, na = "")
    }
    if (i == j) next
  A3S_indiv_cell_types_list[[i]][[j]] <- A3S
  A3S_indiv_cell_types_list[[j]][[i]] <- A3S
  }
}
```

# Now we need to combine all A3S events from each cell type so that we get an event ID, gene name, and average % spliced in for each cell type:
```{r}
setwd("E:/Zach Wolfe's JUM analysis/JUM_results")

A3S_indiv_usage_combined <- vector("list", length(singletypes))
names(A3S_indiv_usage_combined) <- singletypes

A3S_indiv_usage_by_rep <- vector("list", length(singletypes_withreps))
names(A3S_indiv_usage_by_rep) <- singletypes_withreps

# rbind all files for cell type XXX:
for (i in singletypes) {
  A3S_indiv_usage_combined[[i]] <- data.table()
  A3S_indiv_usage_by_rep[[i]] <- data.table()

  matches <- grep(paste0("^", i), singletypes_withreps, value = TRUE)
  
  for (j in singletypes) {
    if (i != j) {
      A3S <- A3S_indiv_cell_types_list[[i]][[j]]
      A3S <- A3S[!is.na(A3S[[paste0("average_", i)]]),]
      A3S_by_rep <- A3S[, c("Gene", "AS_event_ID", paste0("percentage_usage_mean.", matches)), with = FALSE]
      A3S <- A3S[, c("Gene", "AS_event_ID", paste0("average_", i)), with = FALSE]
      
      # Use rbindlist to combine data tables with different numbers of columns
      A3S_indiv_usage_combined[[i]] <- rbindlist(list(A3S_indiv_usage_combined[[i]], A3S), fill = TRUE)
      A3S_indiv_usage_combined[[i]] <- unique(A3S_indiv_usage_combined[[i]])
      
      A3S_indiv_usage_by_rep[[i]] <- rbindlist(list(A3S_indiv_usage_by_rep[[i]], A3S_by_rep), fill = TRUE)
      A3S_indiv_usage_by_rep[[i]] <- unique(A3S_indiv_usage_by_rep[[i]])
    }
    if (i == j) next
  }

  A3S_indiv_usage_combined[[i]] <- A3S_indiv_usage_combined[[i]][rowSums(is.na(A3S_indiv_usage_combined[[i]])) != ncol(A3S_indiv_usage_combined[[i]]) - 2, ]
  A3S_indiv_usage_by_rep[[i]] <- A3S_indiv_usage_by_rep[[i]][rowSums(is.na(A3S_indiv_usage_by_rep[[i]])) != ncol(A3S_indiv_usage_by_rep[[i]]) - 2, ]
  
  write.csv(A3S_indiv_usage_combined[[i]], file = paste0("combined_A3S_usage_", i, ".csv"), row.names = F)
  write.csv(A3S_indiv_usage_by_rep[[i]], file = paste0("combined_A3S_usage_by_rep_", i, ".csv"), row.names = F)
}

# Lastly, remove the sublists which are labeled by individual replicate (while keeping the sublists representing whole cell types)
A3S_indiv_usage_by_rep <- Filter(function(x) !is.null(x), A3S_indiv_usage_by_rep)
```

# Make a heatmap for all A.S. events represented by % usage for each cell type:
```{r}
mega_A3S <- data.frame(Gene = character(), AS_event_ID = character())
mega_A3S_by_rep <- data.frame(Gene = character(), AS_event_ID = character())

for (i in rev(singletypes)) {
A3S <- A3S_indiv_usage_combined[[i]]
A3S <- as.data.frame(A3S)
A3S[,3] <- as.numeric(sub("%", "", A3S[,3])) / 100
A3S_mean <- aggregate(. ~ Gene + AS_event_ID, data = A3S, FUN = function(x) mean(as.numeric(x)), na.action = NULL)
A3S_mean[,3] <- paste0(A3S_mean[,3] * 100, "%")
mega_A3S <- merge(A3S_mean, mega_A3S, by = c("Gene", "AS_event_ID"), all = TRUE)
rownames(mega_A3S) <- paste(mega_A3S$Gene, mega_A3S$AS_event_ID) # You can label both the gene name and the AS_event_ID as rownames, or edit this line to label rownames as just one or the other
}

for (i in singletypes) {
A3S <- A3S_indiv_usage_by_rep[[i]]
A3S <- as.data.frame(A3S)
 numeric_cols <- 3:ncol(A3S)
 A3S[, numeric_cols] <- lapply(A3S[, numeric_cols], function(x) {
   ifelse(grepl("%", x, fixed = TRUE), sub("%", "", x), x)
 })
A3S_mean <- aggregate(. ~ Gene + AS_event_ID, data = A3S, FUN = function(x) mean(as.numeric(x)), na.action = NULL)
mega_A3S_by_rep <- full_join(mega_A3S_by_rep, A3S_mean, by = c("Gene", "AS_event_ID"))
mega_A3S_by_rep <- mega_A3S_by_rep %>% arrange(tolower(Gene))
rownames(mega_A3S_by_rep) <- paste(mega_A3S_by_rep$Gene, mega_A3S_by_rep$AS_event_ID) # You can label both the gene name and the AS_event_ID as rownames, or edit this line to label rownames as just one or the other
}

mega_A3S <- mega_A3S[,c(3:ncol(mega_A3S))]
colnames(mega_A3S) <- gsub("^average_", "", colnames(mega_A3S))

mega_A3S_by_rep <- mega_A3S_by_rep[,c(3:ncol(mega_A3S_by_rep))]
colnames(mega_A3S_by_rep) <- gsub("^percentage_usage_mean.", "", colnames(mega_A3S_by_rep))

mega_A3S_matrix <- apply(mega_A3S, 2, convert_to_numeric)
rownames(mega_A3S_matrix) <- rownames(mega_A3S)

mega_A3S_matrix_by_rep <- apply(mega_A3S_by_rep, 2, convert_to_numeric)
rownames(mega_A3S_matrix_by_rep) <- rownames(mega_A3S_by_rep)

mega_A3S_matrix <- mega_A3S_matrix[rowSums(!is.na(mega_A3S_matrix)) > 0,]

mega_A3S_matrix_by_rep <- mega_A3S_matrix_by_rep[rowSums(!is.na(mega_A3S_matrix_by_rep)) > 0,]

colors <- colorRampPalette(brewer.pal(6, "Blues"))(250)
pheatmap(mega_A3S_matrix,
         na_col = "#FFFFFF",
         border_color = "black",
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = F,
         number_format = "%.0f",
         number_color = "black",
         fontsize = 7,
         fontsize_number = 6.5,
         show_rownames = T,
         show_colnames = T,
         angle_col = 90,
         col = colors,
         main = "Percent Usage of A3S Events by Cell Type")

pheatmap(mega_A3S_matrix_by_rep,
         na_col = "#FFFFFF",
         border_color = "black",
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = F,
         number_format = "%.0f",
         number_color = "black",
         fontsize = 7,
         fontsize_number = 6.5,
         show_rownames = T,
         show_colnames = T,
         angle_col = 90,
         col = colors,
         main = "Percent Usage of A3S Events by Replicate")

write.csv(mega_A3S_matrix, file = "mega_A3S_matrix.csv")
write.csv(mega_A3S_matrix_by_rep, file = "mega_A3S_matrix_by_rep.csv")
```

# Individual cell type-A5S results:
```{r}
setwd("E:/Zach Wolfe's JUM analysis/JUM_results")

A5S_indiv_cell_types_list <- vector("list", length(singletypes))
names(A5S_indiv_cell_types_list) <- singletypes

for (i in singletypes) {
  for (j in singletypes) {
    if (i != j) {
      A5S_indiv_cell_types_list[[i]][j] <- vector("list", length(singletypes) - 1)
    }
  }
}

for (i in singletypes){
  for (j in singletypes){
    if (i != j){
    if (j > i){
A5S <- fread(paste0("FINAL_JUM_OUTPUT_pvalue_1", i, "vs", j, "/AS_differential_JUM_output_A5SS_events_pvalue_1_final_detailed.txt"))
    }
    if (i > j) next

if (nrow(A5S) != 0){

# Find the columns with names starting with "percentage_usage."
percentage_usage_cols <- grep("^percentage_usage\\.", colnames(A5S), value = TRUE)

A5S_plus <- A5S[A5S$sub_junction_strand == "+"]
A5S_minus <- A5S[A5S$sub_junction_strand == "-"]

# add a function to eliminate all A5S values that include more than two 3' splice site options:
A5S_plus <- A5S_plus %>%
  group_by(AS_event_ID) %>%
  filter(n() < 3) %>%
  ungroup()

A5S_minus <- A5S_minus %>%
  group_by(AS_event_ID) %>%
  filter(n() < 3) %>%
  ungroup()

# Remove the deltaPSI columns (for now):
A5S_deltaPSI_col_plus <- A5S_plus$deltaPSI   # Save the column for re-addition later
A5S_plus$deltaPSI <- NULL

# Remove the deltaPSI columns (for now):
A5S_deltaPSI_col_minus <- A5S_minus$deltaPSI   # Save the column for re-addition later
A5S_minus$deltaPSI <- NULL

if (nrow(A5S_plus) > 0){
# Create indices for the 1st coordinate every 2 rows
row_indices <- seq(1, nrow(A5S_plus), by = 2)  # Start from row 1 and increment by 2

for (m in percentage_usage_cols){
  
# Define the new column name
new_column_name <- paste0("percentage_usage_mean.", gsub("percentage_usage\\.", "", m))
  
# Create the new column filled with NAs
A5S_plus[, new_column_name] <- rep(NA_real_, nrow(A5S_plus))
  
  for (n in 1:length(row_indices)){
    
  percentage_string <- as.character(A5S_plus[row_indices[n],m])

  # Remove the percentage symbol using the sub() function and convert to numeric
  percentage_numeric <- as.numeric(sub("%", "", percentage_string))

  # Assign the calculated value to the new column in the current row
  A5S_plus[row_indices[n], new_column_name] <- percentage_numeric
 }

# Define the column names to reorder
col_to_move <- new_column_name
target_col <- remove_mean_suffix(col_to_move)

# Get the current column order
current_cols <- colnames(A5S_plus)

# Find the index of the target column
target_index <- which(current_cols == target_col)

# Reorder the columns
new_col_order <- c(current_cols[1:target_index], col_to_move, current_cols[(target_index + 2):length(current_cols) -1])
A5S_plus <- A5S_plus[, new_col_order]
}

# The very last new_column_name and target_column will be duplicated; remove these columns for simplicity:
A5S_plus <- subset(A5S_plus, select = 1:(ncol(A5S_plus)-2))

# We can re-add our deltaPSI column:
A5S_plus$deltaPSI <- A5S_deltaPSI_col_plus

# And finally we can add columns that represent the means of the percentage_usage_means for each cell type (so long as there is 1 or more A5S event):
i_cols <- grep(paste0("^percentage_usage_mean.", i), colnames(A5S_plus), value = TRUE)
i_vals <- subset(A5S_plus, select = i_cols)
i_vals[] <- lapply(i_vals, replace_non_finite_with_na)
i_vals <- as_tibble(i_vals)
i_vals <- i_vals %>%
  mutate(mean = numeric(n()))
j_cols <- grep(paste0("^percentage_usage_mean.", j), colnames(A5S_plus), value = TRUE)
j_vals <- subset(A5S_plus, select = j_cols)
j_vals[] <- lapply(j_vals, replace_non_finite_with_na)
j_vals <- as_tibble(j_vals)
j_vals <- j_vals %>%
  mutate(mean = numeric(n()))
for (k in 1:nrow(i_vals)) {
  row_values <- as.numeric(i_vals[k, -which(names(i_vals) == "mean")])
  if (all(is.na(row_values)) == TRUE) {
    i_vals[k, "mean"] <- NA_real_
  } else {
    i_vals[k, "mean"] <- mean(row_values, na.rm = TRUE)
  }
  row_values <- as.numeric(j_vals[k, -which(names(j_vals) == "mean")])
  if (all(is.na(row_values))== TRUE) {
    j_vals[k, "mean"] <- NA_real_
  } else {
    j_vals[k, "mean"] <- mean(row_values, na.rm = TRUE)
  }
}

# Define the column names for i and j
col_name_i <- paste0("average_", i)
col_name_j <- paste0("average_", j)

# Add or modify columns
A5S_plus <- A5S_plus %>%
  mutate(!!col_name_i := i_vals$mean,
         !!col_name_j := j_vals$mean)

}
if (nrow(A5S_plus) == 0) next
    
if (nrow(A5S_minus) > 0){
# Create indices for the 2nd coordinate every 2 rows
row_indices <- seq(2, nrow(A5S_minus), by = 2)  # Start from row 2 and increment by 2

for (m in percentage_usage_cols){
  
# Define the new column name
new_column_name <- paste0("percentage_usage_mean.", gsub("percentage_usage\\.", "", m))
  
# Create the new column filled with NAs
A5S_minus[, new_column_name] <- rep(NA_real_, nrow(A5S_minus))

  for (n in 1:length(row_indices)){
    
  percentage_string <- as.character(A5S_minus[row_indices[n],m])

  # Remove the percentage symbol using the sub() function and convert to numeric
  percentage_numeric <- as.numeric(sub("%", "", percentage_string))

  # Assign the calculated value to the new column in the current row
  A5S_minus[row_indices[n], new_column_name] <- percentage_numeric
 }

# Define the column names to reorder
col_to_move <- new_column_name
target_col <- remove_mean_suffix(col_to_move)

# Get the current column order
current_cols <- colnames(A5S_minus)

# Find the index of the target column
target_index <- which(current_cols == target_col)

# Reorder the columns
new_col_order <- c(current_cols[1:target_index], col_to_move, current_cols[(target_index + 2):length(current_cols) -1])
A5S_minus <- A5S_minus[, new_col_order]
 }

A5S_minus <- subset(A5S_minus, select = 1:(ncol(A5S_minus)-2))

# We can re-add our deltaPSI column:
A5S_minus$deltaPSI <- A5S_deltaPSI_col_minus

# And finally we can add columns that represent the percentage_usage_means for each cell type (so long as there is 1 or more A5S event):
i_cols <- grep(paste0("^percentage_usage_mean.", i), colnames(A5S_minus), value = TRUE)
i_vals <- subset(A5S_minus, select = i_cols)
i_vals[] <- lapply(i_vals, replace_non_finite_with_na)
i_vals <- as_tibble(i_vals)
i_vals <- i_vals %>%
  mutate(mean = numeric(n()))
j_cols <- grep(paste0("^percentage_usage_mean.", j), colnames(A5S_minus), value = TRUE)
j_vals <- subset(A5S_minus, select = j_cols)
j_vals[] <- lapply(j_vals, replace_non_finite_with_na)
j_vals <- as_tibble(j_vals)
j_vals <- j_vals %>%
  mutate(mean = numeric(n()))
for (k in 1:nrow(i_vals)) {
  row_values <- as.numeric(i_vals[k, -which(names(i_vals) == "mean")])
  if (all(is.na(row_values)) == TRUE) {
    i_vals[k, "mean"] <- NA_real_
  } else {
    i_vals[k, "mean"] <- mean(row_values, na.rm = TRUE)
  }
  row_values <- as.numeric(j_vals[k, -which(names(j_vals) == "mean")])
  if (all(is.na(row_values))== TRUE) {
    j_vals[k, "mean"] <- NA_real_
  } else {
    j_vals[k, "mean"] <- mean(row_values, na.rm = TRUE)
  }
}

# Define the column names for i and j
col_name_i <- paste0("average_", i)
col_name_j <- paste0("average_", j)

# Add or modify columns
A5S_minus <- A5S_minus %>%
  mutate(!!col_name_i := i_vals$mean,
         !!col_name_j := j_vals$mean)

}
if (nrow(A5S_minus) == 0) next

A5S <- rbind(A5S_plus, A5S_minus)
  
for (p in i_cols){
  A5S[[p]] <- paste0(A5S[[p]], "%")
}

for (q in j_cols){
  A5S[[q]] <- paste0(A5S[[q]], "%")
}

average_cols <- grep("^average_", colnames(A5S), value = TRUE)

for (p in average_cols){
  A5S[[p]] <- paste0(A5S[[p]], "%")
}

A5S[] <- lapply(A5S, replace_na_percent)
A5S[] <- lapply(A5S, replace_inf_percent)
    
}
  if (nrow(A5S) == 0) next
# Save as a .csv:
write.csv(A5S, file = paste0(i, "_", j, "_detailed_A5S_R.csv"), row.names = F, na = "")
write.csv(A5S, file = paste0(j, "_", i, "_detailed_A5S_R.csv"), row.names = F, na = "")
    }
    if (i == j) next
  A5S_indiv_cell_types_list[[i]][[j]] <- A5S
  A5S_indiv_cell_types_list[[j]][[i]] <- A5S
  }
}
```

# Now we need to combine all A5S events from each cell type so that we get an event ID, gene name, and average % spliced in for each cell type:
```{r}
setwd("E:/Zach Wolfe's JUM analysis/JUM_results")

A5S_indiv_usage_combined <- vector("list", length(singletypes))
names(A5S_indiv_usage_combined) <- singletypes

A5S_indiv_usage_by_rep <- vector("list", length(singletypes_withreps))
names(A5S_indiv_usage_by_rep) <- singletypes_withreps

# rbind all files for cell type XXX:
for (i in singletypes) {
  A5S_indiv_usage_combined[[i]] <- data.table()
  A5S_indiv_usage_by_rep[[i]] <- data.table()

  matches <- grep(paste0("^", i), singletypes_withreps, value = TRUE)
  
  for (j in singletypes) {
    if (i != j) {
      A5S <- A5S_indiv_cell_types_list[[i]][[j]]
      A5S <- A5S[!is.na(A5S[[paste0("average_", i)]]),]
      A5S_by_rep <- A5S[, c("Gene", "AS_event_ID", paste0("percentage_usage_mean.", matches)), with = FALSE]
      A5S <- A5S[, c("Gene", "AS_event_ID", paste0("average_", i)), with = FALSE]
      
      # Use rbindlist to combine data tables with different numbers of columns
      A5S_indiv_usage_combined[[i]] <- rbindlist(list(A5S_indiv_usage_combined[[i]], A5S), fill = TRUE)
      A5S_indiv_usage_combined[[i]] <- unique(A5S_indiv_usage_combined[[i]])
      
      A5S_indiv_usage_by_rep[[i]] <- rbindlist(list(A5S_indiv_usage_by_rep[[i]], A5S_by_rep), fill = TRUE)
      A5S_indiv_usage_by_rep[[i]] <- unique(A5S_indiv_usage_by_rep[[i]])
    }
    if (i == j) next
  }

  A5S_indiv_usage_combined[[i]] <- A5S_indiv_usage_combined[[i]][rowSums(is.na(A5S_indiv_usage_combined[[i]])) != ncol(A5S_indiv_usage_combined[[i]]) - 2, ]
  A5S_indiv_usage_by_rep[[i]] <- A5S_indiv_usage_by_rep[[i]][rowSums(is.na(A5S_indiv_usage_by_rep[[i]])) != ncol(A5S_indiv_usage_by_rep[[i]]) - 2, ]
  
  write.csv(A5S_indiv_usage_combined[[i]], file = paste0("combined_A5S_usage_", i, ".csv"), row.names = F)
  write.csv(A5S_indiv_usage_by_rep[[i]], file = paste0("combined_A5S_usage_by_rep_", i, ".csv"), row.names = F)
}

# Lastly, remove the sublists which are labeled by individual replicate (while keeping the sublists representing whole cell types)
A5S_indiv_usage_by_rep <- Filter(function(x) !is.null(x), A5S_indiv_usage_by_rep)
```

# Make a heatmap for all A.S. events represented by % usage for each cell type:
```{r}
mega_A5S <- data.frame(Gene = character(), AS_event_ID = character())
mega_A5S_by_rep <- data.frame(Gene = character(), AS_event_ID = character())

for (i in rev(singletypes)) {
A5S <- A5S_indiv_usage_combined[[i]]
A5S <- as.data.frame(A5S)
A5S[,3] <- as.numeric(sub("%", "", A5S[,3])) / 100
A5S_mean <- aggregate(. ~ Gene + AS_event_ID, data = A5S, FUN = function(x) mean(as.numeric(x)), na.action = NULL)
A5S_mean[,3] <- paste0(A5S_mean[,3] * 100, "%")
mega_A5S <- merge(A5S_mean, mega_A5S, by = c("Gene", "AS_event_ID"), all = TRUE)
rownames(mega_A5S) <- paste(mega_A5S$Gene, mega_A5S$AS_event_ID) # You can label both the gene name and the AS_event_ID as rownames, or edit this line to label rownames as just one or the other
}

for (i in singletypes) {
A5S <- A5S_indiv_usage_by_rep[[i]]
A5S <- as.data.frame(A5S)
 numeric_cols <- 3:ncol(A5S)
 A5S[, numeric_cols] <- lapply(A5S[, numeric_cols], function(x) {
   ifelse(grepl("%", x, fixed = TRUE), sub("%", "", x), x)
 })
A5S_mean <- aggregate(. ~ Gene + AS_event_ID, data = A5S, FUN = function(x) mean(as.numeric(x)), na.action = NULL)
mega_A5S_by_rep <- full_join(mega_A5S_by_rep, A5S_mean, by = c("Gene", "AS_event_ID"))
mega_A5S_by_rep <- mega_A5S_by_rep %>% arrange(tolower(Gene))
rownames(mega_A5S_by_rep) <- paste(mega_A5S_by_rep$Gene, mega_A5S_by_rep$AS_event_ID) # You can label both the gene name and the AS_event_ID as rownames, or edit this line to label rownames as just one or the other
}

mega_A5S <- mega_A5S[,c(3:ncol(mega_A5S))]
colnames(mega_A5S) <- gsub("^average_", "", colnames(mega_A5S))

mega_A5S_by_rep <- mega_A5S_by_rep[,c(3:ncol(mega_A5S_by_rep))]
colnames(mega_A5S_by_rep) <- gsub("^percentage_usage_mean.", "", colnames(mega_A5S_by_rep))

mega_A5S_matrix <- apply(mega_A5S, 2, convert_to_numeric)
rownames(mega_A5S_matrix) <- rownames(mega_A5S)

mega_A5S_matrix_by_rep <- apply(mega_A5S_by_rep, 2, convert_to_numeric)
rownames(mega_A5S_matrix_by_rep) <- rownames(mega_A5S_by_rep)

mega_A5S_matrix <- mega_A5S_matrix[rowSums(!is.na(mega_A5S_matrix)) > 0,]

mega_A5S_matrix_by_rep <- mega_A5S_matrix_by_rep[rowSums(!is.na(mega_A5S_matrix_by_rep)) > 0,]

colors <- colorRampPalette(brewer.pal(6, "Reds"))(250)
pheatmap(mega_A5S_matrix,
         na_col = "#FFFFFF",
         border_color = "black",
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = F,
         number_format = "%.0f",
         number_color = "black",
         fontsize = 7,
         fontsize_number = 6.5,
         show_rownames = T,
         show_colnames = T,
         angle_col = 90,
         col = colors,
         main = "Percent Usage of A5S Events by Cell Type")

pheatmap(mega_A5S_matrix_by_rep,
         na_col = "#FFFFFF",
         border_color = "black",
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = F,
         number_format = "%.0f",
         number_color = "black",
         fontsize = 7,
         fontsize_number = 6.5,
         show_rownames = T,
         show_colnames = T,
         angle_col = 90,
         col = colors,
         main = "Percent Usage of A5S Events by Replicate")

write.csv(mega_A5S_matrix, file = "mega_A5S_matrix.csv")
write.csv(mega_A5S_matrix_by_rep, file = "mega_A5S_matrix_by_rep.csv")
```

# Individual cell type-composite results:
```{r}
setwd("E:/Zach Wolfe's JUM analysis/JUM_results")

composite_indiv_cell_types_list <- vector("list", length(singletypes))
names(composite_indiv_cell_types_list) <- singletypes

for (i in singletypes) {
  for (j in singletypes) {
    if (i != j) {
      composite_indiv_cell_types_list[[i]][j] <- vector("list", length(singletypes) - 1)
    }
  }
}

for (i in singletypes) {
  for (j in singletypes){
    if (i != j){
    if (j > i){
composite <- fread(paste0("FINAL_JUM_OUTPUT_pvalue_1", i, "vs", j, "/AS_differential_JUM_output_composite_events_pvalue_1_final_detailed.txt"))
    }
    if (i > j) next

# Remove the deltaPSI column (for now):
composite_deltaPSI_col <- composite$deltaPSI   # Save the column for re-addition later
composite$deltaPSI <- NULL

# Find the columns with names starting with "percentage_usage."
percentage_usage_cols <- grep("^percentage_usage\\.", colnames(composite), value = TRUE)

# List the replicates in this comparison:
rep_names <- sub("^percentage_usage\\.", "", percentage_usage_cols)

# initialize a new data frame which will be the combined results of every composite event in a comparison:
composite_final <- data.frame()

for (group_id in unique(composite$AS_event_ID)) {
  group_df <- filter(composite, AS_event_ID == group_id)
  sub_group <- list()
  for (structure_ID in unique(group_df$AS_structure_ID)){
   sub_group[[structure_ID]] <- group_df[group_df$AS_structure_ID == structure_ID]
  }

  sub_group_length <- length(sub_group)
 
for (q in 1:sub_group_length){
  percent <- sub_group[[q]][, ..percentage_usage_cols]
  percent <- as.data.frame(lapply(percent, function(x) as.numeric(gsub("%", "", x)) / 100))
  sub_group[[q]] <- data.frame(sub_group[[q]]$AS_structure_ID, sub_group[[q]]$sub_junction_ID, percent)
  names(sub_group[[q]])[1:2] <- c("AS_structure_ID", "sub_junction_ID")
  }

end_group <- list()
  
if (sub_group_length == 2){
for (q1 in 1:sub_group_length) {
  for (q2 in 1:sub_group_length) {
    if (q1 != q2) {
      subgroup1 <- sub_group[[q1]]
      subgroup2 <- sub_group[[q2]]
      
      for (col in percentage_usage_cols) {
        
        # Multiply the values and create new columns
        for (r1 in 1:nrow(subgroup1)){
          for (r2 in 1:nrow(subgroup2)){
               multiplied_col <- subgroup1[r1, col] * subgroup2[r2, col]
               subgroup1[r1, paste("multiplied", col, subgroup2$AS_structure_ID[r2], subgroup2$sub_junction_ID[r2], sep = " ")] <- multiplied_col
               subgroup2[r2, paste("multiplied", col, subgroup1$AS_structure_ID[r1],  subgroup1$sub_junction_ID[r1], sep = " ")] <- multiplied_col
          }
        }
      }
      end_group <- append(end_group, list(subgroup1))
    }
  }
}

  counter <- 1
  for (q1 in 1:sub_group_length){
     for (q2 in 1:sub_group_length){
         if (q1 != q2){
         names(end_group)[counter] <- paste(names(sub_group)[q1], "x", names(sub_group)[q2])
         counter <- counter + 1
         }
     }
  }
}

if (sub_group_length == 3){
for (q1 in 1:sub_group_length) {
  for (q2 in 1:sub_group_length) {
    for (q3 in 1:sub_group_length) {
     if (length(unique(c(q1, q2, q3))) == 3) {
      subgroup1 <- sub_group[[q1]]
      subgroup2 <- sub_group[[q2]]
      subgroup3 <- sub_group[[q3]]
      
      for (col in percentage_usage_cols) {
        
        # Multiply the values and create new columns
        for (r1 in 1:nrow(subgroup1)){
          for (r2 in 1:nrow(subgroup2)){
            for (r3 in 1:nrow(subgroup3)){
               multiplied_col <- subgroup1[r1, col] * subgroup2[r2, col] * subgroup3[r3, col]
               subgroup1[r1, paste("multiplied", col, subgroup2$AS_structure_ID[r2],  subgroup2$sub_junction_ID[r2], subgroup3$AS_structure_ID[r3],  subgroup3$sub_junction_ID[r3], sep = " ")] <- multiplied_col
               subgroup2[r2, paste("multiplied", col, subgroup1$AS_structure_ID[r1],  subgroup1$sub_junction_ID[r1], subgroup3$AS_structure_ID[r3],  subgroup3$sub_junction_ID[r3], sep = " ")] <- multiplied_col
               subgroup3[r3, paste("multiplied", col, subgroup1$AS_structure_ID[r1],  subgroup1$sub_junction_ID[r1], subgroup2$AS_structure_ID[r2],  subgroup2$sub_junction_ID[r2], sep = " ")] <- multiplied_col
          }
        }
        }
      }
      end_group <- append(end_group, list(subgroup1))
      }
     }
  }
 }

  counter <- 1
  for (q1 in 1:sub_group_length){
     for (q2 in 1:sub_group_length){
       for (q3 in 1:sub_group_length){
         if (length(unique(c(q1, q2, q3))) == 3){
         names(end_group)[counter] <- paste(names(sub_group)[q1], "x", names(sub_group)[q2], "x", names(sub_group)[q3])
         counter <- counter + 1
             }
           }
         }
         }
}

if (sub_group_length == 4){
for (q1 in 1:sub_group_length) {
  for (q2 in 1:sub_group_length) {
    for (q3 in 1:sub_group_length) {
      for (q4 in 1:sub_group_length){
     if (length(unique(c(q1, q2, q3, q4))) == 4) {
      subgroup1 <- sub_group[[q1]]
      subgroup2 <- sub_group[[q2]]
      subgroup3 <- sub_group[[q3]]
      subgroup4 <- sub_group[[q4]]
      
      for (col in percentage_usage_cols) {
        
        # Multiply the values and create new columns
        for (r1 in 1:nrow(subgroup1)){
          for (r2 in 1:nrow(subgroup2)){
            for (r3 in 1:nrow(subgroup3)){
              for (r4 in 1:nrow(subgroup4)){
               multiplied_col <- subgroup1[r1, col] * subgroup2[r2, col] * subgroup3[r3, col] * subgroup4[r4, col]
               subgroup1[r1, paste("multiplied", col, subgroup2$AS_structure_ID[r2],  subgroup2$sub_junction_ID[r2], subgroup3$AS_structure_ID[r3],  subgroup3$sub_junction_ID[r3], subgroup4$AS_structure_ID[r4],  subgroup4$sub_junction_ID[r4], sep = " ")] <- multiplied_col
               subgroup2[r2, paste("multiplied", col, subgroup1$AS_structure_ID[r1],  subgroup1$sub_junction_ID[r1], subgroup3$AS_structure_ID[r3],  subgroup3$sub_junction_ID[r3], subgroup4$AS_structure_ID[r4],  subgroup4$sub_junction_ID[r4], sep = " ")] <- multiplied_col
               subgroup3[r3, paste("multiplied", col, subgroup1$AS_structure_ID[r1],  subgroup1$sub_junction_ID[r1], subgroup2$AS_structure_ID[r2],  subgroup2$sub_junction_ID[r2], subgroup4$AS_structure_ID[r4],  subgroup4$sub_junction_ID[r4], sep = " ")] <- multiplied_col
               subgroup4[r4, paste("multiplied_", col, subgroup1$AS_structure_ID[r1],  subgroup1$sub_junction_ID[r1], subgroup2$AS_structure_ID[r2],  subgroup2$sub_junction_ID[r2], subgroup3$AS_structure_ID[r3],  subgroup3$sub_junction_ID[r3], sep = " ")] <- multiplied_col
          }
        }
          }
        }
      }
      end_group <- append(end_group, list(subgroup1))
      }
     }
  }
  }
}

  counter <- 1
  for (q1 in 1:sub_group_length){
     for (q2 in 1:sub_group_length){
       for (q3 in 1:sub_group_length){
         for (q4 in 1:sub_group_length){
         if (length(unique(c(q1, q2, q3, q4))) == 4){
         names(end_group)[counter] <- paste(names(sub_group)[q1], "x", names(sub_group)[q2], "x", names(sub_group)[q3], "x", names(sub_group)[q4])
         counter <- counter + 1
         }
             }
           }
         }
         }
}

if (sub_group_length == 5){
for (q1 in 1:sub_group_length) {
  for (q2 in 1:sub_group_length) {
    for (q3 in 1:sub_group_length) {
      for (q4 in 1:sub_group_length) {
        for (q5 in 1:sub_group_length) {
          if (length(unique(c(q1, q2, q3, q4, q5))) == 5) {
            subgroup1 <- sub_group[[q1]]
            subgroup2 <- sub_group[[q2]]
            subgroup3 <- sub_group[[q3]]
            subgroup4 <- sub_group[[q4]]
            subgroup5 <- sub_group[[q5]]
            
            for (col in percentage_usage_cols) {
              for (r1 in 1:nrow(subgroup1)) {
                for (r2 in 1:nrow(subgroup2)) {
                  for (r3 in 1:nrow(subgroup3)) {
                    for (r4 in 1:nrow(subgroup4)) {
                      for (r5 in 1:nrow(subgroup5)) {
                        multiplied_col <- subgroup1[r1, col] * subgroup2[r2, col] * subgroup3[r3, col] * subgroup4[r4, col] * subgroup5[r5, col]

               subgroup1[r1, paste("multiplied", col, subgroup2$AS_structure_ID[r2],  subgroup2$sub_junction_ID[r2], subgroup3$AS_structure_ID[r3],  subgroup3$sub_junction_ID[r3], subgroup4$AS_structure_ID[r4],  subgroup4$sub_junction_ID[r4], subgroup5$AS_structure_ID[r5],  subgroup5$sub_junction_ID[r5], sep = " ")] <- multiplied_col
               subgroup2[r2, paste("multiplied", col, subgroup1$AS_structure_ID[r1],  subgroup1$sub_junction_ID[r1], subgroup3$AS_structure_ID[r3],  subgroup3$sub_junction_ID[r3], subgroup4$AS_structure_ID[r4],  subgroup4$sub_junction_ID[r4], subgroup5$AS_structure_ID[r5],  subgroup5$sub_junction_ID[r5], sep = " ")] <- multiplied_col
               subgroup3[r3, paste("multiplied", col, subgroup1$AS_structure_ID[r1],  subgroup1$sub_junction_ID[r1], subgroup2$AS_structure_ID[r2],  subgroup2$sub_junction_ID[r2], subgroup4$AS_structure_ID[r4],  subgroup4$sub_junction_ID[r4], subgroup5$AS_structure_ID[r5],  subgroup5$sub_junction_ID[r5], sep = " ")] <- multiplied_col
               subgroup4[r4, paste("multiplied", col, subgroup1$AS_structure_ID[r1],  subgroup1$sub_junction_ID[r1], subgroup2$AS_structure_ID[r2],  subgroup2$sub_junction_ID[r2], subgroup3$AS_structure_ID[r3],  subgroup3$sub_junction_ID[r3], subgroup5$AS_structure_ID[r5],  subgroup5$sub_junction_ID[r5], sep = " ")] <- multiplied_col  
               subgroup5[r5, paste("multiplied", col, subgroup1$AS_structure_ID[r1],  subgroup1$sub_junction_ID[r1], subgroup2$AS_structure_ID[r2],  subgroup2$sub_junction_ID[r2], subgroup3$AS_structure_ID[r3],  subgroup3$sub_junction_ID[r3], subgroup4$AS_structure_ID[r4],  subgroup4$sub_junction_ID[r4], sep = " ")] <- multiplied_col        
                      }
                    }
                  }
                }
              }
            }
   end_group <- append(end_group, list(subgroup1))
          }
        }
      }
    }
  }
}
    counter <- 1
  for (q1 in 1:sub_group_length){
     for (q2 in 1:sub_group_length){
       for (q3 in 1:sub_group_length){
         for (q4 in 1:sub_group_length){
           for (q5 in 1:sub_group_length){
         if (length(unique(c(q1, q2, q3, q4, q5))) == 5){
         names(end_group)[counter] <- paste(names(sub_group)[q1], "x", names(sub_group)[q2], "x", names(sub_group)[q3], "x", names(sub_group)[q4], "x", names(sub_group)[q5])
         counter <- counter + 1
         }
             }
           }
         }
     }
  }
}

if (sub_group_length == 6){
for (q1 in seq_len(sub_group_length)) {
  for (q2 in seq_len(sub_group_length)) {
    for (q3 in seq_len(sub_group_length)) {
      for (q4 in seq_len(sub_group_length)) {
        for (q5 in seq_len(sub_group_length)) {
          for (q6 in seq_len(sub_group_length)) {
            if (length(unique(c(q1, q2, q3, q4, q5, q6))) == 6) {
              subgroup1 <- sub_group[[q1]]
              subgroup2 <- sub_group[[q2]]
              subgroup3 <- sub_group[[q3]]
              subgroup4 <- sub_group[[q4]]
              subgroup5 <- sub_group[[q5]]
              subgroup6 <- sub_group[[q6]]

              for (col in percentage_usage_cols) {
                for (r1 in 1:nrow(subgroup1)) {
                  for (r2 in 1:nrow(subgroup2)) {
                    for (r3 in 1:nrow(subgroup3)) {
                      for (r4 in 1:nrow(subgroup4)) {
                        for (r5 in 1:nrow(subgroup5)) {
                          for (r6 in 1:nrow(subgroup6)) {
                            multiplied_col <- subgroup1[r1, col] * subgroup2[r2, col] * subgroup3[r3, col] * subgroup4[r4, col] * subgroup5[r5, col] * subgroup6[r6, col]

subgroup1[r1, paste("multiplied", col, subgroup2$AS_structure_ID[r2], subgroup2$sub_junction_ID[r2], subgroup3$AS_structure_ID[r3], subgroup3$sub_junction_ID[r3], subgroup4$AS_structure_ID[r4], subgroup4$sub_junction_ID[r4], subgroup5$AS_structure_ID[r5], subgroup5$sub_junction_ID[r5], subgroup6$AS_structure_ID[r6], subgroup6$sub_junction_ID[r6], sep = " ")] <- multiplied_col

subgroup2[r2, paste("multiplied_", col, subgroup1$AS_structure_ID[r1], subgroup1$sub_junction_ID[r1], subgroup3$AS_structure_ID[r3], subgroup3$sub_junction_ID[r3], subgroup4$AS_structure_ID[r4], subgroup4$sub_junction_ID[r4], subgroup5$AS_structure_ID[r5], subgroup5$sub_junction_ID[r5], subgroup6$AS_structure_ID[r6], subgroup6$sub_junction_ID[r6], sep = " ")] <- multiplied_col

subgroup3[r3, paste("multiplied_", col, subgroup1$AS_structure_ID[r1], subgroup1$sub_junction_ID[r1], subgroup2$AS_structure_ID[r2], subgroup2$sub_junction_ID[r2], subgroup4$AS_structure_ID[r4], subgroup4$sub_junction_ID[r4], subgroup5$AS_structure_ID[r5], subgroup5$sub_junction_ID[r5], subgroup6$AS_structure_ID[r6], subgroup6$sub_junction_ID[r6], sep = " ")] <- multiplied_col

subgroup4[r4, paste("multiplied_", col, subgroup1$AS_structure_ID[r1], subgroup1$sub_junction_ID[r1], subgroup2$AS_structure_ID[r2], subgroup2$sub_junction_ID[r2], subgroup3$AS_structure_ID[r3], subgroup3$sub_junction_ID[r3], subgroup5$AS_structure_ID[r5], subgroup5$sub_junction_ID[r5], subgroup6$AS_structure_ID[r6], subgroup6$sub_junction_ID[r6], sep = " ")] <- multiplied_col

subgroup5[r5, paste("multiplied_", col, subgroup1$AS_structure_ID[r1], subgroup1$sub_junction_ID[r1], subgroup2$AS_structure_ID[r2], subgroup2$sub_junction_ID[r2], subgroup3$AS_structure_ID[r3], subgroup3$sub_junction_ID[r3], subgroup4$AS_structure_ID[r4], subgroup4$sub_junction_ID[r4], subgroup6$AS_structure_ID[r6], subgroup6$sub_junction_ID[r6], sep = " ")] <- multiplied_col

subgroup6[r6, paste("multiplied_", col, subgroup1$AS_structure_ID[r1], subgroup1$sub_junction_ID[r1], subgroup2$AS_structure_ID[r2], subgroup2$sub_junction_ID[r2], subgroup3$AS_structure_ID[r3], subgroup3$sub_junction_ID[r3], subgroup4$AS_structure_ID[r4], subgroup4$sub_junction_ID[r4], subgroup5$AS_structure_ID[r5], subgroup5$sub_junction_ID[r5], sep = " ")] <- multiplied_col

                          }
                        }
                      }
                    }
                  }
                }
              }
              end_group <- append(end_group, list(subgroup1))
            }
          }
        }
      }
    }
  }
}

counter <- 1
for (q1 in seq_len(sub_group_length)) {
  for (q2 in seq_len(sub_group_length)) {
    for (q3 in seq_len(sub_group_length)) {
      for (q4 in seq_len(sub_group_length)) {
        for (q5 in seq_len(sub_group_length)) {
          for (q6 in seq_len(sub_group_length)) {
            if (length(unique(c(q1, q2, q3, q4, q5, q6))) == 6) {
              names(end_group)[counter] <- paste(names(sub_group)[q1], "x", names(sub_group)[q2], "x", names(sub_group)[q3], "x", names(sub_group)[q4], "x", names(sub_group)[q5], "x", names(sub_group)[q6])
              counter <- counter + 1
            }
          }
        }
      }
    }
  }
}
}

### Running these loops may be computationally difficult. We have decided to only examine composite events with 6 or fewer junctions, as the benefit of examining 7+ sub-junctions does not outweigh the computational cost. If you want to examine composite events with 7 or more junctions, you can use the following loop/template (which is similar to the previous loops)

# if (sub_group_length == 7){
# for (q1 in seq_len(sub_group_length)) {
#   for (q2 in seq_len(sub_group_length)) {
#     for (q3 in seq_len(sub_group_length)) {
#       for (q4 in seq_len(sub_group_length)) {
#         for (q5 in seq_len(sub_group_length)) {
#           for (q6 in seq_len(sub_group_length)) {
#             for (q7 in seq_len(sub_group_length)) {
#               if (length(unique(c(q1, q2, q3, q4, q5, q6, q7))) == 7) {
#                 subgroup1 <- sub_group[[q1]]
#                 subgroup2 <- sub_group[[q2]]
#                 subgroup3 <- sub_group[[q3]]
#                 subgroup4 <- sub_group[[q4]]
#                 subgroup5 <- sub_group[[q5]]
#                 subgroup6 <- sub_group[[q6]]
#                 subgroup7 <- sub_group[[q7]]
# 
#                 for (col in percentage_usage_cols) {
#                   for (r1 in 1:nrow(subgroup1)) {
#                     for (r2 in 1:nrow(subgroup2)) {
#                       for (r3 in 1:nrow(subgroup3)) {
#                         for (r4 in 1:nrow(subgroup4)) {
#                           for (r5 in 1:nrow(subgroup5)) {
#                             for (r6 in 1:nrow(subgroup6)) {
#                               for (r7 in 1:nrow(subgroup7)) {
#                                 multiplied_col <- subgroup1[r1, col] * subgroup2[r2, col] * subgroup3[r3, col] * subgroup4[r4, col] * subgroup5[r5, col] * subgroup6[r6, col] * subgroup7[r7, col]
# 
#                                 subgroup1[r1, paste("multiplied", col, subgroup2$AS_structure_ID[r2], subgroup2$sub_junction_ID[r2], subgroup3$AS_structure_ID[r3], subgroup3$sub_junction_ID[r3], subgroup4$AS_structure_ID[r4], subgroup4$sub_junction_ID[r4], subgroup5$AS_structure_ID[r5], subgroup5$sub_junction_ID[r5], subgroup6$AS_structure_ID[r6], subgroup6$sub_junction_ID[r6], subgroup7$AS_structure_ID[r7], subgroup7$sub_junction_ID[r7], sep = " ")] <- multiplied_col
#                                 subgroup2[r2, paste("multiplied_", col, subgroup1$AS_structure_ID[r1], subgroup1$sub_junction_ID[r1], subgroup3$AS_structure_ID[r3], subgroup3$sub_junction_ID[r3], subgroup4$AS_structure_ID[r4], subgroup4$sub_junction_ID[r4], subgroup5$AS_structure_ID[r5], subgroup5$sub_junction_ID[r5], subgroup6$AS_structure_ID[r6], subgroup6$sub_junction_ID[r6], subgroup7$AS_structure_ID[r7], subgroup7$sub_junction_ID[r7], sep = " ")] <- multiplied_col
#                                 subgroup3[r3, paste("multiplied_", col, subgroup1$AS_structure_ID[r1], subgroup1$sub_junction_ID[r1], subgroup2$AS_structure_ID[r2], subgroup2$sub_junction_ID[r2], subgroup4$AS_structure_ID[r4], subgroup4$sub_junction_ID[r4], subgroup5$AS_structure_ID[r5], subgroup5$sub_junction_ID[r5], subgroup6$AS_structure_ID[r6], subgroup6$sub_junction_ID[r6], subgroup7$AS_structure_ID[r7], subgroup7$sub_junction_ID[r7], sep = " ")] <- multiplied_col
#                                 subgroup4[r4, paste("multiplied_", col, subgroup1$AS_structure_ID[r1], subgroup1$sub_junction_ID[r1], subgroup2$AS_structure_ID[r2], subgroup2$sub_junction_ID[r2], subgroup3$AS_structure_ID[r3], subgroup3$sub_junction_ID[r3], subgroup5$AS_structure_ID[r5], subgroup5$sub_junction_ID[r5], subgroup6$AS_structure_ID[r6], subgroup6$sub_junction_ID[r6], subgroup7$AS_structure_ID[r7], subgroup7$sub_junction_ID[r7], sep = " ")] <- multiplied_col
#                                 subgroup5[r5, paste("multiplied_", col, subgroup1$AS_structure_ID[r1], subgroup1$sub_junction_ID[r1], subgroup2$AS_structure_ID[r2], subgroup2$sub_junction_ID[r2], subgroup3$AS_structure_ID[r3], subgroup3$sub_junction_ID[r3], subgroup4$AS_structure_ID[r4], subgroup4$sub_junction_ID[r4], subgroup6$AS_structure_ID[r6], subgroup6$sub_junction_ID[r6], subgroup7$AS_structure_ID[r7], subgroup7$sub_junction_ID[r7], sep = " ")] <- multiplied_col
#                                 subgroup6[r6, paste("multiplied_", col, subgroup1$AS_structure_ID[r1], subgroup1$sub_junction_ID[r1], subgroup2$AS_structure_ID[r2], subgroup2$sub_junction_ID[r2], subgroup3$AS_structure_ID[r3], subgroup3$sub_junction_ID[r3], subgroup4$AS_structure_ID[r4], subgroup4$sub_junction_ID[r4], subgroup5$AS_structure_ID[r5], subgroup5$sub_junction_ID[r5], subgroup7$AS_structure_ID[r7], subgroup7$sub_junction_ID[r7], sep = " ")] <- multiplied_col
#                                 subgroup7[r7, paste("multiplied_", col, subgroup1$AS_structure_ID[r1], subgroup1$sub_junction_ID[r1], subgroup2$AS_structure_ID[r2], subgroup2$sub_junction_ID[r2], subgroup3$AS_structure_ID[r3], subgroup3$sub_junction_ID[r3], subgroup4$AS_structure_ID[r4], subgroup4$sub_junction_ID[r4], subgroup5$AS_structure_ID[r5], subgroup5$sub_junction_ID[r5], subgroup6$AS_structure_ID[r6], subgroup6$sub_junction_ID[r6], sep = " ")] <- multiplied_col
#                               }
#                             }
#                           }
#                         }
#                       }
#                     }
#                   }
#                 }
#                 end_group <- append(end_group, list(subgroup1))
#               }
#             }
#           }
#         }
#       }
#     }
#   }
# }
# 
# counter <- 1
# for (q1 in seq_len(sub_group_length)) {
#   for (q2 in seq_len(sub_group_length)) {
#     for (q3 in seq_len(sub_group_length)) {
#       for (q4 in seq_len(sub_group_length)) {
#         for (q5 in seq_len(sub_group_length)) {
#           for (q6 in seq_len(sub_group_length)) {
#             for (q7 in seq_len(sub_group_length)) {
#               if (length(unique(c(q1, q2, q3, q4, q5, q6, q7))) == 7) {
#                 names(end_group)[counter] <- paste(names(sub_group)[q1], "x", names(sub_group)[q2], "x", names(sub_group)[q3], "x", names(sub_group)[q4], "x", names(sub_group)[q5], "x", names(sub_group)[q6], "x", names(sub_group)[q7])
#                 counter <- counter + 1
#               }
#             }
#           }
#         }
#       }
#     }
#   }
# }
# }


## For now, we have to remove all AS_events with more than 6 AS_structure_IDs because of compute time. The :
if (sub_group_length > 6){
 composite <- composite[!composite$AS_event_ID == group_id,]
}

# bind_rows:
if (sub_group_length <= 6){
  comp_with_juncs <- data.frame()
  for (k in 1:length(end_group)){
     comp_with_juncs <- bind_rows(end_group[[k]], comp_with_juncs)      
  }
  
  rownames(comp_with_juncs) <- paste(comp_with_juncs$AS_structure_ID, comp_with_juncs$sub_junction_ID, rownames(comp_with_juncs))
 
  comp_with_juncs <- comp_with_juncs %>%
  mutate(row_names = rownames(comp_with_juncs))
   
  comp_with_juncs$row_names <- sub(" \\d+$", "", comp_with_juncs$row_names)

  # Find rows with duplicate "row_names"
  duplicate_row_names <- comp_with_juncs$row_names[duplicated(comp_with_juncs$row_names)]

# Loop through each unique duplicate row name
for (name in unique(duplicate_row_names)) {
  # Find rows with the same duplicate row name
  matching_rows <- which(comp_with_juncs$row_names == name)
  
  # Check if there are exactly two matching rows
  if (length(matching_rows) == 2) {
    # Iterate through columns
    for (col_idx in (3+length(rep_names)):(ncol(comp_with_juncs) - 1)) {    # Exclude the original replicate names and row_name columns
      value1 <- comp_with_juncs[matching_rows[1], col_idx]
      value2 <- comp_with_juncs[matching_rows[2], col_idx]
      
      # Replace NA values in the first row with values from the second row
      if (is.na(value1) && !is.na(value2)) {
        comp_with_juncs[matching_rows[1], col_idx] <- value2
      } else if (is.na(value2) && !is.na(value1)) {
        comp_with_juncs[matching_rows[2], col_idx] <- value1
      }
    }
  }
  comp_with_juncs <- comp_with_juncs[!duplicated(comp_with_juncs$row_names),]
}
  
  rownames(comp_with_juncs) <- sub(" \\d+$", "", rownames(comp_with_juncs))

  # Get rid of the row_names column we created earlier:
  comp_with_juncs <- comp_with_juncs[, -ncol(comp_with_juncs)]

  # Some values of multiplied percentage_usage for some replicates is 'NaN'. This will be a problem for us; I am going to temporarily resolve this issue by replacing all 'NaN's with 1000 - a number much higher than any percentage value could account for in this pipeline. I will replace the 1000's with NaN's again later in the pipeline:
  is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
  comp_with_juncs[is.nan(comp_with_juncs)] <- 1000
  
  # Remove columns which contain all NA values but keep columns with NaN values
  comp_with_juncs <- comp_with_juncs[, colSums(is.na(comp_with_juncs)) != nrow(comp_with_juncs)]
  
  # Initialize an empty list to store the new columns
  new_columns <- list()

# Loop through each column
for (col_name in colnames(comp_with_juncs)[(3+length(rep_names)):(ncol(comp_with_juncs))]) {
  # Get the non-NA values in the column
  non_na_values <- comp_with_juncs[!is.na(comp_with_juncs[, col_name]), col_name]
  
  # Create new columns based on the number of non-NA values
  for (n in 1:length(non_na_values)) {
    new_col_name <- paste(col_name, rownames(comp_with_juncs)[n], sep = " ")
    new_col <- rep(NA, nrow(comp_with_juncs))
    new_col[which(!is.na(comp_with_juncs[, col_name]))[n]] <- non_na_values[n]
    new_columns[[new_col_name]] <- new_col
  }
}

# Remove NA values:
  for (n in 1:length(new_columns)){
    new_columns[[n]] <- new_columns[[n]][!is.na(new_columns[[n]])]
  }
  
  # Remove empty elements from the list
  new_columns <- new_columns[sapply(new_columns, function(x) length(x) > 0)]

  # Combine the new columns into a dataframe:
  comp_with_juncs_final <- t(data.frame(new_columns))
  # Replace "." with space in row names
  rownames(comp_with_juncs_final) <- gsub("\\.", " ", rownames(comp_with_juncs_final))
  
  # Count how many rows are in every sub_group, then multiply the row counts by each other:
  total_rows <- prod(sapply(sub_group, function(x) nrow(x)))
  # Multiply the total_rows by the number of replicates. This is how many valid values our comp_with_juncs_final should have:
  rows_to_keep <- as.numeric(total_rows * length(rep_names))
  comp_with_juncs_final <- as.data.frame(comp_with_juncs_final[1:rows_to_keep,])

  comp_with_juncs_final[c(rep_names)] <- NA
  
  # Copy data from V1 to the new column for matching rows:
  for (r in 1:length(rep_names)) {
  col_name <- rep_names[r]
  matching_rows <- grepl(col_name, rownames(comp_with_juncs_final))
  comp_with_juncs_final[matching_rows, col_name] <- comp_with_juncs_final[matching_rows, 1]
  }
  
  rownames(comp_with_juncs_final)[1:total_rows] <- sub(paste0("^multiplied percentage_usage ", rep_names[1], " "), "", rownames(comp_with_juncs_final)[1:total_rows])
 
  for (r in rep_names){
   col_values <- comp_with_juncs_final[[r]]

   # Find the non-NA values:
   non_na_values <- col_values[!is.na(col_values)]

   # Replace the first xxx rows with non-NA values:
   comp_with_juncs_final[1:total_rows, r] <- non_na_values
  }
  
  # Lastly, remove the remaining rows which are not completely populated, as well as the original "V1" column:
  comp_with_juncs_final <- comp_with_juncs_final[1:total_rows,-1]
  
  # We still need to convert values = 1000 back to 'NaN':
  comp_with_juncs_final[comp_with_juncs_final == 1000] <- NaN
  
  rownames(comp_with_juncs_final) <- paste(group_df$Gene[1], group_df$AS_event_ID[1], rownames(comp_with_juncs_final))
  
  # combine comp_with_juncs_final with our previously initialized data frame:
  composite_final <- rbind(composite_final, comp_with_juncs_final)
 }
}

# And finally we can add columns that represent the means of the percentage_usage_means for each cell type:
i_cols <- grep(i, colnames(composite_final), value = TRUE)
i_vals <- subset(composite_final, select = i_cols)
i_vals[] <- lapply(i_vals, replace_non_finite_with_na)
i_vals$mean <- as.numeric(0)
j_cols <- grep(j, colnames(composite_final), value = TRUE)
j_vals <- subset(composite_final, select = j_cols)
j_vals[] <- lapply(j_vals, replace_non_finite_with_na)
j_vals$mean <- as.numeric(0)
for (k in 1:nrow(i_vals)){
i_vals[k,"mean"] <- mean(as.numeric(i_vals[k, -which(names(i_vals) == "mean")]), na.rm = T)
j_vals[k,"mean"] <- mean(as.numeric(j_vals[k, -which(names(j_vals) == "mean")]), na.rm = T)
}

# Define the column names for i and j
col_name_i <- paste0("average_", i)
col_name_j <- paste0("average_", j)

  # Modify mean columns
  composite_final[, (col_name_i)] <- i_vals$mean
  composite_final[, (col_name_j)] <- j_vals$mean

  # Multuiply all values by 100 before adding "%":
  composite_final <- composite_final * 100
  
for (p in i_cols){
  composite_final[[p]] <- paste0(composite_final[[p]], "%")
}

for (q in j_cols){
  composite_final[[q]] <- paste0(composite_final[[q]], "%")
}

average_cols <- grep("^average_", colnames(composite_final), value = TRUE)

for (p in average_cols){
  composite_final[[p]] <- paste0(composite_final[[p]], "%")
}

composite_final[] <- lapply(composite_final, replace_na_percent)
composite_final[] <- lapply(composite_final, replace_inf_percent)

# Save as a .csv:
write.csv(composite_final, file = paste0(i, "_", j, "_detailed_composite_R.csv"), row.names = T, na = "")
write.csv(composite_final, file = paste0(j, "_", i, "_detailed_composite_R.csv"), row.names = T, na = "")
    }
    if (i == j) next
  }
}
```

# Now we need to combine all composite events from each cell type so that we get an event ID, gene name, and average % spliced in for each cell type:
```{r}
setwd("E:/Zach Wolfe's JUM analysis/JUM_results")

composite_indiv_usage_combined <- vector("list", length(singletypes))
names(composite_indiv_usage_combined) <- singletypes

composite_indiv_usage_by_rep <- vector("list", length(singletypes_withreps))
names(composite_indiv_usage_by_rep) <- singletypes_withreps

# rbind all files for cell type XXX:
for (i in singletypes) {
  composite_indiv_usage_combined[[i]] <- data.table()
  composite_indiv_usage_by_rep[[i]] <- data.table()

  matches <- grep(paste0("^", i), singletypes_withreps, value = TRUE)
  
  for (j in singletypes) {
    if (i != j) {
      composite <- fread(paste0(i, "_", j, "_detailed_composite_R.csv"), header = T)
      split_first_col <- strsplit(as.character(composite$V1), " ", fixed = TRUE)
      composite$Gene <- sapply(split_first_col, `[`, 1)
      composite$AS_event_ID <- sapply(split_first_col, function(x) paste(x[-1], collapse = " "))

      composite <- composite[!is.na(composite[[paste0("average_", i)]]),]
      composite_by_rep <- composite[, c("Gene", "AS_event_ID", matches), with = FALSE]
      composite <- composite[, c("Gene", "AS_event_ID", paste0("average_", i)), with = FALSE]
      
      # Use rbindlist to combine data tables with different numbers of columns
      composite_indiv_usage_combined[[i]] <- rbindlist(list(composite_indiv_usage_combined[[i]], composite), fill = TRUE)
      composite_indiv_usage_combined[[i]] <- unique(composite_indiv_usage_combined[[i]])
      
      composite_indiv_usage_by_rep[[i]] <- rbindlist(list(composite_indiv_usage_by_rep[[i]], composite_by_rep), fill = TRUE)
      composite_indiv_usage_by_rep[[i]] <- unique(composite_indiv_usage_by_rep[[i]])
    }
    if (i == j) next
  }

  composite_indiv_usage_combined[[i]] <- composite_indiv_usage_combined[[i]][rowSums(is.na(composite_indiv_usage_combined[[i]])) != ncol(composite_indiv_usage_combined[[i]]) - 2, ]
  composite_indiv_usage_by_rep[[i]] <- composite_indiv_usage_by_rep[[i]][rowSums(is.na(composite_indiv_usage_by_rep[[i]])) != ncol(composite_indiv_usage_by_rep[[i]]) - 2, ]
  
  write.csv(composite_indiv_usage_combined[[i]], file = paste0("combined_composite_usage_", i, ".csv"), row.names = F)
  write.csv(composite_indiv_usage_by_rep[[i]], file = paste0("combined_composite_usage_by_rep_", i, ".csv"), row.names = F)
}

# Lastly, remove the sublists which are labeled by individual replicate (while keeping the sublists representing whole cell types)
composite_indiv_usage_by_rep <- Filter(function(x) !is.null(x), composite_indiv_usage_by_rep)
```

# optional - eliminate all composite events that are not the max % usage for each A.S. event ID:
```{r}
setwd("E:/Zach Wolfe's JUM analysis/JUM_results")

composite_indiv_usage_combined_max <- vector("list", length(singletypes))
names(composite_indiv_usage_combined_max) <- singletypes

for (i in singletypes){
composite_indiv_usage_combined_max[[i]] <- fread(paste0("combined_composite_usage_", i, ".csv"))

setDT(composite_indiv_usage_combined_max[[i]])

# Extract the AS_event_ID without J00x information:
composite_indiv_usage_combined_max[[i]][, AS_event_ID_base := gsub("\\s\\d_Junction.*", "", AS_event_ID)]

# Convert % usage column to numeric:
composite_indiv_usage_combined_max[[i]][, paste0("average_", i)] <- convert_to_numeric(composite_indiv_usage_combined_max[[i]][, paste0("average_", i)])

# Order the data by % usage in descending order:
composite_indiv_usage_combined_max[[i]] <- composite_indiv_usage_combined_max[[i]][order(-paste0("average_", i))]

# Keep only the row with the max % usage for each unique AS_event_ID_base:
composite_indiv_usage_combined_max[[i]] <- composite_indiv_usage_combined_max[[i]][, .SD[which.max(paste0("average_", i))], by = AS_event_ID_base]

# Remove the intermediate column AS_event_ID_base:
composite_indiv_usage_combined_max[[i]][, AS_event_ID_base := NULL]

write.csv(composite_indiv_usage_combined_max[[i]], file = paste0("composite_events_", i, "_max_%_usage_per_row.csv"))
}
```

# Make a heatmap for all A.S. events represented by % usage for each cell type:
```{r}
mega_composite <- data.frame(Gene = character(), AS_event_ID = character())
mega_composite_by_rep <- data.frame(Gene = character(), AS_event_ID = character())

for (i in rev(singletypes)) {
composite <- composite_indiv_usage_combined[[i]]
composite <- as.data.frame(composite)
composite[,3] <- as.numeric(sub("%", "", composite[,3])) / 100
composite_mean <- aggregate(. ~ Gene + AS_event_ID, data = composite, FUN = function(x) mean(as.numeric(x)), na.action = NULL)
composite_mean[,3] <- paste0(composite_mean[,3] * 100, "%")
mega_composite <- merge(composite_mean, mega_composite, by = c("Gene", "AS_event_ID"), all = TRUE)
rownames(mega_composite) <- paste(mega_composite$Gene, mega_composite$AS_event_ID) # You can label both the gene name and the AS_event_ID as rownames, or edit this line to label rownames as just one or the other
}

for (i in singletypes) {
composite <- composite_indiv_usage_by_rep[[i]]
composite <- as.data.frame(composite)
 numeric_cols <- 3:ncol(composite)
 composite[, numeric_cols] <- lapply(composite[, numeric_cols], function(x) {
   ifelse(grepl("%", x, fixed = TRUE), sub("%", "", x), x)
 })
composite_mean <- aggregate(. ~ Gene + AS_event_ID, data = composite, FUN = function(x) mean(as.numeric(x)), na.action = NULL)
mega_composite_by_rep <- full_join(mega_composite_by_rep, composite_mean, by = c("Gene", "AS_event_ID"))
mega_composite_by_rep <- mega_composite_by_rep %>% arrange(tolower(Gene))
rownames(mega_composite_by_rep) <- paste(mega_composite_by_rep$Gene, mega_composite_by_rep$AS_event_ID) # You can label both the gene name and the AS_event_ID as rownames, or edit this line to label rownames as just one or the other
}

mega_composite <- mega_composite[,c(3:ncol(mega_composite))]
colnames(mega_composite) <- gsub("^average_", "", colnames(mega_composite))

mega_composite_by_rep <- mega_composite_by_rep[,c(3:ncol(mega_composite_by_rep))]
colnames(mega_composite_by_rep) <- gsub("^percentage_usage_mean.", "", colnames(mega_composite_by_rep))

mega_composite_matrix <- apply(mega_composite, 2, convert_to_numeric)
rownames(mega_composite_matrix) <- rownames(mega_composite)

mega_composite_matrix_by_rep <- apply(mega_composite_by_rep, 2, convert_to_numeric)
rownames(mega_composite_matrix_by_rep) <- rownames(mega_composite_by_rep)

mega_composite_matrix <- mega_composite_matrix[rowSums(!is.na(mega_composite_matrix)) > 0,]

mega_composite_matrix_by_rep <- mega_composite_matrix_by_rep[rowSums(!is.na(mega_composite_matrix_by_rep)) > 0,]

# organize the composite matrices by junctions with >10% usage:
mega_composite_matrix_over_10 <- mega_composite_matrix[which(apply(mega_composite_matrix, 1, max, na.rm = TRUE) > 0.1),]
mega_composite_matrix_by_rep_over_10 <- mega_composite_matrix_by_rep[which(apply(mega_composite_matrix_by_rep, 1, max, na.rm = TRUE) > 0.1),]

# order by a cell type (in this case, using AFD as an example):
mega_composite_AFD <- mega_composite_matrix[order(-mega_composite_matrix[,"AFD"]),]

colors <- colorRampPalette(brewer.pal(6, "BuPu"))(250)
pheatmap(mega_composite_AFD[1:30,],
         na_col = "#FFFFFF",
         border_color = "black",
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = F,
         number_format = "%.0f",
         number_color = "black",
         fontsize = 7,
         fontsize_number = 6.5,
         show_rownames = F,
         show_colnames = T,
         angle_col = 90,
         col = colors,
         main = "Percent Usage of Composite Events by Cell Type (ordered by AFD)")

pheatmap(mega_composite_matrix_by_rep,
         na_col = "#FFFFFF",
         border_color = "black",
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = F,
         number_format = "%.0f",
         number_color = "black",
         fontsize = 7,
         fontsize_number = 6.5,
         show_rownames = T,
         show_colnames = T,
         angle_col = 90,
         col = colors,
         main = "Percent Usage of Composite Events by Replicate")

write.csv(mega_composite_matrix, file = "mega_composite_matrix.csv")
write.csv(mega_composite_matrix_by_rep, file = "mega_composite_matrix_by_rep.csv")
```

# Composite heatmap: Choose one gene and print all AS events/junctions for that gene
```{r}
# Create a new column for gene name:
mega_composite_matrix$Gene <- sub(" .*", "", rownames(mega_composite_matrix))

# Subset the data for rows of a given gene name:
subset_matrix <- mega_composite_matrix[mega_composite_matrix$Gene == "unc-75",1:46]
#subset_matrix <- subset_matrix[order(-subset_matrix[, "AVM"]), ]   # optional - order by a cell type of your choice

concatenate_non_na <- function(x) paste(na.omit(x), collapse = "_")

for (i in 1:nrow(subset_matrix)){
elements <- unlist(strsplit(rownames(subset_matrix)[i], " ", fixed = TRUE))
print(length(elements))
if (length(elements) == 6){
subset_matrix[i,"Gene"] <- elements[1]
subset_matrix[i,"AS_event_ID"] <- elements[2]
subset_matrix[i,"Junc_number_1"] <- elements[4]
subset_matrix[i,"Junc_number_2"] <- elements[6]
}
if (length(elements) == 8){
subset_matrix[i,"Gene"] <- elements[1]
subset_matrix[i,"AS_event_ID"] <- elements[2]
subset_matrix[i,"Junc_number_1"] <- elements[4]
subset_matrix[i,"Junc_number_2"] <- elements[6]
subset_matrix[i,"Junc_number_3"] <- elements[8]
}
if (length(elements) == 10){
subset_matrix[i,"Gene"] <- elements[1]
subset_matrix[i,"AS_event_ID"] <- elements[2]
subset_matrix[i,"Junc_number_1"] <- elements[4]
subset_matrix[i,"Junc_number_2"] <- elements[6]
subset_matrix[i,"Junc_number_3"] <- elements[8]
subset_matrix[i,"Junc_number_4"] <- elements[10]
}
if (length(elements) == 12){
subset_matrix[i,"Gene"] <- elements[1]
subset_matrix[i,"AS_event_ID"] <- elements[2]
subset_matrix[i,"Junc_number_1"] <- elements[4]
subset_matrix[i,"Junc_number_2"] <- elements[6]
subset_matrix[i,"Junc_number_3"] <- elements[8]
subset_matrix[i,"Junc_number_4"] <- elements[10]
subset_matrix[i,"Junc_number_5"] <- elements[12]
}
if (length(elements) == 14){
subset_matrix[i,"Gene"] <- elements[1]
subset_matrix[i,"AS_event_ID"] <- elements[2]
subset_matrix[i,"Junc_number_1"] <- elements[4]
subset_matrix[i,"Junc_number_2"] <- elements[6]
subset_matrix[i,"Junc_number_3"] <- elements[8]
subset_matrix[i,"Junc_number_4"] <- elements[10]
subset_matrix[i,"Junc_number_5"] <- elements[12]
subset_matrix[i,"Junc_number_6"] <- elements[14]
}
}

# Loop over each column and apply the concatenate function
for (col in singletypes) {
  subset_matrix[[paste0("collapsed_", col)]] <- apply(subset_matrix[, col, drop = FALSE], 1, concatenate_non_na)
}

### STOP: What is the max # of junctions for the gene name you chose? This will affect the following steps - read carefully! ###

# Create a new data frame with collapsed values:
collapsed_matrix <- subset_matrix[c("AS_event_ID", "Junc_number_1", "Junc_number_2", "Junc_number_3", "Junc_number_4", "Junc_number_5", "Junc_number_6", paste0("collapsed_", singletypes), "Gene")]   ## add on Junc_number_4, Junc_number_5, Junc_number_6 if applicable
collapsed_matrix$Junc_number_3[is.na(collapsed_matrix$Junc_number_3)] <- "JXXX"   ## repeat this step for each number of unused junctions in your dataframe
collapsed_matrix$Junc_number_4[is.na(collapsed_matrix$Junc_number_4)] <- "JXXX"
collapsed_matrix$Junc_number_5[is.na(collapsed_matrix$Junc_number_5)] <- "JXXX"
collapsed_matrix$Junc_number_6[is.na(collapsed_matrix$Junc_number_6)] <- "JXXX"

# Use aggregate to combine collapsed values
collapsed_matrix <- aggregate(. ~ AS_event_ID + Junc_number_1 + Junc_number_2 + Junc_number_3 + Junc_number_4 + Junc_number_5 + Junc_number_6 + Gene, data = collapsed_matrix, FUN = function(x) paste(unique(na.omit(x)), collapse = " "))   ## add on Junc_number_4, Junc_number_5, Junc_number_6 if applicable
collapsed_matrix[, 9:ncol(collapsed_matrix)] <- apply(collapsed_matrix[, 9:ncol(collapsed_matrix)], 2, as.numeric)   ## change the column start based on your max # of junctions
collapsed_matrix <- collapsed_matrix[, c(8, 1:7, 9:ncol(collapsed_matrix))]
names(collapsed_matrix) <- sub("^collapsed_", "", names(collapsed_matrix))

for (i in 1:nrow(collapsed_matrix)) {
  if (collapsed_matrix$Junc_number_6[i] != "JXXX") {   # here, you will have to repeat the if else loop for every unused junction
    rownames(collapsed_matrix)[i] <- paste(collapsed_matrix$AS_event_ID[i], collapsed_matrix$Junc_number_1[i], collapsed_matrix$Junc_number_2[i], collapsed_matrix$Junc_number_3[i], collapsed_matrix$Junc_number_4[i], collapsed_matrix$Junc_number_5[i], collapsed_matrix$Junc_number_6[i])
  } 
  else if (collapsed_matrix$Junc_number_5[i] != "JXXX" && collapsed_matrix$Junc_number_6[i] == "JXXX") {   
    rownames(collapsed_matrix)[i] <- paste(collapsed_matrix$AS_event_ID[i], collapsed_matrix$Junc_number_1[i], collapsed_matrix$Junc_number_2[i], collapsed_matrix$Junc_number_3[i], collapsed_matrix$Junc_number_4[i], collapsed_matrix$Junc_number_5[i])
  }
  else if (collapsed_matrix$Junc_number_4[i] != "JXXX" && collapsed_matrix$Junc_number_5[i] == "JXXX") {  
    rownames(collapsed_matrix)[i] <- paste(collapsed_matrix$AS_event_ID[i], collapsed_matrix$Junc_number_1[i], collapsed_matrix$Junc_number_2[i], collapsed_matrix$Junc_number_3[i], collapsed_matrix$Junc_number_4[i])
  }
  else if (collapsed_matrix$Junc_number_3[i] != "JXXX" && collapsed_matrix$Junc_number_4[i] == "JXXX") {  
    rownames(collapsed_matrix)[i] <- paste(collapsed_matrix$AS_event_ID[i], collapsed_matrix$Junc_number_1[i], collapsed_matrix$Junc_number_2[i], collapsed_matrix$Junc_number_3[i])
  } 
  else {
    rownames(collapsed_matrix)[i] <- paste(collapsed_matrix$AS_event_ID[i], collapsed_matrix$Junc_number_1[i], collapsed_matrix$Junc_number_2[i])
  }
}

colors <- colorRampPalette(brewer.pal(6, "BuPu"))(250)

# Split the data based on unique values in AS_event_ID
unique_as_event_ids <- unique(collapsed_matrix$AS_event_ID)
heatmaps <- list()

for (as_event_id in unique_as_event_ids) {
  subset_data <- subset(collapsed_matrix, AS_event_ID == as_event_id)
  subset_heatmap <- collapsed_matrix[collapsed_matrix$AS_event_ID == as_event_id, ]
  write.csv(subset_heatmap, paste0("heatmap ", subset_matrix$Gene[1], " ", gsub(" ", "_", as_event_id), ".csv"), row.names = TRUE)

  
## If your gene does not have a different number of junctions in a given AS_event_ID, you can run this simpler heatmap function:
# heatmap <- pheatmap(
#      subset_data[, 9:ncol(subset_data)],
#      na_col = "#FFFFFF",
#      border_color = "black",
#      cluster_rows = F,
#      cluster_cols = F,
#      display_numbers = F,
#      number_format = "%.0f",
#      number_color = "black",
#      fontsize = 7,
#      show_rownames = T,
#      show_colnames = T,
#      angle_col = 90,
#      col = colors,
#      main = paste("Percent Usage of unc-75 Junctions by Cell Type")
#    )
#    heatmaps[[as_event_id]] <- heatmap

  
## ... otherwise you will have to separate each AS_event_ID heatmap into multiple heatmaps based on the number of junctions:
# Declare "JXXX" values that are in columns 3 through X:
junc_columns <- subset_data[, 3:8]

# Count the number of "JXXX" occurrences in each row
num_JXXX <- rowSums(junc_columns == "JXXX")

# Create a list to store heatmaps
heatmaps_list <- vector("list", length = max(num_JXXX) + 1)

# Split data based on the number of "JXXX" occurrences
split_data <- split(subset_data[, 9:ncol(subset_data)], num_JXXX)

# Create a heatmap for each subset
for (i in seq_along(split_data)) {
  heatmap <- pheatmap(
    split_data[[i]],
    na_col = "#FFFFFF",
    border_color = "black",
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    display_numbers = FALSE,
    number_format = "%.0f",
    number_color = "black",
    fontsize = 7,
    show_rownames = TRUE,
    show_colnames = TRUE,
    angle_col = 90,
    col = colors,
    main = paste("Percent Usage of unc-75 Junctions by Cell Type:", as_event_id, "Subset", i)
  )
  heatmaps_list[[i]] <- heatmap
}
}
```

# Sometimes RefSeq and WormBase gene names do not always align. Due to this inconsistency, we should read in an updated and correct gtf file (from WormBase):
```{r}
gtf <- fread("caenorhabditis_elegans.PRJNA13758.WBPS18.canonical_geneset.gtf")
colnames(gtf) <- c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
# replace the combined attribute column with multiple columns that are more appropriate for naming:
gtf$gene_id <- sub('.*gene_id "(.*?)".*', '\\1', gtf$attribute)
gtf$gene_version <- sub('.*gene_version "(.*?)".*', '\\1', gtf$attribute)
gtf$gene_biotype <- sub('.*gene_biotype "(.*?)".*', '\\1', gtf$attribute)
gtf$gene_name <- sub('.*gene_name "(.*?)".*', '\\1', gtf$attribute)

# Remove the original 'attribute' column:
gtf <- gtf[, c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "gene_id", "gene_version", "gene_biotype", "gene_name")]
```

# Extract the coordinates of each AS event. This will be used to produce the correct and updated gene names for AS events with more than one gene name:
```{r}
## A3S events:
for (i in which(grepl(";", mega_A3S$Gene))) { 
  as_event_id <- mega_A3S$AS_event_ID[i]

  # Extract chromosome, strand, start, and end coordinates from the AS event ID:
  as_parts <- unlist(strsplit(as_event_id, "_"))
  chr_var <- as_parts[1]
  strand_var <- as_parts[2]
  start_var <- as.numeric(as_parts[3])
  end_var <- as.numeric(as_parts[length(as_parts)])

  # Find the corresponding gene in the GTF dataframe:
  matching_gene <- subset(gtf, chr == chr_var & feature == "gene" & start <= start_var & end >= end_var)
  
  if (nrow(matching_gene) == 0){
  matching_gene_start <- subset(gtf, chr == chr_var & feature == "gene" & start <= start_var)
  matching_gene_end <- subset(gtf, chr == chr_var & feature == "gene" & end >= end_var)
  matching_gene_start <- matching_gene_start[nrow(matching_gene_start),]
  matching_gene_end <- matching_gene_end[1,]
  ## consider if matching_gene_start$gene_name == matching_gene_end$gene_name, then paste gene_name
  if (matching_gene_start$end >= start_var & matching_gene_end$start <= end_var){
  gene_name <- paste0(matching_gene_start$gene_name, "; ", matching_gene_end$gene_name)
  WB_gene_name <- paste0(matching_gene_start$gene_id, "; ", matching_gene_end$gene_id)
  }
  else if (matching_gene_start$end >= start_var & matching_gene_end$start >= end_var){
  gene_name <- matching_gene_start$gene_name
  WB_gene_name <- matching_gene_start$gene_id
  }
  else if (matching_gene_start$end <= start_var & matching_gene_end$start <= end_var){
  gene_name <- matching_gene_end$gene_name
  WB_gene_name <- matching_gene_end$gene_id
  }
  else {
  gene_name <- NA
  WB_gene_name <- NA
  }

  # Replace the original gene name with the updated one:
  mega_A3S$Gene[i] <- gene_name
  mega_A3S$WormBase_Gene_ID[i] <- WB_gene_name
  }
  
  if (nrow(matching_gene) > 0){
  gene_name <- matching_gene$gene_name
  WB_gene_name <- matching_gene$gene_id

  # Replace the original gene name with the updated one:
  mega_A3S$Gene[i] <- gene_name
  mega_A3S$WormBase_Gene_ID[i] <- WB_gene_name
  }
}

# Don't forget to change the sequence name to match the correct and updated version!
for (g in 1:nrow(mega_A3S)) {
  gene_ids <- unlist(strsplit(as.character(mega_A3S$WormBase_Gene_ID[g]), "; "))
  matching_sequence_names <- character(length(gene_ids))
  
  for (i in seq_along(gene_ids)) {
    matching_row <- all_genes_conversion[all_genes_conversion$`WormBase Gene ID` == gene_ids[i], ]
    if (nrow(matching_row) > 0) {
      matching_sequence_names[i] <- matching_row$`Sequence Name`
    }
  }
  mega_A3S$Sequence_Name[g] <- paste(matching_sequence_names, collapse = "; ")
  rownames(mega_A3S)[g] <- paste(mega_A3S$Gene[g], mega_A3S$AS_event_ID[g])
}

# Save the new, corrected matrix as a .csv:
write.csv(mega_A3S, file = "mega_A3S_matrix.csv")

## A5S events:
for (i in which(grepl(";", mega_A5S$Gene))) {
  as_event_id <- mega_A5S$AS_event_ID[i]

  # Extract chromosome, strand, start, and end coordinates from the AS event ID:
  as_parts <- unlist(strsplit(as_event_id, "_"))
  chr_var <- as_parts[1]
  strand_var <- as_parts[2]
  start_var <- as.numeric(as_parts[3])
  end_var <- as.numeric(as_parts[length(as_parts)])

  # Find the corresponding gene in the GTF dataframe:
  matching_gene <- subset(gtf, chr == chr_var & feature == "gene" & start <= start_var & end >= end_var)
  
  if (nrow(matching_gene) == 0){
  matching_gene_start <- subset(gtf, chr == chr_var & feature == "gene" & start <= start_var)
  matching_gene_end <- subset(gtf, chr == chr_var & feature == "gene" & end >= end_var)
  matching_gene_start <- matching_gene_start[nrow(matching_gene_start),]
  matching_gene_end <- matching_gene_end[1,]
  if (matching_gene_start$end >= start_var & matching_gene_end$start <= end_var){
  gene_name <- paste0(matching_gene_start$gene_name, "; ", matching_gene_end$gene_name)
  WB_gene_name <- paste0(matching_gene_start$gene_id, "; ", matching_gene_end$gene_id)
  }
  else if (matching_gene_start$end >= start_var & matching_gene_end$start >= end_var){
  gene_name <- matching_gene_start$gene_name
  WB_gene_name <- matching_gene_start$gene_id
  }
  else if (matching_gene_start$end <= start_var & matching_gene_end$start <= end_var){
  gene_name <- matching_gene_end$gene_name
  WB_gene_name <- matching_gene_end$gene_id
  }
  else {
  gene_name <- NA
  WB_gene_name <- NA
  }

  # Replace the original gene name with the updated one:
  mega_A5S$Gene[i] <- gene_name
  mega_A5S$WormBase_Gene_ID[i] <- WB_gene_name
  }
  
  if (nrow(matching_gene) > 0){
  gene_name <- matching_gene$gene_name
  WB_gene_name <- matching_gene$gene_id

  # Replace the original gene name with the updated one:
  mega_A5S$Gene[i] <- gene_name
  mega_A5S$WormBase_Gene_ID[i] <- WB_gene_name
  }
}

# Don't forget to change the sequence name to match the correct and updated version!
for (g in 1:nrow(mega_A5S)) {
  gene_ids <- unlist(strsplit(as.character(mega_A5S$WormBase_Gene_ID[g]), "; "))
  matching_sequence_names <- character(length(gene_ids))
  
  for (i in seq_along(gene_ids)) {
    matching_row <- all_genes_conversion[all_genes_conversion$`WormBase Gene ID` == gene_ids[i], ]
    if (nrow(matching_row) > 0) {
      matching_sequence_names[i] <- matching_row$`Sequence Name`
    }
  }
  mega_A5S$Sequence_Name[g] <- paste(matching_sequence_names, collapse = "; ")
  rownames(mega_A5S)[g] <- paste(mega_A5S$Gene[g], mega_A5S$AS_event_ID[g])
}

# Save the new, corrected matrix as a .csv:
write.csv(mega_A5S, file = "mega_A5S_matrix.csv")

## Cassette events:
for (i in which(grepl(";", mega_cassette$Gene))) {
  as_event_id <- mega_cassette$AS_event_ID[i]

  # Extract chromosome, strand, start, and end coordinates from the AS event ID:
  as_parts <- unlist(strsplit(as_event_id, "_"))
  chr_var <- as_parts[1]
  strand_var <- as_parts[2]
  start_var <- as.numeric(as_parts[3])
  end_var <- as.numeric(as_parts[length(as_parts)])

  # Find the corresponding gene in the GTF dataframe:
  matching_gene <- subset(gtf, chr == chr_var & feature == "gene" & start <= start_var & end >= end_var)
  
  if (nrow(matching_gene) == 0){
  matching_gene_start <- subset(gtf, chr == chr_var & feature == "gene" & start <= start_var)
  matching_gene_end <- subset(gtf, chr == chr_var & feature == "gene" & end >= end_var)
  matching_gene_start <- matching_gene_start[nrow(matching_gene_start),]
  matching_gene_end <- matching_gene_end[1,]
  if (matching_gene_start$end >= start_var & matching_gene_end$start <= end_var){
  gene_name <- paste0(matching_gene_start$gene_name, "; ", matching_gene_end$gene_name)
  WB_gene_name <- paste0(matching_gene_start$gene_id, "; ", matching_gene_end$gene_id)
  }
  else if (matching_gene_start$end >= start_var & matching_gene_end$start >= end_var){
  gene_name <- matching_gene_start$gene_name
  WB_gene_name <- matching_gene_start$gene_id
  }
  else if (matching_gene_start$end <= start_var & matching_gene_end$start <= end_var){
  gene_name <- matching_gene_end$gene_name
  WB_gene_name <- matching_gene_end$gene_id
  }
  else {
  gene_name <- NA
  WB_gene_name <- NA
  }

  # Replace the original gene name with the updated one:
  mega_cassette$Gene[i] <- gene_name
  mega_cassette$WormBase_Gene_ID[i] <- WB_gene_name
  }
  
  if (nrow(matching_gene) > 0){
  gene_name <- matching_gene$gene_name
  WB_gene_name <- matching_gene$gene_id

  # Replace the original gene name with the updated one:
  mega_cassette$Gene[i] <- gene_name
  mega_cassette$WormBase_Gene_ID[i] <- WB_gene_name
  }
}

# Don't forget to change the sequence name to match the correct and updated version!
for (g in 1:nrow(mega_cassette)) {
  gene_ids <- unlist(strsplit(as.character(mega_cassette$WormBase_Gene_ID[g]), "; "))
  matching_sequence_names <- character(length(gene_ids))
  
  for (i in seq_along(gene_ids)) {
    matching_row <- all_genes_conversion[all_genes_conversion$`WormBase Gene ID` == gene_ids[i], ]
    if (nrow(matching_row) > 0) {
      matching_sequence_names[i] <- matching_row$`Sequence Name`
    }
  }
  mega_cassette$Sequence_Name[g] <- paste(matching_sequence_names, collapse = "; ")
  rownames(mega_cassette)[g] <- paste(mega_cassette$Gene[g], mega_cassette$AS_event_ID[g])
}

# Save the new, corrected matrix as a .csv:
write.csv(mega_cassette, file = "mega_cassette_matrix.csv")

# ## Composite events:
# for (i in which(grepl(";", mega_composite$Gene))) {
#   as_event_id <- mega_composite$AS_event_ID[i]
# 
#   # Extract chromosome, strand, start, and end coordinates from the AS event ID:
#   as_parts <- unlist(strsplit(as_event_id, "_"))
#   chr_var <- as_parts[1]
#   strand_var <- as_parts[2]
#   start_var <- as.numeric(as_parts[3])
#   end_var <- as.numeric(as_parts[length(as_parts)])
# 
#   # Find the corresponding gene in the GTF dataframe:
#   matching_gene <- subset(gtf, chr == chr_var & feature == "gene" & start <= start_var & end >= end_var)
# 
#   # Extract gene name:
#   gene_name <- matching_gene$gene_name
#   WB_gene_name <- matching_gene$gene_id
#   
#   # Replace the original gene name with the updated one:
#   mega_composite$Gene[i] <- gene_name
#   mega_composite$WormBase_Gene_ID[i] <- WB_gene_name
# }
# 
# # Don't forget to change the sequence name to match the correct and updated version!
#  for (g in 1:nrow(mega_composite)) {
#     matching_row <- all_genes_conversion[all_genes_conversion$Gene == mega_composite$Gene[g], ]
#       mega_composite$Sequence_Name[g] <- matching_row$`Sequence Name`
#       rownames(mega_composite)[g] <- paste(mega_composite$Gene[g], mega_composite$AS_event_ID[g])
#  }
# 
# # Save the new, corrected matrix as a .csv:
# write.csv(mega_composite, file = "mega_composite_matrix.csv")

## Intron retention events:
for (i in which(grepl(";", mega_intron$Gene))) {
  as_event_id <- mega_intron$AS_event_ID[i]

  # Extract chromosome, strand, start, and end coordinates from the AS event ID:
  as_parts <- unlist(strsplit(as_event_id, "_"))
  chr_var <- as_parts[1]
  strand_var <- as_parts[2]
  start_var <- as.numeric(as_parts[3])
  end_var <- as.numeric(as_parts[length(as_parts)])

  # Find the corresponding gene in the GTF dataframe:
  matching_gene <- subset(gtf, chr == chr_var & feature == "gene" & start <= start_var & end >= end_var)
  
  if (nrow(matching_gene) == 0){
  matching_gene_start <- subset(gtf, chr == chr_var & feature == "gene" & start <= start_var)
  matching_gene_end <- subset(gtf, chr == chr_var & feature == "gene" & end >= end_var)
  matching_gene_start <- matching_gene_start[nrow(matching_gene_start),]
  matching_gene_end <- matching_gene_end[1,]
  if (matching_gene_start$end >= start_var & matching_gene_end$start <= end_var){
  gene_name <- paste0(matching_gene_start$gene_name, "; ", matching_gene_end$gene_name)
  WB_gene_name <- paste0(matching_gene_start$gene_id, "; ", matching_gene_end$gene_id)
  }
  else if (matching_gene_start$end >= start_var & matching_gene_end$start >= end_var){
  gene_name <- matching_gene_start$gene_name
  WB_gene_name <- matching_gene_start$gene_id
  }
  else if (matching_gene_start$end <= start_var & matching_gene_end$start <= end_var){
  gene_name <- matching_gene_end$gene_name
  WB_gene_name <- matching_gene_end$gene_id
  }
  else {
  gene_name <- NA
  WB_gene_name <- NA
  }

  # Replace the original gene name with the updated one:
  mega_intron$Gene[i] <- gene_name
  mega_intron$WormBase_Gene_ID[i] <- WB_gene_name
  }
  
  if (nrow(matching_gene) > 0){
  gene_name <- matching_gene$gene_name
  WB_gene_name <- matching_gene$gene_id

  # Replace the original gene name with the updated one:
  mega_intron$Gene[i] <- gene_name
  mega_intron$WormBase_Gene_ID[i] <- WB_gene_name
  }
}

# Don't forget to change the sequence name to match the correct and updated version!
for (g in 1:nrow(mega_intron)) {
  gene_ids <- unlist(strsplit(as.character(mega_intron$WormBase_Gene_ID[g]), "; "))
  matching_sequence_names <- character(length(gene_ids))
  
  for (i in seq_along(gene_ids)) {
    matching_row <- all_genes_conversion[all_genes_conversion$`WormBase Gene ID` == gene_ids[i], ]
    if (nrow(matching_row) > 0) {
      matching_sequence_names[i] <- matching_row$`Sequence Name`
    }
  }
  mega_intron$Sequence_Name[g] <- paste(matching_sequence_names, collapse = "; ")
  rownames(mega_intron)[g] <- paste(mega_intron$Gene[g], mega_intron$AS_event_ID[g])
}

# Save the new, corrected matrix as a .csv:
write.csv(mega_intron, file = "mega_intron_matrix.csv")

## MXE events:
for (i in which(grepl(";", mega_MXE$Gene))) {
  as_event_id <- mega_MXE$AS_event_ID[i]

  # Extract chromosome, strand, start, and end coordinates from the AS event ID:
  as_parts <- unlist(strsplit(as_event_id, "_"))
  chr_var <- as_parts[1]
  strand_var <- as_parts[2]
  start_var <- as.numeric(as_parts[3])
  end_var <- as.numeric(as_parts[length(as_parts)])

  # Find the corresponding gene in the GTF dataframe:
  matching_gene <- subset(gtf, chr == chr_var & feature == "gene" & start <= start_var & end >= end_var)
  
  if (nrow(matching_gene) == 0){
  matching_gene_start <- subset(gtf, chr == chr_var & feature == "gene" & start <= start_var)
  matching_gene_end <- subset(gtf, chr == chr_var & feature == "gene" & end >= end_var)
  matching_gene_start <- matching_gene_start[nrow(matching_gene_start),]
  matching_gene_end <- matching_gene_end[1,]
  if (matching_gene_start$end >= start_var & matching_gene_end$start <= end_var){
  gene_name <- paste0(matching_gene_start$gene_name, "; ", matching_gene_end$gene_name)
  WB_gene_name <- paste0(matching_gene_start$gene_id, "; ", matching_gene_end$gene_id)
  }
  else if (matching_gene_start$end >= start_var & matching_gene_end$start >= end_var){
  gene_name <- matching_gene_start$gene_name
  WB_gene_name <- matching_gene_start$gene_id
  }
  else if (matching_gene_start$end <= start_var & matching_gene_end$start <= end_var){
  gene_name <- matching_gene_end$gene_name
  WB_gene_name <- matching_gene_end$gene_id
  }
  else {
  gene_name <- NA
  WB_gene_name <- NA
  }

  # Replace the original gene name with the updated one:
  mega_MXE$Gene[i] <- gene_name
  mega_MXE$WormBase_Gene_ID[i] <- WB_gene_name
  }
  
  if (nrow(matching_gene) > 0){
  gene_name <- matching_gene$gene_name
  WB_gene_name <- matching_gene$gene_id

  # Replace the original gene name with the updated one:
  mega_MXE$Gene[i] <- gene_name
  mega_MXE$WormBase_Gene_ID[i] <- WB_gene_name
  }
}

# Don't forget to change the sequence name to match the correct and updated version!
for (g in 1:nrow(mega_MXE)) {
  gene_ids <- unlist(strsplit(as.character(mega_MXE$WormBase_Gene_ID[g]), "; "))
  matching_sequence_names <- character(length(gene_ids))
  
  for (i in seq_along(gene_ids)) {
    matching_row <- all_genes_conversion[all_genes_conversion$`WormBase Gene ID` == gene_ids[i], ]
    if (nrow(matching_row) > 0) {
      matching_sequence_names[i] <- matching_row$`Sequence Name`
    }
  }
  mega_MXE$Sequence_Name[g] <- paste(matching_sequence_names, collapse = "; ")
  rownames(mega_MXE)[g] <- paste(mega_MXE$Gene[g], mega_MXE$AS_event_ID[g])
}

# Save the new, corrected matrix as a .csv:
write.csv(mega_MXE, file = "mega_MXE_matrix.csv")
```

# PCAs of different splice types:
```{r}
condition$replicate <- rownames(condition)
### A3S events:
mega_A3S_PCA <- sapply(mega_A3S_matrix_by_rep[,5:ncol(mega_A3S_matrix_by_rep)], function(x) as.numeric(gsub("%", "", x)))
mega_A3S_PCA <- mega_A3S_PCA[rowSums(is.na(mega_A3S_PCA)) != ncol(mega_A3S_PCA), ]

for(i in 1:ncol(mega_A3S_PCA)) {
  mega_A3S_PCA[is.na(mega_A3S_PCA[,i]), i] <- mean(mega_A3S_PCA[,i], na.rm = TRUE)
}

data_t <- t(mega_A3S_PCA)
data_scaled <- scale(data_t)
pca_res <- prcomp(data_scaled, center = TRUE, scale. = TRUE)
pca_var <- pca_res$sdev^2
pca_var_perc <- pca_var / sum(pca_var) * 100
pca_data <- data.frame(pca_res$x)
pca_data$CellType <- rownames(pca_data)

pca_data <- pca_data %>%
  mutate(CellTypeBase = ifelse(grepl("\\d$", CellType), CellType, sub(".*(\\d+)$", "\\1", CellType)))

pca_data <- left_join(pca_data, condition, by = c("CellTypeBase" = "replicate"))

# plot the PCA:
ggplot(pca_data, aes(x = PC1, y = PC2, label = CellType, color = condition)) +
  geom_point(size = 2) +
  geom_text(aes(label = CellType), hjust = 1.25, vjust = 1.25, size = 3) +
  labs(
    title = "PCA of A3S Events (upstreamness)",
    x = paste0("Principal Component 1 (", round(pca_var_perc[1], 2), "%)"),
    y = paste0("Principal Component 2 (", round(pca_var_perc[2], 2), "%)")
  ) +
  theme_minimal() +
  scale_color_manual(values = c("Sensory" = "dodgerblue", "Motor" = "gold3", "Interneuron" = "salmon", "Polymodal" = "darkgreen")) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

### A5S events:
mega_A5S_PCA <- sapply(mega_A5S_matrix_by_rep[,5:ncol(mega_A5S_matrix_by_rep)], function(x) as.numeric(gsub("%", "", x)))
mega_A5S_PCA <- mega_A5S_PCA[rowSums(is.na(mega_A5S_PCA)) != ncol(mega_A5S_PCA), ]

for(i in 1:ncol(mega_A5S_PCA)) {
  mega_A5S_PCA[is.na(mega_A5S_PCA[,i]), i] <- mean(mega_A5S_PCA[,i], na.rm = TRUE)
}

data_t <- t(mega_A5S_PCA)
data_scaled <- scale(data_t)
pca_res <- prcomp(data_scaled, center = TRUE, scale. = TRUE)
pca_var <- pca_res$sdev^2
pca_var_perc <- pca_var / sum(pca_var) * 100
pca_data <- data.frame(pca_res$x)
pca_data$CellType <- rownames(pca_data)

pca_data <- pca_data %>%
  mutate(CellTypeBase = ifelse(grepl("\\d$", CellType), CellType, sub(".*(\\d+)$", "\\1", CellType)))

pca_data <- left_join(pca_data, condition, by = c("CellTypeBase" = "replicate"))

# plot the PCA:
ggplot(pca_data, aes(x = PC1, y = PC2, label = CellType, color = condition)) +
  geom_point(size = 2) +
  geom_text(aes(label = CellType), hjust = 1.25, vjust = 1.25, size = 3) +
  labs(
    title = "PCA of A5S Events (upstreamness)",
    x = paste0("Principal Component 1 (", round(pca_var_perc[1], 2), "%)"),
    y = paste0("Principal Component 2 (", round(pca_var_perc[2], 2), "%)")
  ) +
  theme_minimal() +
  scale_color_manual(values = c("Sensory" = "dodgerblue", "Motor" = "gold3", "Interneuron" = "salmon", "Polymodal" = "darkgreen")) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

### cassette events:
mega_cassette_PCA <- sapply(mega_cassette_matrix_by_rep[,5:ncol(mega_cassette_matrix_by_rep)], function(x) as.numeric(gsub("%", "", x)))
mega_cassette_PCA <- mega_cassette_PCA[rowSums(is.na(mega_cassette_PCA)) != ncol(mega_cassette_PCA), ]

for(i in 1:ncol(mega_cassette_PCA)) {
  mega_cassette_PCA[is.na(mega_cassette_PCA[,i]), i] <- mean(mega_cassette_PCA[,i], na.rm = TRUE)
}

data_t <- t(mega_cassette_PCA)
data_scaled <- scale(data_t)
pca_res <- prcomp(data_scaled, center = TRUE, scale. = TRUE)
pca_var <- pca_res$sdev^2
pca_var_perc <- pca_var / sum(pca_var) * 100
pca_data <- data.frame(pca_res$x)
pca_data$CellType <- rownames(pca_data)

pca_data <- pca_data %>%
  mutate(CellTypeBase = ifelse(grepl("\\d$", CellType), CellType, sub(".*(\\d+)$", "\\1", CellType)))

pca_data <- left_join(pca_data, condition, by = c("CellTypeBase" = "replicate"))

# plot the PCA:
ggplot(pca_data, aes(x = PC1, y = PC2, label = CellType, color = condition)) +
  geom_point(size = 2) +
  geom_text(aes(label = CellType), hjust = 1.25, vjust = 1.25, size = 3) +
  labs(
    title = "PCA of Cassette Events",
    x = paste0("Principal Component 1 (", round(pca_var_perc[1], 2), "%)"),
    y = paste0("Principal Component 2 (", round(pca_var_perc[2], 2), "%)")
  ) +
  theme_minimal() +
  scale_color_manual(values = c("Sensory" = "dodgerblue", "Motor" = "gold3", "Interneuron" = "salmon", "Polymodal" = "darkgreen")) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

### composite events:
mega_composite_PCA <- sapply(mega_composite_matrix_by_rep[,5:ncol(mega_composite_matrix_by_rep)], function(x) as.numeric(gsub("%", "", x)))
mega_composite_PCA <- mega_composite_PCA[rowSums(is.na(mega_composite_PCA)) != ncol(mega_composite_PCA), ]

for(i in 1:ncol(mega_composite_PCA)) {
  mega_composite_PCA[is.na(mega_composite_PCA[,i]), i] <- mean(mega_composite_PCA[,i], na.rm = TRUE)
}

data_t <- t(mega_composite_PCA)
data_scaled <- scale(data_t)
pca_res <- prcomp(data_scaled, center = TRUE, scale. = TRUE)
pca_var <- pca_res$sdev^2
pca_var_perc <- pca_var / sum(pca_var) * 100
pca_data <- data.frame(pca_res$x)
pca_data$CellType <- rownames(pca_data)

pca_data <- pca_data %>%
  mutate(CellTypeBase = ifelse(grepl("\\d$", CellType), CellType, sub(".*(\\d+)$", "\\1", CellType)))

pca_data <- left_join(pca_data, condition, by = c("CellTypeBase" = "replicate"))

# plot the PCA:
ggplot(pca_data, aes(x = PC1, y = PC2, label = CellType, color = condition)) +
  geom_point(size = 2) +
  geom_text(aes(label = CellType), hjust = 1.25, vjust = 1.25, size = 3) +
  labs(
    title = "PCA of Composite Events",
    x = paste0("Principal Component 1 (", round(pca_var_perc[1], 2), "%)"),
    y = paste0("Principal Component 2 (", round(pca_var_perc[2], 2), "%)")
  ) +
  theme_minimal() +
  scale_color_manual(values = c("Sensory" = "dodgerblue", "Motor" = "gold3", "Interneuron" = "salmon", "Polymodal" = "darkgreen")) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

### intron events:
mega_intron_PCA <- sapply(mega_intron_matrix_by_rep[,5:ncol(mega_intron_matrix_by_rep)], function(x) as.numeric(gsub("%", "", x)))
mega_intron_PCA <- mega_intron_PCA[rowSums(is.na(mega_intron_PCA)) != ncol(mega_intron_PCA), ]

for(i in 1:ncol(mega_intron_PCA)) {
  mega_intron_PCA[is.na(mega_intron_PCA[,i]), i] <- mean(mega_intron_PCA[,i], na.rm = TRUE)
}

data_t <- t(mega_intron_PCA)
data_scaled <- scale(data_t)
pca_res <- prcomp(data_scaled, center = TRUE, scale. = TRUE)
pca_var <- pca_res$sdev^2
pca_var_perc <- pca_var / sum(pca_var) * 100
pca_data <- data.frame(pca_res$x)
pca_data$CellType <- rownames(pca_data)

pca_data <- pca_data %>%
  mutate(CellTypeBase = ifelse(grepl("\\d$", CellType), CellType, sub(".*(\\d+)$", "\\1", CellType)))

pca_data <- left_join(pca_data, condition, by = c("CellTypeBase" = "replicate"))

# plot the PCA:
ggplot(pca_data, aes(x = PC1, y = PC2, label = CellType, color = condition)) +
  geom_point(size = 2) +
  geom_text(aes(label = CellType), hjust = 1.25, vjust = 1.25, size = 3) +
  labs(
    title = "PCA of Intron Retention Events",
    x = paste0("Principal Component 1 (", round(pca_var_perc[1], 2), "%)"),
    y = paste0("Principal Component 2 (", round(pca_var_perc[2], 2), "%)")
  ) +
  theme_minimal() +
  scale_color_manual(values = c("Sensory" = "dodgerblue", "Motor" = "gold3", "Interneuron" = "salmon", "Polymodal" = "darkgreen")) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

### MXE events:
mega_MXE_PCA <- sapply(mega_MXE_matrix_by_rep[,5:ncol(mega_MXE_matrix_by_rep)], function(x) as.numeric(gsub("%", "", x)))
mega_MXE_PCA <- mega_MXE_PCA[rowSums(is.na(mega_MXE_PCA)) != ncol(mega_MXE_PCA), ]

for(i in 1:ncol(mega_MXE_PCA)) {
  mega_MXE_PCA[is.na(mega_MXE_PCA[,i]), i] <- mean(mega_MXE_PCA[,i], na.rm = TRUE)
}

data_t <- t(mega_MXE_PCA)
data_scaled <- scale(data_t)
pca_res <- prcomp(data_scaled, center = TRUE, scale. = TRUE)
pca_var <- pca_res$sdev^2
pca_var_perc <- pca_var / sum(pca_var) * 100
pca_data <- data.frame(pca_res$x)
pca_data$CellType <- rownames(pca_data)

pca_data <- pca_data %>%
  mutate(CellTypeBase = ifelse(grepl("\\d$", CellType), CellType, sub(".*(\\d+)$", "\\1", CellType)))

pca_data <- left_join(pca_data, condition, by = c("CellTypeBase" = "replicate"))

# plot the PCA:
ggplot(pca_data, aes(x = PC1, y = PC2, label = CellType, color = condition)) +
  geom_point(size = 2) +
  geom_text(aes(label = CellType), hjust = 1.25, vjust = 1.25, size = 3) +
  labs(
    title = "PCA of MXE Events (upstreamness)",
    x = paste0("Principal Component 1 (", round(pca_var_perc[1], 2), "%)"),
    y = paste0("Principal Component 2 (", round(pca_var_perc[2], 2), "%)")
  ) +
  theme_minimal() +
  scale_color_manual(values = c("Sensory" = "dodgerblue", "Motor" = "gold3", "Interneuron" = "salmon", "Polymodal" = "darkgreen")) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

### All AS events combined:
all_datasets <- list(mega_A3S_matrix_by_rep, mega_A5S_matrix_by_rep, mega_cassette_matrix_by_rep, mega_composite_matrix_by_rep, mega_intron_matrix_by_rep, mega_MXE_matrix_by_rep)
combined_data <- do.call(rbind, all_datasets)

# Convert necessary columns to numeric, replacing '%' and handling NA values
mega_ALL_PCA <- sapply(combined_data[, 5:ncol(combined_data)], function(x) as.numeric(gsub("%", "", x)))

mega_ALL_PCA <- mega_ALL_PCA[rowSums(is.na(mega_ALL_PCA)) != ncol(mega_ALL_PCA), ]

for(i in 1:ncol(mega_ALL_PCA)) {
  mega_ALL_PCA[is.na(mega_ALL_PCA[,i]), i] <- mean(mega_ALL_PCA[,i], na.rm = TRUE)
}

data_t <- t(mega_ALL_PCA)
data_scaled <- scale(data_t)
pca_res <- prcomp(data_scaled, center = TRUE, scale. = TRUE)
pca_var <- pca_res$sdev^2
pca_var_perc <- pca_var / sum(pca_var) * 100
pca_data <- data.frame(pca_res$x)
pca_data$CellType <- rownames(pca_data)

pca_data <- pca_data %>%
  mutate(CellTypeBase = ifelse(grepl("\\d$", CellType), CellType, sub(".*(\\d+)$", "\\1", CellType)))

pca_data <- left_join(pca_data, condition, by = c("CellTypeBase" = "replicate"))

# plot the PCA:
ggplot(pca_data, aes(x = PC1, y = PC2, label = CellType, color = condition)) +
  geom_point(size = 2) +
  geom_text(aes(label = CellType), hjust = 1.25, vjust = 1.25, size = 3) +
  labs(
    title = "PCA of ALL AS Events Combined",
    x = paste0("Principal Component 1 (", round(pca_var_perc[1], 2), "%)"),
    y = paste0("Principal Component 2 (", round(pca_var_perc[2], 2), "%)")
  ) +
  theme_minimal() +
  scale_color_manual(values = c("Sensory" = "dodgerblue", "Motor" = "gold3", "Interneuron" = "salmon", "Polymodal" = "darkgreen")) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

### All AS events combined (minus composite events):
all_datasets <- list(mega_A3S_matrix_by_rep, mega_A5S_matrix_by_rep, mega_cassette_matrix_by_rep, mega_intron_matrix_by_rep, mega_MXE_matrix_by_rep)
combined_data <- do.call(rbind, all_datasets)

# Convert necessary columns to numeric, replacing '%' and handling NA values
mega_ALL_PCA <- sapply(combined_data[, 5:ncol(combined_data)], function(x) as.numeric(gsub("%", "", x)))

mega_ALL_PCA <- mega_ALL_PCA[rowSums(is.na(mega_ALL_PCA)) != ncol(mega_ALL_PCA), ]

for(i in 1:ncol(mega_ALL_PCA)) {
  mega_ALL_PCA[is.na(mega_ALL_PCA[,i]), i] <- mean(mega_ALL_PCA[,i], na.rm = TRUE)
}

data_t <- t(mega_ALL_PCA)
data_scaled <- scale(data_t)
pca_res <- prcomp(data_scaled, center = TRUE, scale. = TRUE)
pca_var <- pca_res$sdev^2
pca_var_perc <- pca_var / sum(pca_var) * 100
pca_data <- data.frame(pca_res$x)
pca_data$CellType <- rownames(pca_data)

pca_data <- pca_data %>%
  mutate(CellTypeBase = ifelse(grepl("\\d$", CellType), CellType, sub(".*(\\d+)$", "\\1", CellType)))

pca_data <- left_join(pca_data, condition, by = c("CellTypeBase" = "replicate"))

# plot the PCA:
ggplot(pca_data, aes(x = PC1, y = PC2, label = CellType, color = condition)) +
  geom_point(size = 2) +
  geom_text(aes(label = CellType), hjust = 1.25, vjust = 1.25, size = 3) +
  labs(
    title = "PCA of ALL AS Events Combined (minus composite events)",
    x = paste0("Principal Component 1 (", round(pca_var_perc[1], 2), "%)"),
    y = paste0("Principal Component 2 (", round(pca_var_perc[2], 2), "%)")
  ) +
  theme_minimal() +
  scale_color_manual(values = c("Sensory" = "dodgerblue", "Motor" = "gold3", "Interneuron" = "salmon", "Polymodal" = "darkgreen")) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))
```

# unc-31 plots:
```{r}
data_frame_unc31 <- data.frame(`Cell Type` = names(mega_cassette_matrix[1058,]), `Percent Usage` = mega_cassette_matrix[1058,])

# Barplots:
ggplot(data_frame_unc31, aes(x = Cell.Type, y = Percent.Usage, fill = factor(Percent.Usage))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = rainbow(46)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Percent Usage of unc-31 for All Cell Types",
       x = "Cell Type",
       y = "Percent Usage of unc-31") +
  guides(fill = FALSE)  

ggplot(data_frame_unc31, aes(x = reorder(Cell.Type, -Percent.Usage), y = Percent.Usage, fill = factor(Percent.Usage))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = rainbow(46)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Percent Usage of unc-31 for All Cell Types",
       x = "Cell Type",
       y = "Percent Usage of unc-31") +
  guides(fill = FALSE)

data_frame_unc31 <- data.frame(`Cell Type` = names(mega_cassette_matrix_by_rep[1058,]), `Percent Usage` = mega_cassette_matrix_by_rep[1058,])

data_frame_unc31$`Cell Type` <- sub("\\d+", "", data_frame_unc31$Cell.Type)

data_frame_unc31$`Cell Type`[100:110] <- c("I5", "I5", "I5", "I5","IL1", "IL1", "IL1", "IL2","IL2", "IL2", "IL2")

colnames(data_frame_unc31)[1] <- c("Replicate")
data_frame_unc31 <- na.omit(data_frame_unc31)

# Boxplots:
ggplot(data_frame_unc31, aes(x = `Cell Type`, y = Percent.Usage)) +
  geom_boxplot(fill = "lightblue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Percent Usage for unc-31",
       x = "Cell Type",
       y = "Percent Usage") +
  guides(fill = FALSE)

ggplot(data_frame_unc31, aes(x = `Cell Type`, y = Percent.Usage)) +
  stat_boxplot(fill = NA, color = "black", outlier.colour = "black", outlier.shape = NA) +
  geom_jitter(position = position_jitter(width = 0.3), color = "dodgerblue", alpha = 1.2) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Percent Usage for unc-31 Cassette Event",
       x = "Cell Type",
       y = "Percent Usage") +
  guides(fill = FALSE)

ggplot(data_frame_unc31, aes(x = `Cell Type`, y = Percent.Usage, color = `Cell Type`)) +
  stat_boxplot(fill = NA, color = "black", outlier.colour = "black", outlier.shape = NA) +
  geom_jitter(position = position_jitter(width = 0.3), alpha = 1.2) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Percent Usage for unc-31 Cassette Event",
       x = "Cell Type",
       y = "Percent Usage") +
  guides(fill = FALSE)

# Strip plots:
# ggplot(data_frame_unc31, aes(x = `Cell Type`, y = Percent.Usage, color = `Cell Type`)) +
#   geom_jitter(position = position_jitter(width = 0.1), alpha = 1.7) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7.5)) +
#   labs(title = "Strip Plot for unc-31 Cassette Event",
#        x = "Cell Type",
#        y = "Percent Usage") +
#   guides(fill = FALSE)

# create a dataframe of the median percentage usage values of each cell type
summary_stats <- data_frame_unc31 %>%
  group_by(`Cell Type`) %>%
  summarise(median = median(Percent.Usage),
            SEM = sd(Percent.Usage) / sqrt(n()))

# get the order of cell types based on median values:
order_levels <- summary_stats$`Cell Type`[order(-summary_stats$median)]

# reorder the cell types:
data_frame_unc31$`Cell Type` <- factor(data_frame_unc31$`Cell Type`, levels = order_levels)

# Plot:
p <- ggplot(data_frame_unc31, aes(x = `Cell Type`, y = Percent.Usage, color = `Cell Type`)) +
  geom_jitter(position = position_jitter(width = 0.1), alpha = 1, size = 0.55) +
  geom_crossbar(data = summary_stats, aes(x = `Cell Type`, y = median, ymin = median - SEM, ymax = median + SEM), width = 0.5, color = "black") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7.5)) +
  labs(title = "Percent Inclusion of unc-31 Cassette Event",
       x = "Cell Type",
       y = "Percent Inclusion") +
  guides(fill = FALSE) +
  theme(legend.position = "none")

print(p)
ggsave("unc-31_plot.svg", plot = p)
```

Copyright 2024 The Regents of the University of California

All Rights Reserved

Created by Zachery Wolfe

Department of Biochemistry

This file is part of Differential Expression in C. elegans. \
Differential Expression in C. elegans is a free pipeline: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. \
Differential Expression in C. elegans is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. \
You should have received a copy of the GNU General Public License along with Differential Expression in C. elegans. If not, see <https://www.gnu.org/licenses/>.
