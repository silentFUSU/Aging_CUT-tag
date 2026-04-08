check_ages <- function(df) {
  # Define the valid ages for "Old" and "Young"
  valid_old_ages <- c("24m", "25m", "26m", "27m")
  valid_young_ages <- c("2m", "3m")
  
  # Find indices of incorrect ages for "Old"
  wrong_old <- which(grepl("Old", df$sample) & !df$age %in% valid_old_ages)
  # Find indices of incorrect ages for "Young"
  wrong_young <- which(grepl("Young", df$sample) & !df$age %in% valid_young_ages)
  
  # Combine the indices
  wrong_indices <- c(wrong_old, wrong_young)
  
  if (length(wrong_indices) == 0) {
    message("age right")
  } else {
    for (i in wrong_indices) {
      message(sprintf("Row %d: age wrong (%s)", i, df$age[i]))
    }
  }
}
check_sample_uniqueness <- function(df) {
  # Check for non-unique samples
  non_unique_samples <- which(duplicated(df$sample))
  
  # Sample uniqueness check
  if (length(non_unique_samples) == 0) {
    message("sample right")
  } else {
    for (i in non_unique_samples) {
      message(sprintf("Row %d: sample not unique (%s)", i, df$sample[i]))
    }
  }
}
check_exact_tissue_in_sample <- function(df) {
  # Define the mappings between tissue names and their corresponding sample names
  tissue_to_sample_map <- list(
    "Frontal Cortex" = "Frontal_Cortex",
    "Bone Marrow" = "Bone_Marrow",
    "Mammary Gland" = "Mammary_Gland"
  )
  
  # Initialize a flag to check if all tissues have corresponding samples
  all_matched <- TRUE
  
  # Iterate over each row to check if the mapped tissue exists in the sample
  for (i in seq_len(nrow(df))) {
    tissue <- df$organ[i]
    
    # Determine the expected sample substring using mappings or transformations
    expected_sample <- ifelse(tissue %in% names(tissue_to_sample_map), 
                              tissue_to_sample_map[[tissue]], 
                              gsub(" ", "_", tissue))
    
    # Use grep to check for direct match of the expected sample name in the sample column
    if (!grepl(paste0("\\b", expected_sample, "\\b"), df$sample[i])) {
      message(sprintf("Row %d: tissue '%s' does not match sample '%s'", i, tissue, df$sample[i]))
      all_matched <- FALSE
    }
  }
  
  if (all_matched) {
    message("All tissues matched in samples")
  }
}
for(i in c(1:5)){
  df <- readxl::read_excel("data/samples/all/mouse_information.xlsx",sheet = i)
  print(i)
  df <- as.data.frame(df)
  check_ages(df)
  check_sample_uniqueness(df)
  check_exact_tissue_in_sample(df)
}


