# Read data (assumed to be available in filepath)
data <- readRDS("C:/path/to/data.rds")

#---------------------------------------------#
# SENSITIVITY ANALYSIS missing data handling  #
#---------------------------------------------#

# ------- PREP FOR IMPUTATION ------- #

# Generate a complete time grid for each participant (max. 105 EMA questionnaires)
data_full <- data %>%
  group_by(label) %>%
  complete(Daynumber = 1:21, Bleepnumber = 1:5) %>%  # fill missing timepoints
  fill(label, .direction = "downup") %>%             # fill label for the new rows
  ungroup()

# Calculate TimeIndex
data_full <- data_full %>%
  mutate(TimeIndex = (Daynumber - 1) * 5 + Bleepnumber)

# Check if the TimeIndex is correct
print(range(data_full$TimeIndex))                    # Should be 1 to 105

# Save the expanded dataset
# saveRDS(data_full, "C:/path/to/data_full.rds")


# --- Helper function to filter participants by compliance threshold --- #

# Computes the compliance rate based on a chosen variable ("Fat")
# and returns the participant labels meeting the threshold.
filter_compliance <- function(data_full, threshold, var = "Fat") {
  compliant_ids <- data_full %>%
    group_by(label) %>%
    summarise(rate = sum(!is.na(.data[[var]])) / n()) %>%
    filter(rate >= threshold) %>%
    pull(label)
  return(compliant_ids)
}

# Reverse code concentration and motivation vars
data_full <- data_full %>%
  mutate(
    ConR = 10 - Con,  # Higher TdoR = more concentration problems
    TdoR = 10 - Tdo   # Higher TdoR = more motivation problems
  )

# Define thresholds
pts_00 <- unique(data_full$label)
pts_25 <- filter_compliance(data_full, 0.25, "Fat")
pts_50 <- filter_compliance(data_full, 0.50, "Fat")
pts_75 <- filter_compliance(data_full, 0.75, "Fat")

# For easier looping later, create a named list of participant IDs per threshold
compliance_lists <- list("00" = pts_00, "25" = pts_25, "50" = pts_50, "75" = pts_75)


# ------- IMPUTATION METHODS ------- #

# --- (1) no imputation reference datasets --- #
# For each threshold, simply filter data_full
noimp <- list()
for (thr in names(compliance_lists)) {
  noimp[[paste0("noimp_", thr)]] <- data_full %>% 
    filter(label %in% compliance_lists[[thr]])
}
saveRDS(noimp, "C:/path/to/noimp.rds")

# ---- (2) kalman filtering imputation datasets ----
# Define the Kalman filtering function
kalman_filtering <- function(data_imp, participant_ids) {
  # Initialize the imputed data frame
  participant_ids <- as.character(participant_ids)
  data_klmn <- data_full %>%
    mutate(label = as.character(label)) %>%
    filter(label %in% participant_ids)
  # Define optimization control parameters
  optim_control_params = list(maxit=100)
  # Apply normalization, imputation, and denormalization for each var and pts
  for (var in vars) {
    for (participant in participant_ids) {
      participant_rows <- data_klmn$label == participant
      time_series <- data_klmn[participant_rows, var]
      # check if there are non-missing values to calculate min and max
      if (all(is.na(time_series))) {
        next  # skip to the next iteration if all values are NA
      }
      min_val <- min(time_series, na.rm = TRUE)
      max_val <- max(time_series, na.rm = TRUE)
      # Avoid division by zero in case all non-NA values are the same
      if (min_val == max_val) {
        # Impute missing values with the constant value
        time_series[is.na(time_series)] <- min_val
      } else {
        # Normalize the data
        normalized_series <- (time_series - min_val) / (max_val - min_val)
        
        # Perform imputation
        imputed_series <- na_kalman(normalized_series, model = "StructTS", 
                                    optim.control = optim_control_params)
        # Denormalize the data
        time_series <- imputed_series * (max_val - min_val) + min_val
      }
      # Assign the imputed series back to the data frame
      data_klmn[participant_rows, var] <- time_series
    }
    # Clip values to be within the VAS scale range [0, 10]
    data_klmn[[var]] <- pmax(pmin(data_klmn[[var]], 10), 0)
  }
  return(data_klmn)
}

# Create scenario list for Kalman and store the datasets
klmn <- list()
for (thr in names(compliance_lists)) {
  klmn[[paste0("klmn_", thr)]] <- kalman_filtering(data_full, compliance_lists[[thr]])
}
saveRDS(klmn, "C:/path/to/klmn.rds")

#---- (3) mlMICE with TimeIndex imputation datasets ----
# TimeIndex is used instead of Daynumber, hence the "ti"
Cnow <- Sys.time()

run_mlMICE <- function(data_subset) {
  data_miceti <- data_subset %>%
    mutate(TimeIndex = (Daynumber - 1) * 5 + Bleepnumber) %>%
    select(label, TimeIndex, Bleepnumber, all_of(vars))
  
  label_lookup <- data_miceti %>% distinct(label) %>% mutate(label_factor = as.integer(as.factor(label)))
  
  data_miceti <- data_miceti %>%
    mutate(label_factor = as.integer(as.factor(label))) %>%
    select(label_factor, TimeIndex, Bleepnumber, all_of(vars))
  
  pred_matrix <- make.predictorMatrix(data_miceti)
  pred_matrix[,] <- 0
  pred_matrix[vars, vars] <- 1
  pred_matrix[vars, "label_factor"] <- -2
  pred_matrix[vars, "TimeIndex"] <- 1
  pred_matrix <- pred_matrix[vars, ]
  
  bl <- setNames(as.list(vars), vars)
  meth <- rep("2l.pan", length(vars))
  names(meth) <- vars
  
  impti <- mice(data_miceti,
                method = meth,
                predictorMatrix = pred_matrix,
                blocks = bl,
                m = 5,
                maxit = 10,
                seed = 123,
                printFlag = FALSE)
  
  data_micedti <- complete(impti, "long") %>%
    group_by(.id) %>%
    summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) %>%
    ungroup()
  
  data_micedti <- data_micedti %>%
    mutate(label_factor = ceiling(.id / 105),
           Bleepnumber = rep(1:5, times = 21, length.out = n()),
           Daynumber = rep(1:21, each = 5, length.out = n()),
           TimeIndex = (Daynumber - 1) * 5 + Bleepnumber) %>%
    left_join(label_lookup, by = "label_factor") %>%
    select(-label_factor, -.id) %>%
    relocate(label)
  
  for (v in vars) {
    data_micedti[[v]] <- pmax(pmin(data_micedti[[v]], 10), 0)
  }
  return(as.data.frame(data_micedti))
}

pb <- txtProgressBar(min = 0, max = length(compliance_lists), style = 3)
i <- 0
miceti_list <- list()
for (thr in names(compliance_lists)) {
  i <- i + 1
  data_subset <- data_full %>% filter(label %in% compliance_lists[[thr]])
  miceti_list[[paste0("mice_", thr)]] <- run_mlMICE(data_subset)
  setTxtProgressBar(pb, i)
}
close(pb)
saveRDS(miceti_list, "C:/path/to/miceti_list.rds")

# Takes some time so we keep track of it
Cres <- difftime(Sys.time(), Cnow, units = "secs")
total_seconds <- as.numeric(Cres)
minutes <- floor(total_seconds / 60)
seconds <- round(total_seconds %% 60, 2)
paste("Time for this code chunk to run:", minutes, "minutes and", seconds, "seconds")

# --- PLOTTING NETWORKS FROM THE 12 IMPUTED DATASETS ---

# Helper function to run mlVAR on a single dataset
mlvar_helper <- function(dataset,
                             vars,
                             label = "label",
                             dayvar = "Daynumber",
                             beepvar = "Bleepnumber",
                             lags = 1,
                             estimator = "lmer",
                             temporal = "orthogonal",
                             plot = TRUE,
                             name = "Dataset") {
  missing_vars <- setdiff(vars, names(dataset))
  if (length(missing_vars) > 0) {
    stop("Error: Dataset is missing required variables: ", paste(missing_vars, collapse = ", "))
  }

  if (is.null(dataset)) {
    stop("Error: Dataset is NULL.")
  }
  
  # Run mlVAR
  result <- mlVAR(
    data = dataset,
    vars = vars,
    idvar = label,
    lags = lags,
    dayvar = dayvar,
    beepvar = beepvar,
    estimator = estimator,
    temporal = temporal,
    verbose = FALSE,
    scaleWithin = TRUE
  )
  
  summary_obj <- summary(result)
  
  return(list(
    result = result,
    temporal = result$results$Beta,
    contemporaneous = summary_obj$contemporaneous
  ))
}


# --- Streamlined mlVAR pipeline ---

# Load imputed datasets
noimp <- readRDS("C:/path/to/noimp.rds")
klmn <- readRDS("C:/path/to/klmn.rds")
mice <- readRDS("C:/path/to/miceti_list.rds")

# Define method and threshold lists
methods <- list(
  noimp = noimp,
  klmn = klmn,
  mice = mice
)
method_names <- names(methods)
thresholds <- c("00", "25", "50", "75")

# Store results
networks_temp <- list()
networks_contemp <- list()

# Add a container
model_summaries <- list()

# Loop to run mlVAR on each dataset
for (method in method_names) {
  networks_temp[[method]] <- list()
  networks_contemp[[method]] <- list()

  model_summaries[[method]] <- list()
  
  for (thr in thresholds) {
    data <- methods[[method]][[paste0(method, "_", thr)]]
    
    if (is.null(data)) {
      warning(paste("Missing dataset:", method, thr))
      next
    }

    if (!"ConR" %in% names(data)) data$ConR <- 10 - data$Con
    if (!"TdoR" %in% names(data)) data$TdoR <- 10 - data$Tdo

    model_out <- mlvar_helper(
      dataset = data,
      vars = vars,
      plot = FALSE,
      name = paste(toupper(method), "-", thr, "%")
    )
    
    # Store network components
    model_summaries[[method]][[thr]] <- model_out
    networks_temp[[method]][[thr]] <- getNet(model_out$result, "temporal")
    networks_contemp[[method]][[thr]] <- getNet(model_out$result, "contemporaneous")
  }
}
saveRDS(model_summaries, "C:/path/to/model_summaries.rds")

# Pick any available model as reference
ref_model <- model_summaries[["noimp"]][["00"]]$result

# Extract fixed layout from the contemporaneous network
layout_contemp <- qgraph(getNet(ref_model, "contemporaneous"), layout = "circle", DoNotPlot = TRUE)$layout
layout_temporal <- qgraph(getNet(ref_model, "temporal"), layout = "circle", DoNotPlot = TRUE)$layout


# --- Plot all contemporaneous networks in a 3x4 grid ---
node_labels <- c(
  Fat       = "F1",        # Physical fatigue
  ConR      = "F2",        # Concentration problems
  TdoR      = "F3",        # Motivation problems
  PoAf      = "A+",        # Positive affect
  NeAf      = "A-",        # Negative affect
  Pain      = "PA",        # Pain
  PA        = "PH"         # Physical activity
)

par(mfrow = c(3, 4), mar = c(2, 2, 3, 1))

for (method in method_names) {
  for (thr in thresholds) {
    model <- model_summaries[[method]][[thr]]$result
    
    if (is.null(model)) next
    
    plot_title <- paste(toupper(method), "-", thr, "%")
    
    cmat <- getNet(model, "contemporaneous")  # Extract matrix
    
    qgraph(cmat,
           layout = layout_contemp,
           labels = node_labels[colnames(cmat)],
           title = plot_title,
           theme = "classic",
           label.prop = 1,
           fade = FALSE,
           vsize = 12,
           shape = "circle",
           edge.width = 2,
           edge.labels = FALSE)
  }
}
grid_plot_contemp <- recordPlot()

# --- Plot all temporal networks in a 3x4 grid ---
par(mfrow = c(3, 4), mar = c(2, 2, 3, 1))

for (method in method_names) {
  for (thr in thresholds) {
    model <- model_summaries[[method]][[thr]]$result
    
    if (is.null(model)) next
    
    plot_title <- paste(toupper(method), "-", thr, "%")
    
    tmat <- getNet(model, "temporal")  # Extract temporal network matrix
    
    qgraph(tmat,
           layout = layout_temporal,
           labels = node_labels[colnames(tmat)],
           title = plot_title,
           theme = "classic",
           label.prop = 1,
           fade = FALSE,
           edge.labels = FALSE,
           vsize = 12,
           shape = "circle",
           edge.width = 3)
  }
}
grid_plot_temp <- recordPlot()

file_out <- function(type) paste0("C:/path/to/networks_", type, "_", Sys.Date(), ".png")

# Save contemporaneous networks
pdf(file_out_pdf("contemporaneous"), width = 12, height = 7)
replayPlot(grid_plot_contemp)
dev.off()

# Save temporal networks
pdf(file_out_pdf("temporal"), width = 12, height = 7)
replayPlot(grid_plot_temp)
dev.off()

