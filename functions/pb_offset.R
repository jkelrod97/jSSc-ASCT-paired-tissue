# This script contains the following functions
#  1) pb_offset - runs pseudobulk offset model
#  2) plot_pb_offset - produces longitudinal gene expression plots 
#     with empirical and fitted values for specified genes
#  3) plot_gene - helper function for plot_pb_offset
#  4) pbo_pred_by_gene - gets pb_offset predictions for a given gene
#  5) check_pbo_inputs - makes sure inputs to pb_offset are valid

################################################################################
# (1) PB_OFFSET ################################################################
################################################################################

# Load dependencies
library(Matrix)
library(SingleCellExperiment)
library(glmGamPoi)
library(dplyr)
library(rlang)
library(patchwork)

# Function
pb_offset <- function(seurat_obj,
                      sample_id_variable,
                      design_formula,
                      subset = NULL) {
  # PB_OFFSET:
  # Function that runs pseudobulk-offset algorithm (based off of Lee & Han, 2024)
  # for a given Seurat object, sample ID variable (for pseudobulking) and design formula
  # 
  # Arguments:
  # (1) seurat_obj (Seurat object): Seurat object.
  # (2) sample_id_variable (character): name of sample ID variable in seurat_obj.
  # (3) design_formula (formula): Formula with full model to test. Variables
  #     must be contained in seurat_obj.
  # (4) subset (call, optional): Subsetting to be applied to Seurat object. For example, 
  #  subset = Tissue == "PBMC"
  #
  # Outputs:
  # out_list (list), containing:
  #   (1) fit (glmGamPoi): entire glmGamPoi output,
  #   (2) pb (matrix): pseudobulk matrix, 
  #   (3) input (list): all function inputs
  
  # Check that inputs are valid
  check_pbo_inputs(seurat_obj,
                  sample_id_variable,
                  design_formula,
                  subset)
  
  # Subset Seurat object if requested
  original_seurat_obj <- seurat_obj
  if (!is.null(subset)) {
    seurat_obj <- subset(seurat_obj, subset = !!rlang::f_rhs(subset))
  } 
  
  # Extract metadata from seurat_obj
  meta <- seurat_obj@meta.data
  
  # Give sample_id_variable an order (convert to factor)
  meta[[sample_id_variable]] <- factor(meta[[sample_id_variable]])
  
  # Store original column names (PseudobulkExpression calls CreateCategoryMatrix, which
  # converts underscores to dashes. We do not want any column names to be altered)
  original_colnames <- levels(meta[[sample_id_variable]])
  
  # Get pseudobulk
  pb <- PseudobulkExpression(seurat_obj, 
                             group.by = sample_id_variable, 
                             return.seurat = F,
                             layer = "counts",
                             method = "aggregate") %>% 
    .[["RNA"]] %>% suppressWarnings()
  
  
  # Restore original column names
  colnames(pb) <- original_colnames

  # Create colData
  covariates <- all.vars(design_formula)
  colData <- meta %>%
    dplyr::select(all_of(c(sample_id_variable, covariates))) %>%
    distinct()
  
  # Reorder PB columns to match order in colData
  pb <- pb[, colData[[sample_id_variable]]]
  
  # Store data in a SingleCellExperiment object
  sce <- SingleCellExperiment(assays = list(counts = pb))
  # Add colData to sce object (and transform to DataFrame object)
  colData(sce) <- DataFrame(colData)
  
  # Create design matrix
  design_matrix <- model.matrix(design_formula, data = colData)
  
  # Fit model (using glmGamPoi's Gamma-Poisson GLM algorithm)
  fit <- glm_gp(sce, design = design_matrix, on_disk = F) %>% suppressWarnings()
  
  # Return list containing 
  #   (1) fit (glmGamPoi): entire glmGamPoi output,
  #   (2) pb (matrix): pseudobulk matrix, 
  #   (3) seurat_obj (Seurat): Seurat object, with subsetting applied
  #   (4) input (list): 
  #       (a) original_seurat_obj (Seurat object--no subsetting)
  #       (b) sample_id_variable (string)
  #       (c) design_formula (formula)
  #       (d) subset (call)
  
  out_list <- list(fit = fit,
                   pb = pb,
                   seurat_obj = seurat_obj,
                   input = list(original_seurat_obj = original_seurat_obj,
                                sample_id_variable = sample_id_variable,
                                design_formula = design_formula,
                                subset = subset))
  return(out_list)
}

################################################################################
# (2) PLOT_PB_OFFSET ###########################################################
################################################################################

plot_pb_offset <- function(pbo, 
                           genes_to_plot,
                           separate_plots_by = NULL,
                           groupBy = NULL,
                           new_timepoints = NULL,
                           scale = "log2",
                           linewidth = 2,
                           legend.pos = "bottom",
                           title_text_size = 15,
                           axis_text_size = 20,
                           axis_label_size = 17,
                           legend_text_size = 15,
                           legend_title_size = 20
                           ) {
  # PLOT_PB_OFFSET
  # Plots empirical gene expression by subject over time for top genes
  # identified by pb_offset Can also plot expression predicted
  # by the model. 
  # 
  # Arguments:
  # (1) pbo (list): output from pb_offset
  # (2) genes_to_plot (vector of strings): list of genes to plot. 
  # (3) separate_plots_by (string): variable to separate plots by
  # (3) groupBy (string): variable to group by
  # (4) new_timepoints (numeric): new x axis values for predictions
  # (5) scale (string): log2 (default), log (ln), or original
  # (6) linewidth (numeric): linewidth parameter in ggplot
  # (7) legend.pos (string): legend.position in ggplot
  # (8) title_text_size (numeric): title text size
  # (9) axis_text_size (numeric): text size for axis values
  # (10) axis_label_size (numeric): text size for axis labels
  # (11) legend_text_size (numeric): text size for legend values
  # (12) legend_title_size (numeric): text size for legend title

  # Output:
  # plot_list (list): list containing longitudinal expression plots 
  # for specified genes
  
  # Extract pb_offset output
  # (1) pseudobulk matrix
  pb <- pbo$pb
  # (2) sample_id_variable
  sample_id_variable <- pbo$input$sample_id_variable
  # (3) seurat_obj
  seurat_obj <- pbo$seurat_obj
  # (4) metadata
  meta <- seurat_obj@meta.data
  
  # Subset pseudobulk matrix to only include desired genes
  pb_sub <- pb[genes_to_plot, ]
  
  # If only plotting one gene, make sure pb_sub is still a matrix with 1 column per sample
  if (length(genes_to_plot) == 1) {
    pb_sub <- t(as.matrix(pb_sub))
  }

  # Get longitudinal plot list, running plot_gene for all genes in genes_to_plot.
  long_plot_list <- mapply(plot_gene, 
                      genes_to_plot, 
                      MoreArgs = list(pbo = pbo, 
                                      separate_plots_by = separate_plots_by,
                                      new_timepoints = new_timepoints,
                                      scale = scale, 
                                      groupBy = groupBy,
                                      linewidth = linewidth,
                                      legend.pos = legend.pos,
                                      title_text_size = title_text_size,
                                      axis_text_size = axis_text_size,
                                      axis_label_size = axis_label_size,
                                      legend_text_size = legend_text_size,
                                      legend_title_size = legend_title_size),
                      SIMPLIFY = FALSE)
  # Get dot plot list
  # We precompute the split based on separate_plots_by
  seurat_list <- SplitObject(pbo$input$original_seurat_obj, split.by = separate_plots_by)
  dot_plot_list <- lapply(genes_to_plot,
                          dot_by_gene,
                          pbo = pbo,
                          precomputed_split = seurat_list)
  names(dot_plot_list) <- genes_to_plot
  
  # Put plots together--1 for each gene
  # (1) Combine longitudinal plots (for all values of separate_plots_by. Right
  # now can only accomodate 2 values (e.g., skin and PBMC for separate_plot_by = Tissue))
  # Display plots together
  long_plots_combined <-
    lapply(genes_to_plot,
           function (gene) {
             # One shared legend.
             long_plot_list[[gene]]$PBMC + theme(legend.position = "none") + long_plot_list[[gene]]$skin + theme(legend.position = "right")
           }
    )
  # Set genes as plot names
  names(long_plots_combined) <- genes_to_plot
  # (2) Combined dot plots
  # Dot size legend data - horizontal layout
  dot_legend <- data.frame(
    # x values of dots
    x = c(0.5, 1.5, 2.5),
    # y values of dots
    y = rep(7.5, 3),
    # size of dots
    size = c(2, 5, 9),
    # label smallest and largest dot
    label = c("small %", "", "large %")
  )
  # Color gradient bar - horizontal
  color_bar <- data.frame(
    # x values of color bar
    x = seq(0.2, 2.8, length.out = 100),
    # y values of color bar
    y = rep(1.5, 100),
    fill = seq(0, 1, length.out = 100)
  )
  # Create the legend plot
  legend_plot <- ggplot() +
    # Dot size legend
    geom_point(data = dot_legend, aes(x = x, y = y, size = size), 
               color = "black", stroke = 0.5, shape = 21, fill = "grey80") +
    ylim(0, 10) +
    # Label dots
    geom_text(data = dot_legend, aes(x = x, y = y - 2, label = label), 
              hjust = 0.5, size = 5) +
    # Legend title
    annotate("text", x = 1.5, y = 10, label = "% cells expressing", 
             hjust = 0.5, size = 7, fontface = "bold") +
    
    # Color gradient bar
    geom_tile(data = color_bar, aes(x = x, y = y, fill = fill), 
              height = 0.75, width = 0.03) +
    
    # Color bar labels
    annotate("text", x = 0.2, y = 0.1, label = "Low", hjust = 0.5, size = 5) +
    annotate("text", x = 2.8, y = 0.1, label = "High", hjust = 0.5, size = 5) +
    # Legend title
    annotate("text", x = 1.5, y = 3, label = "Average expression", 
             hjust = 0.5, size = 7, fontface = "bold") +
    # Formatting
    scale_size_identity() +
    # Set colors
    scale_fill_gradient(low = "yellow", high = "red", guide = "none") +
    # Set xlim
    xlim(-0.2, 3.2) +
    # Set theme
    theme_void() +
    # White box for legend
    theme(
      panel.background = element_rect(fill = "white", color = "white")
    )
  # Combine plots
  dot_plots_combined <- lapply(dot_plot_list, function(plots){
    combined <- plots$PBMC + plots$skin + legend_plot +
      # Specify proportions
      plot_layout(ncol = 3, widths = c(5, 5, 3))})
  # (3) Combine longitudinal and dot plots
    combined_plots <- 
      mapply(function(gene){
        # Combine plots
        combined_plot <- long_plots_combined[[gene]] / plot_spacer() / dot_plots_combined[[gene]] +
          # Set proportinos
          plot_layout(heights = c(4, 0.2, 1.5)) +
          # Set overall title
          plot_annotation(paste0(gene, " expression in jSSc patients following transplant"), 
                          theme = theme(plot.title = element_text(hjust = 0.5,
                                                                  size = 33,
                                                                  face = "bold")))
        # Return plot
        return(combined_plot)
      }, genes_to_plot)

  # Return combined dot plots
  return(combined_plots)
  
}

################################################################################
# (3) PLOT_GENE ################################################################
################################################################################

plot_gene <- function(gene,
                      pbo,
                      separate_plots_by = NULL,
                      new_timepoints = NULL,
                      scale = "log2", 
                      groupBy = NULL,
                      linewidth = 2,
                      legend.pos = legend.pos,
                      title_text_size = 20,
                      axis_text_size = 10,
                      axis_label_size = 10,
                      legend_text_size = 15,
                      legend_title_size = 20) {
  # PLOT_GENE
  # Helper for plot_pb_offset
  # Plots empirical gene expression by subject over time for a single gene,
  # along with fitted values from pb_offset.
  #
  # Arguments:
  #  (1) gene (string): gene name
  #  (2) pbo (list): output from pb_offset
  #  (3) separate_plots_by (string): variable to separate plots by (e.g., Tissue)
  #  (4) new_timepoints (numeric vector): new values to use for predictions
  #  (5) scale (string): log2 (default) or original
  #  (6) groupBy: Variable to group by (separate predictions for each group)
  #  (7) linewidth (numeric): linewidth parameter in ggplot
  #  (8) legend.pos (string): legend.position in ggplot
  #  (9) title_text_size (numeric): title text size
  #  (10) axis_text_size (numeric): text size for axis values
  #  (11) axis_label_size (numeric): text size for axis labels
  #  (12) legend_text_size (numeric): text size for legend values
  #  (13) legend_title_size (numeric): text size for legend title
  # Output:
  # gene_plot (ggplot): longitudinal expression plot for a given gene
  
  # Extract objects from pbo output.
  # (1) Data for plotting
  plot_data <- pbo$fit$data@colData
  # (2) Pseudobulk matrix
  pb <- pbo$pb
  # (3) Cell type
  cell_type <- pbo$input$cell_type
  # (4) Scale factor
  scale_factor <- exp(mean(log(colSums(pb))))
  # (5) Sample ID variable
  sample_id_variable <- pbo$input$sample_id_variable
  # (6) Seurat object
  seurat_obj <- pbo$seurat_obj
  # (7) Seurat metadata
  meta <- seurat_obj@meta.data
  # (8) Model covariates
  covariates <- all.vars(pbo$input$design_formula)
  
  # Extract row from pb
  pb_row <- pb[gene, ]
  
  # Make sure that the order of the sample IDs match
  if (any(plot_data[[sample_id_variable]] != names(pb_row))) {
    stop("Sample ID orders do not match.")
  }
  
  # Convert to CPM
  cpm <- pb_row / colSums(pb) * 1e6
  
  # Add counts (CPM scale) to plot data
  plot_data$counts <- cpm
  
  # Adjust y axis label to reflect scale
  if (scale == "log2") {
    y_axis_label <- "Counts per million (log2 scale)"
  } else if (scale == "log") {
    y_axis_label <- "Counts per million (natural logrithm scale)"
  } else {
    y_axis_label <- "Counts per million"
  }
  
  # If new_timepoints have not been specified, autogenerate them
  if (is.null(new_timepoints)) {
    # Lowest observed value
    time_min <- min(na.omit(unique(meta$timepoint)))
    # Highest observed value
    time_max <- max(na.omit(unique(meta$timepoint)))
    # New values are generated as increments of 1 between lowest and highest observed value.
    new_timepoints <- seq(time_min, time_max)
  }
  
  # Before creating the plots, we ensure consistent ordering of the groups
  plot_data[[groupBy]] <- factor(plot_data[[groupBy]], levels = sort(unique(plot_data[[groupBy]])))
  
  # Mark plot_data (only includes empirical values) as "measured"
  plot_data$type <- "measured"
  
  # Generate fitted values
  predicted_df <- pbo_pred_by_gene(gene, pbo, new_timepoints)
  
  # Convert fitted values to CPM scale
  pred_cpm_df <- predicted_df %>%
    mutate(pred_cpm = (counts / scale_factor) * 1e6)
  
  # Store CPM in counts variable
  pred_cpm_df$counts <- pred_cpm_df$pred_cpm; pred_cpm_df$pred_cpm <- NULL
  # Mark pred_cpm_df (only includes fitted values) as "predicted"
  pred_cpm_df$type <- "predicted"
  # Remove variables not in pred_cpm_df
  plot_data <- plot_data[, colnames(pred_cpm_df)]
  # Combine measured and predicted values into one data frame for plotting
  all_df <- rbind(plot_data, pred_cpm_df)
  
  # Separate plots?
  if (!is.null(separate_plots_by)) {
    all_df_list <- split(all_df, all_df[[separate_plots_by]])
  } else {
    all_df_list <- list(all_df)
  }
  
  # Plot it (one plot for each data frame)
  plots <-
    lapply(all_df_list, function(df){
      gene_plot <- ggplot(df, aes(x = timepoint, y = counts, col = !!sym(groupBy), linetype = type, alpha = type)) + 
        geom_line(linewidth = linewidth) +
        theme_minimal() +
        # Set plot labels.
        labs(y = y_axis_label, x = "Months since transplant") +
        # Set x axis breaks
        scale_x_continuous(breaks = na.omit(unique(meta$timepoint)) %>% sort) +
        # Set legend position
        theme(legend.position = legend.pos) +
        # Set linetype values by type.
        scale_linetype_manual(values = c(measured = "solid", predicted = "22"),
                              labels = c(measured = "measured", predicted = "estimated")) +
        # Set alpha values by type. We want fitted values to be darker than observed values.
        scale_alpha_manual(values = c(measured = 0.4, predicted = 1),
                           labels = c(measured = "measured", predicted = "estimated"),
                           guide = "none") +
        # Set legend label
        labs(linetype = "count type") +
        # Make linewidth in legend smaller (so that difference between dashed and solid line is visible)
        guides(linetype = guide_legend(override.aes = list(linewidth = 1))) 
      
      # Adjust scale if specified.
      if (scale == "log2") {
        # Custom scale for log2p1
        log2p1_trans <- scales::trans_new(
          name = "log2p1",
          transform = function(x) log2(x + 1),
          inverse = function(x) 2^x - 1
        )
        # Set breaks
        log2p1_breaks <- function(limits) {
          breaks <- scales::log_breaks(base = 2)(limits + 1) - 1
          breaks[breaks > 0]
        }
        # Set labels
        log2p1_labels <- function(breaks) {
          sprintf("%.0f", breaks)
        }
        # Add scale to plot
        gene_plot <- gene_plot + scale_y_continuous(trans = log2p1_trans,
                                                    breaks = log2p1_breaks,
                                                    labels = log2p1_labels)
      } else if (scale == "log") {
        # Use built-in log1p scale
        gene_plot <- gene_plot + scale_y_continuous(trans = "log1p")
      }
      
      # Adjust aesthetics entered by user
      gene_plot <- gene_plot + 
        theme(plot.title = element_text(size = title_text_size,
                                        face = "bold"),
              axis.text = element_text(size = axis_text_size),
              axis.title = element_text(size = axis_label_size),
              legend.text = element_text(size = legend_text_size),
              legend.title = element_text(size = legend_title_size,
                                          face = "bold"))
      
      
      gene_plot
    })
  
  # Add titles to plots
  if (length(all_df_list) > 1) {
    titles <- paste0(gene, ", ", names(plots))
    plots <- 
      mapply(function(plot, title){plot + ggtitle(title)},
             plots,
             titles,
             SIMPLIFY = F)
  } else {
    plots <- 
      mapply(function(plot){plot + ggtitle(gene)},
             plots,
             SIMPLIFY = F)
  }

  # Return plots
  return(plots)
}

################################################################################
# (4) DOT_BY_GENE ##############################################################
################################################################################

dot_by_gene <- function(gene, 
                        pbo,
                        subset = NULL,
                        separate_plots_by = NULL,
                        precomputed_split = NULL) {
  # DOT_BY_GENE
  # Function to produce dot plots showing avg. expression & pct expressed
  # for a given gene. Optionally, separate plots will be created according
  # to "separate_plots_by" variable (for example, Tissue)
  #
  # Arguments:
  #  (1) gene (string): gene name
  #  (2) pbo (list): output from pb_offset
  #  (3) subset (call): call for subsetting Seurat object
  #  (4) separate_plots_by (string): variable to separate plots by
  #  (5) precomputed_split (list of Seurat objects): for efficiency,
  #      you can split the Seurat object ahead of time if running for many genes
  
  # Only do these steps if split has not happened yet
  if (is.null(precomputed_split)) {
    # Get original seurat object (no subsetting) from pbo output
    seurat_obj <- pbo$input$original_seurat_obj
    
    # Now, subset if specified
    if (!is.null(subset)) {
      seurat_obj <- subset(seurat_obj, subset = !!rlang::f_rhs(subset))
    } 
    
    # Check separate_plots_by and separate accordingly
    if (!is.null(separate_plots_by)) {
      seurat_list <- SplitObject(seurat_obj, split.by = separate_plots_by)
    } else {
      seurat_list <- list(seurat_obj)
    }
  } else {
    seurat_list <- precomputed_split
  }
    
  # Create dot plots
  dot_plots <- lapply(seurat_list,
                      function(obj){
                        DotPlot(obj, 
                                features = gene, 
                                group.by = "group") + 
                          # This line makes sure the smallest dots are still visible
                          scale_size(range = c(5, 25)) +
                          # Add border around dots (helps to see yellow dots more easily)
                          geom_point(aes(size = pct.exp, color = avg.exp.scaled), 
                                     shape = 21, fill = NA, stroke = 0.5, color = "black") + 
                          # Add dashed line for separation between healthy controls and patients
                          geom_hline(yintercept = 1.5, linetype = "dashed") + 
                          # Remove x and y labels
                          labs(y = "", x = "") + 
                          # Flip coordinates
                          coord_flip() +
                          # Set shared color scale
                          scale_colour_gradient(low = "yellow", high = "red", 
                                                name = "Average expression\n(scaled)") +
                          # Set y-axis labels
                          scale_y_discrete(labels = c(
                            "healthy" = "healthy \ncontrols",
                            "baseline" = "SSc \nbaseline",
                            "six" = "SSc \n6 mo.",
                            "twelve" = "SSc \n12 mo.",
                            "twen4" = "SSc \n24 mo."
                          ), expand = expansion(mult = c(0.2, 0.2))) +
                          # Improve layout with a single theme call
                          theme(
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            axis.text.y = element_blank(),
                            axis.text.x = element_text(margin = margin(t = 5), size = 15),
                            axis.ticks.length = unit(2, "pt"),
                            legend.position = "none"
                          )
                      })
  # Add titles
  SEPVARS <- names(dot_plots)
  dot_plots <- lapply(SEPVARS, function(SEPVAR){dot_plots[[SEPVAR]] +
                        ggtitle(paste0(gene, ", ", SEPVAR))})
  names(dot_plots) <- SEPVARS
  dot_plots
}
  
  


################################################################################
# (5) PBO_PRED_BY_GENE #########################################################
################################################################################

pbo_pred_by_gene <- function(gene, pbo, new_timepoints) {
  # PBO_PRED_BY_GENE
  # Function to make predictions for unmeasured time points with 
  # pb_offset output
  #
  # Arguments:
  # (1) gene (string): gene of interest
  # (2) pbo (list): output from pb_offset
  # (3) new_timepoints (numeric vector): new x values to use for predictions
  # pred_Data (data.frame): data frame with predictions for new x values
  
  # Extract pb_offset output
  # (1) colData
  colData <- pbo$fit$data@colData
  # (2) beta
  beta <- pbo$fit$Beta
  # (3) pb
  pb <- pbo$pb
  # (4) design formula
  design_formula <- pbo$input$design_formula

  # Generate colData for all covariate and new x value combos
  covariates <- all.vars(design_formula)
  # Get covariates without time (since we are generating new values)
  covariates_minus_time <- covariates[covariates != "timepoint"]
  pred_Data <- expand.grid(
    lapply(covariates_minus_time, function(x) {
      # 1. Get the original column
      original_col <- colData[[x]]
      
      # 2. Get the unique values from that column
      unique_values <- unique(original_col)
      
    }) %>% 
      setNames(covariates_minus_time) %>% 
      append(list(timepoint = new_timepoints))
  )

  # Get beta for current gene
  beta_gene <- beta[gene, ]
  
  # Prediction design matrix
  pred_design_matrix <- model.matrix(design_formula, data = pred_Data)
  
  # Compute predictions in log space (matrix multiplication)
  # Offset term is zero on average, so we set to zero for an unobserved new sample
  # Size factors are normalized as in glm_gp ("normalized_sum")
  log_predictions <- pred_design_matrix %*% matrix(beta_gene)
  
  # Convert back to original scale
  predictions <- exp(log_predictions)
  
  # Combine into a data frame
  pred_Data$counts <- predictions
  
  # Make note that these counts are predictions
  pred_Data$type <- "predicted"
  
  # Console output
  print(paste0("Generated predictions for gene ", gene))
  
  # Return result
  return(pred_Data)
}



################################################################################
# (6) CHECK_PBO_INPUTS #########################################################
################################################################################

check_pbo_inputs <- function(seurat_obj,
                            sample_id_variable,
                            design_formula,
                            subset = NULL) {
  # CHECK_PBO_INPUTS
  # Helper function to check that pb_offset inputs are valid
  # 
  # Returns error if any inputs are not as expected.
  
  # (1) Check seurat_obj
  if (!(inherits(seurat_obj, "Seurat"))) {
    stop("pb_offset requires a Seurat object as input.")
  }
  # Pull out the metadata
  meta <- seurat_obj@meta.data
  # (2) Check sample_id_variable
  if (!(is.character(sample_id_variable) & length(sample_id_variable) == 1)) {
    stop("sample_id_variable must be a character of length 1.")
  }
  if (!(sample_id_variable %in% names(meta))) {
    stop("sample_id_variable not found in names(seurat_obj@meta.data)")
  }  
  # (3) Check design formula
  if (class(design_formula) != "formula") {
    stop("design_formula must be a formula")
  }
  design_variables <- all.vars(design_formula)
  if (any(!(design_variables %in% names(meta)))) {
    stop("At least one design variable not found in seurat metadata.")
  }
  # Make sure each design variable has only one unique value per sample
  # UNLESS it's tissue (we can have two unique values per sample)
  for (design_variable in design_variables) {
    if (design_variable != "Tissue") {
    result <- meta %>%
      group_by(.data[[sample_id_variable]]) %>%
      summarise(n_unique = n_distinct(.data[[design_variable]])) %>%
      filter(n_unique > 1)
    
    if (nrow(result) > 0) {
      stop(paste0(design_variable, " has more than one unique value per sample"))
    }
    }
  }
  # (4) Check subset variable
  # Only check if not null (since variable is optional)
  if (!is.null(subset)) {
    if (class(subset) != "formula") {
      stop("subset must be a formula. For example, subset = ~health == 'SSC'")
    }
  }
}

