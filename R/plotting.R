#' @import ggplot2
#' @import data.table
#' @import gUtils
#' @import gTrack
#' @importFrom pals glasbey

#' @name plot_multiplicity_histogram
#' @title Plot Multiplicity Histograms
#' @description
#' Creates histograms to visualize the distribution of estimated SNV copy numbers
#' from the output of the `multiplicity()` function.
#'
#' @param multiplicity_gr A GRanges object from the output of `multiplicity()`.
#' @param field The field to plot, e.g., "total_snv_copies", "major_snv_copies", "minor_snv_copies"
#' @param facet_by_cn Logical, whether to facet the histograms by integer copy number
#' @param title The title for the plot
#' @param max_cn Maximum copy number to include in the plot
#' @param binwidth Width of histogram bins
#'
#' @return A ggplot object
#' @export
plot_multiplicity_histogram <- function(multiplicity_gr,
                     field = "total_snv_copies",
                     facet_by_cn = TRUE,
                     title = "SNV Copy Number Distribution",
                     max_cn = 10,
                     binwidth = 0.1) {
  if (!field %in% names(mcols(multiplicity_gr))) {
  stop(paste("Field '", field, "' not found in the input GRanges.", sep = ""))
  }

  # Convert GRanges to data.table for plotting
  plot_dt <- gUtils::gr2dt(multiplicity_gr)
  plot_dt <- plot_dt[!is.na(cn) & !is.na(get(field))]
  plot_dt[, cn_round := round(cn)]
  plot_dt <- plot_dt[cn_round <= max_cn & get(field) <= max_cn]
  
  # Generate integer lines for reference
  integer_lines <- seq(0, max_cn, by = 1)
  
  # Create the plot
  p <- ggplot(plot_dt, aes(x = .data[[field]], fill = factor(cn_round))) +
  geom_histogram(binwidth = binwidth, alpha = 0.8, position = "identity") +
  geom_vline(xintercept = integer_lines, linetype = "dashed", color = "grey40") +
  labs(
    title = title,
    x = gsub("_", " ", field),
    y = "Frequency",
    fill = "Integer CN"
  ) +
  scale_fill_manual(values = as.vector(pals::glasbey(n = length(unique(plot_dt$cn_round))))) +
  xlim(0, max_cn) +
  theme_bw() +
  theme(legend.position = "bottom")
  
  if (facet_by_cn) {
  p <- p + facet_wrap(~cn_round, scales = "free_y")
  }
  
  return(p)
}

#' @name plot_vaf_cn_scatter
#' @title Plot VAF vs. SNV Copy Number
#' @description
#' Creates a scatter plot to visualize the relationship between Variant Allele
#' Frequency (VAF) and estimated SNV copy numbers or altered fraction.
#'
#' @param multiplicity_gr A GRanges object from the output of `multiplicity()`.
#' @param plot_type Either "altered_copies" or "altered_fraction" to determine what to plot on y-axis
#' @param color_by_cn Logical, whether to color points by integer copy number
#' @param facet_by_cn Logical, whether to facet the plot by integer copy number
#' @param title The title for the plot (if NULL, a default title based on plot_type will be used)
#' @param max_cn Maximum copy number to include in the plot
#'
#' @return A ggplot object
#' @export
plot_vaf_cn_scatter <- function(multiplicity_gr,
  plot_type = c("altered_copies", "altered_fraction"),
  color_by_cn = TRUE,
  facet_by_cn = TRUE,
  title = NULL,
  max_cn = 10) {
  
  plot_type <- match.arg(plot_type)
  
  # Set default title based on plot_type if not provided
  if (is.null(title)) {
  title <- if (plot_type == "altered_copies") {
    "VAF vs. Altered Allele Copy Number"
  } else {
    "VAF vs. Altered Allele Fraction"
  }
  }
  
  # Check required fields
  required_fields <- c("VAF", "altered_copies", "cn")
  if (plot_type == "altered_fraction") {
  required_fields <- c(required_fields, "total_snv_copies")
  }
  
  if (!all(required_fields %in% names(mcols(multiplicity_gr)))) {
  stop(paste("One or more required fields not found:", paste(required_fields, collapse = ", ")))
  }
  
  # Convert GRanges to data.table for plotting
  plot_dt <- gUtils::gr2dt(multiplicity_gr)
  plot_dt <- plot_dt[!is.na(cn) & !is.na(VAF) & !is.na(altered_copies)]
  
  # Calculate altered_fraction if needed
  if (plot_type == "altered_fraction") {
  plot_dt[, altered_fraction := altered_copies / total_snv_copies]
  plot_dt <- plot_dt[!is.na(altered_fraction)]
  y_field <- "altered_fraction"
  } else {
  y_field <- "altered_copies"
  }
  
  plot_dt[, cn_round := round(cn)]
  plot_dt <- plot_dt[cn_round <= max_cn]
  
  # Calculate correlation for each copy number
  cor_dt <- plot_dt[, .(correlation = round(cor(VAF, get(y_field), use = "complete.obs"), 3)), by = cn_round]
  plot_dt <- merge(plot_dt, cor_dt, by = "cn_round")
  
  # Create the plot
  p <- ggplot(plot_dt, aes(x = VAF, y = .data[[y_field]])) +
  geom_point(alpha = 0.5, size = 0.8) +
  geom_smooth(method = "lm", se = FALSE, size = 0.8) +
  geom_text(data = unique(plot_dt[, .(cn_round, correlation)]), 
     aes(x = 0.2, y = 0.9 * ifelse(y_field == "altered_fraction", 1, max(plot_dt[[y_field]])), 
      label = paste0("r = ", correlation)), 
     hjust = 0, size = 3) +
  labs(
    title = title,
    x = "Variant Allele Frequency (VAF)",
    y = gsub("_", " ", y_field)
  ) +
  theme_bw()
  
  # Apply appropriate scale transformations based on plot type
  if (plot_type == "altered_fraction") {
    p <- p + ylim(0, 1)
  } else if (plot_type == "altered_copies") {
    p <- p + scale_y_continuous(trans = "log1p")
  }
  
  if (color_by_cn) {
  p <- p + aes(color = factor(cn_round)) +
    scale_color_manual(values = as.vector(pals::glasbey(n = length(unique(plot_dt$cn_round)))), 
       name = "Integer CN") +
    guides(color = guide_legend(override.aes = list(size = 3)))
  }
  
  if (facet_by_cn) {
  p <- p + facet_wrap(~cn_round, scales = "free_y")
  }
  
  return(p)
}
#}

#' @name plot_multiplicity_gtrack
#' @title Create a gTrack for Multiplicity Results
#' @description
#' Generates a gTrack object to visualize multiplicity results along the genome,
#' which can be plotted alongside other genomic data.
#'
#' @param multiplicity_gr A GRanges object from the output of `multiplicity()`.
#' @param name The name for the track
#' @param y_fields Fields to plot on the y-axis
#' @param ... Additional arguments passed to gTrack
#'
#' @return A gTrack object
#' @export
plot_multiplicity_gtrack <- function(multiplicity_gr,
                                    name = "Multiplicity",
                                    y_fields = c("total_snv_copies", "major_snv_copies", "minor_snv_copies"),
                                    ...) {
  
  # Check that y_fields exist in the GRanges object
  if (!all(y_fields %in% names(mcols(multiplicity_gr)))) {
    missing_fields <- setdiff(y_fields, names(mcols(multiplicity_gr)))
    stop(paste("Some y_fields not found in GRanges:", paste(missing_fields, collapse = ", ")))
  }
  
  # Create the gTrack
  gt <- gTrack::gTrack(
    multiplicity_gr,
    y.field = y_fields,
    name = name,
    ...
  )
  
  return(gt)
}

#' @name plot_multiplicity_cn_density
#' @title Plot Density of Copy Number Estimates
#' @description
#' Creates a density plot to visualize the distribution of estimated copy numbers
#' across different segments.
#'
#' @param multiplicity_gr A GRanges object from the output of `multiplicity()`.
#' @param field The field to plot, e.g., "total_snv_copies", "major_snv_copies", "minor_snv_copies"
#' @param group_by_class Logical, whether to group by variant class
#' @param title The title for the plot
#' @param max_cn Maximum copy number to include in the plot
#'
#' @return A ggplot object
#' @export
plot_multiplicity_cn_density <- function(multiplicity_gr,
                                         field = "total_snv_copies",
                                         group_by_class = TRUE,
                                         title = "Copy Number Density Distribution",
                                         max_cn = 10) {
  required_fields <- c(field, "cn")
  if (!all(required_fields %in% names(mcols(multiplicity_gr)))) {
    stop(paste("One or more required fields not found:", paste(required_fields, collapse = ", ")))
  }

  # Convert GRanges to data.table for plotting
  plot_dt <- gUtils::gr2dt(multiplicity_gr)
  plot_dt <- plot_dt[!is.na(get(field)) & !is.na(cn) & get(field) <= max_cn]
  plot_dt[, cn_round := round(cn)]
  plot_dt <- plot_dt[cn_round <= max_cn]

  # Create integer lines for reference
  integer_lines <- seq(0, max_cn, by = 1)

  # Create the plot colored by segment CN with contours
  p <- ggplot(plot_dt, aes(x = .data[[field]], fill = factor(cn_round), color = factor(cn_round))) +
    geom_density(adjust = 1.5, alpha = 0.3) +
    stat_density(geom = "line", adjust = 1.5, size = 0.8, alpha = 0.8, position = "identity") +
    geom_vline(xintercept = integer_lines, linetype = "dashed", color = "grey40") +
    labs(
      title = title,
      x = gsub("_", " ", field),
      y = "Density",
      fill = "Segment CN",
      color = "Segment CN"
    ) +
    xlim(0, max_cn) +
    scale_fill_manual(values = as.vector(pals::glasbey(n = length(unique(plot_dt$cn_round))))) +
    scale_color_manual(values = as.vector(pals::glasbey(n = length(unique(plot_dt$cn_round))))) +
    theme_bw() +
    theme(legend.position = "bottom")

  return(p)
}

#' @name plot_multiplicity_cn_vs_segment
#' @title Plot Copy Number Estimates vs. Segment CN
#' @description
#' Creates a violin plot with jittered points and boxplots without outliers to visualize how estimated SNV copy numbers
#' relate to the underlying segment copy number. Also computes and displays the linear correlation.
#'
#' @param multiplicity_gr A GRanges object from the output of `multiplicity()`.
#' @param field The field to plot, e.g., "total_snv_copies"
#' @param title The title for the plot
#' @param color_by_class Logical, whether to color by variant class
#' @param alpha Point transparency
#' @param max_cn Maximum copy number to display
#' @param jitter_width Width of the jitter for points
#' @param jitter_alpha Transparency for jittered points
#'
#' @return A ggplot object
#' @export
plot_multiplicity_cn_vs_segment <- function(multiplicity_gr,
            field = "total_snv_copies",
            title = "SNV Copy Number vs. Segment Copy Number",
            color_by_class = FALSE,
            alpha = 0.5,
            max_cn = 10,
            jitter_width = 0.2,
            jitter_alpha = 0.05) {
  
  required_fields <- c(field, "cn")
  if (!all(required_fields %in% names(mcols(multiplicity_gr)))) {
  stop(paste("One or more required fields not found:", paste(required_fields, collapse = ", ")))
  }
  
  plot_dt <- gUtils::gr2dt(multiplicity_gr)
  plot_dt <- plot_dt[!is.na(get(field)) & !is.na(cn)]
  plot_dt <- plot_dt[get(field) <= max_cn & cn <= max_cn]
  
  plot_dt[, cn_discrete := round(cn)]
  
  cor_value <- cor(plot_dt[[field]], plot_dt$cn, use = "complete.obs")
  cor_text <- sprintf("Correlation: r = %.3f", cor_value)
  
  # Create the plot
  p <- ggplot(plot_dt, aes(x = factor(cn_discrete), y = .data[[field]])) +
  geom_violin(alpha = 0.3, scale = "width", fill = "lightblue") +
  geom_jitter(width = jitter_width, alpha = jitter_alpha, size = 0.5, color = "black") +
  geom_boxplot(outlier.shape = NA, alpha = 0, width = 0.3, 
         color = "red", fill = "red", alpha = 0.1,
         outlier.colour = "black", 
         fatten = 2) +
  annotate("text", x = 1, y = max_cn * 0.95, label = cor_text, hjust = 0, size = 3.5) +
  labs(
    title = title,
    x = "Segment Copy Number",
    y = gsub("_", " ", field)
  ) +
  ylim(0, max_cn) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p)
}







