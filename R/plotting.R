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
  plot_dt <- plot_dt[cn_round <= max_cn]
  
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
#' Frequency (VAF) and estimated SNV copy numbers.
#'
#' @param multiplicity_gr A GRanges object from the output of `multiplicity()`.
#' @param y_field The copy number field for y-axis
#' @param color_by_cn Logical, whether to color points by integer copy number
#' @param title The title for the plot
#' @param max_cn Maximum copy number to include in the plot
#'
#' @return A ggplot object
#' @export
plot_vaf_cn_scatter <- function(multiplicity_gr,
                               y_field = "altered_copies",
                               color_by_cn = TRUE,
                               title = "VAF vs. Altered Allele Copy Number",
                               max_cn = 10) {
  
  # Check required fields
  required_fields <- c("VAF", y_field, "cn")
  if (!all(required_fields %in% names(mcols(multiplicity_gr)))) {
    stop(paste("One or more required fields not found:", paste(required_fields, collapse = ", ")))
  }
  
  # Convert GRanges to data.table for plotting
  plot_dt <- gUtils::gr2dt(multiplicity_gr)
  plot_dt <- plot_dt[!is.na(cn) & !is.na(VAF) & !is.na(get(y_field))]
  plot_dt[, cn_round := round(cn)]
  plot_dt <- plot_dt[cn_round <= max_cn]
  
  # Create the plot
  p <- ggplot(plot_dt, aes(x = VAF, y = .data[[y_field]])) +
    geom_point(alpha = 0.5, size = 0.8) +
    labs(
      title = title,
      x = "Variant Allele Frequency (VAF)",
      y = gsub("_", " ", y_field)
    ) +
    theme_bw()
  
  if (color_by_cn) {
    p <- p + aes(color = factor(cn_round)) +
      scale_color_manual(values = as.vector(pals::glasbey(n = length(unique(plot_dt$cn_round)))), 
                          name = "Integer CN") +
      guides(color = guide_legend(override.aes = list(size = 3)))
  }
  
  return(p)
}

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
  
  if (!field %in% names(mcols(multiplicity_gr))) {
    stop(paste("Field '", field, "' not found in the input GRanges.", sep = ""))
  }
  
  # Convert GRanges to data.table for plotting
  plot_dt <- gUtils::gr2dt(multiplicity_gr)
  plot_dt <- plot_dt[!is.na(get(field)) & get(field) <= max_cn]
  
  # Create integer lines for reference
  integer_lines <- seq(0, max_cn, by = 1)
  
  # Create the plot
  p <- ggplot(plot_dt, aes(x = .data[[field]])) +
    geom_density(adjust = 1.5, fill = "steelblue", alpha = 0.5) +
    geom_vline(xintercept = integer_lines, linetype = "dashed", color = "grey40") +
    labs(
      title = title,
      x = gsub("_", " ", field),
      y = "Density"
    ) +
    xlim(0, max_cn) +
    theme_bw()
  
  if (group_by_class && "class" %in% names(plot_dt)) {
    p <- ggplot(plot_dt, aes(x = .data[[field]], fill = class)) +
      geom_density(adjust = 1.5, alpha = 0.5) +
      geom_vline(xintercept = integer_lines, linetype = "dashed", color = "grey40") +
      labs(
        title = title,
        x = gsub("_", " ", field),
        y = "Density",
        fill = "Class"
      ) +
      xlim(0, max_cn) +
      scale_fill_brewer(palette = "Set1") +
      theme_bw()
  }
  
  return(p)
}

#' @name plot_multiplicity_cn_vs_segment
#' @title Plot Copy Number Estimates vs. Segment CN
#' @description
#' Creates a scatter plot to visualize how estimated SNV copy numbers
#' relate to the underlying segment copy number.
#'
#' @param multiplicity_gr A GRanges object from the output of `multiplicity()`.
#' @param field The field to plot, e.g., "total_snv_copies"
#' @param title The title for the plot
#' @param color_by_class Logical, whether to color by variant class
#' @param alpha Point transparency
#' @param max_cn Maximum copy number to display
#'
#' @return A ggplot object
#' @export
plot_multiplicity_cn_vs_segment <- function(multiplicity_gr,
                                          field = "total_snv_copies",
                                          title = "SNV Copy Number vs. Segment Copy Number",
                                          color_by_class = TRUE,
                                          alpha = 0.5,
                                          max_cn = 10) {
  
  required_fields <- c(field, "cn")
  if (!all(required_fields %in% names(mcols(multiplicity_gr)))) {
    stop(paste("One or more required fields not found:", paste(required_fields, collapse = ", ")))
  }
  
  # Convert GRanges to data.table for plotting
  plot_dt <- gUtils::gr2dt(multiplicity_gr)
  plot_dt <- plot_dt[!is.na(get(field)) & !is.na(cn)]
  plot_dt <- plot_dt[get(field) <= max_cn & cn <= max_cn]
  
  # Create the identity line coordinates
  identity_line <- data.table(x = seq(0, max_cn, 0.1),
                             y = seq(0, max_cn, 0.1))
  
  # Create the plot
  p <- ggplot(plot_dt, aes(x = cn, y = .data[[field]])) +
    geom_point(alpha = alpha) +
    geom_line(data = identity_line, aes(x = x, y = y), 
              linetype = "dashed", color = "red") +
    labs(
      title = title,
      x = "Segment Copy Number",
      y = gsub("_", " ", field)
    ) +
    xlim(0, max_cn) +
    ylim(0, max_cn) +
    theme_bw()
  
  if (color_by_class && "class" %in% names(plot_dt)) {
    p <- p + aes(color = class) +
      scale_color_brewer(palette = "Set1")
  }
  
  return(p)
}