#' plot_simulation_outcome
#'
#' Function to generate quality control plots for altered and unaltered counts.
#'
#' This function generates a combined plot showing the distributions of
#' various quality control metrics (such as mitochondrial and ribosomal
#' proportions) before and after damage simulation. It compares unaltered
#' counts against altered counts for a more comprehensive assessment.
#'
#' @param qc_summary A data frame containing the quality control summary for
#'  cells.
#' @param target_damage Numeric vector specifying the target damage levels
#'  for color scaling.
#' @param palette A character vector specifying the color gradient used for
#'  coloring the damage levels.
#' @return A `ggplot2` object representing the combined plot of altered
#' and unaltered counts.
#' @import ggplot2
#' @importFrom dplyr select rename mutate filter bind_rows
#' @importFrom tidyr pivot_longer
#' @importFrom scales rescale
#' @keywords internal
plot_simulation_outcome <- function(
    qc_summary,
    target_damage =  c(0.1, 0.8),
    palette = c("grey", "#7023FD", "#E60006"),
    base_size = 14
) {
  # Generate QC plots for altered & unaltered counts

  altered_counts_plot <- plot_altered_counts(
    qc_summary = qc_summary,
    palette = palette,
    target_damage = target_damage,
    base_size = base_size
  )

  unaltered_counts_plot <- plot_unaltered_counts(
    qc_summary = qc_summary,
    base_size = base_size
  )

  # Combine & return
  final_plot <- unaltered_counts_plot / altered_counts_plot

  return(final_plot)

}

#' plot_detection_outcome
#'
#' Function to generate a plot showing the distribution of quality control
#' metrics across altered data.
#'
#' This function visualizes the distribution of features and proportions of
#'  mitochondrial and ribosomal genes for altered cells, coloring the points by
#'  their damage levels. It helps in assessing how well the damage detection
#'  process has classified cells based on their quality control metrics.
#'
#' @param qc_summary A data frame containing the quality control summary for
#' cells.
#' @param target_damage Numeric vector specifying the target damage levels for
#' color scaling.
#' @param palette A character vector specifying the color gradient used for
#' coloring the damage levels.
#'
#' @return A `ggplot2` object representing the scatter plot of quality
#' control metrics.
#' @importFrom ggplot2 ggplot aes geom_point facet_wrap
#' @importFrom ggplot2 scale_y_continuous labs theme_minimal
#' @importFrom dplyr select rename mutate filter bind_rows
#' @importFrom scales rescale
#' @keywords internal
plot_detection_outcome <- function(
  qc_summary,
  target_damage =  c(0.1, 0.8),
  palette = c("grey", "#7023FD", "#E60006"),
  base_size = 14
) {
  # Isolate columns of interest
  qc_summary_long <- qc_summary %>%
    dplyr::select(features,
                  mt.prop,
                  rb.prop,
                  DamageDetective) %>%
    dplyr::rename(
      Features = features,
      `Mito. prop` = mt.prop,
      `Ribo. prop` = rb.prop
    ) %>%
    tidyr::pivot_longer(
      cols = c(Features, `Ribo. prop`),
      names_to = "X_Variable",
      values_to = "X_Value"
    )

  # Create scatter plot showing QC metric distribution
  final_plot <- ggplot2::ggplot(qc_summary_long,
                                aes(x = .data$X_Value,
                                    y = .data$`Mito. prop`,
                                    color = .data$DamageDetective)) +
    ggplot2::geom_point(alpha = 0.7, size = 0.7) +
    facet_wrap(~ .data$X_Variable,
               scales = "free_x",
               strip.position = "bottom") +
    ggplot2::scale_color_gradientn(
      colours = palette,
      values = scales::rescale(target_damage),
      guide = guide_colorbar(
        title.position = "top",
        barwidth = 10
      )
    ) +
    ggplot2::labs(x = NULL, y = "Mito. prop", color = "Damage score") +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.background = element_rect(fill = "#F5F5F5", color = NA),
      panel.grid.minor = element_blank(),
      axis.ticks = element_blank(),
      axis.title.y = element_text(face = "bold", size = base_size, vjust = 2),
      axis.title.x = element_text(face = "bold", size = base_size),
      strip.text = element_text(face = "bold", size = base_size * 0.9),
      strip.placement = "outside",
      legend.position = "bottom",
      legend.title = element_text(
        face = "bold", hjust = 0.5, vjust = 2, size = base_size),
      legend.margin = margin(t = 0, b = 0),
      legend.key.height = unit(0.5, "cm")
    )

  return(final_plot)
}

#' plot_altered_counts
#'
#' Function to generate a scatter plot of quality control metrics for altered
#' data.
#'
#' This function visualizes the distribution of features and proportions of
#' mitochondrial and ribosomal genes for cells with altered counts. It also
#' adds a reference point to help assess the quality of altered data against
#' the expected distributions.
#'
#' @param qc_summary A data frame containing the quality control summary for
#' cells.
#' @param palette A character vector specifying the color gradient used for
#' coloring the damage levels.
#' @param target_damage Numeric vector specifying the target damage levels for
#' color scaling.
#'
#' @return A list containing the plot object
#' @importFrom ggplot2 ggplot aes geom_point facet_wrap
#' @importFrom ggplot2 scale_y_continuous labs theme_minimal
#' @importFrom dplyr select rename mutate filter bind_rows
#' @importFrom tidyr pivot_longer
#' @importFrom scales rescale
#' @keywords internal
plot_altered_counts <- function(
  qc_summary,
  palette = c("grey", "#7023FD", "#E60006"),
  target_damage = c(0.5, 1),
  base_size = 14
){
  # Isolate altered counts for visualizing QC metrics
  qc_summary_long_filtered <- qc_summary %>%
    tidyr::pivot_longer(
      cols = c(Original_Features, New_Features,
               Original_MitoProp, New_MitoProp,
               Original_RiboProp, New_RiboProp),
      names_to = c("State", ".value"),
      names_pattern = "(Original|New)_(.*)"
    ) %>%
    dplyr::mutate(
      State = factor(ifelse(
        .data$State == "Original", "Original", "Altered"),
        levels = c("Original", "Altered")
      )
    ) %>%
    dplyr::rename(
      `Ribo. prop` = RiboProp
    ) %>%
    dplyr::filter(.data$State == "Altered") %>%
    tidyr::pivot_longer(
      cols = c(Features, `Ribo. prop`),
      names_to = "X_Variable",
      values_to = "X_Value"
    )

  # Create scatter plot showing QC metric distribution
  plot <- ggplot2::ggplot(qc_summary_long_filtered,
                          aes(x = .data$X_Value,
                              y = .data$MitoProp,
                              color = .data$Damaged_Level)) +
    ggplot2::geom_point(alpha = 0.7, size = 0.7) +
    facet_wrap(~ .data$X_Variable,
               scales = "free_x",
               strip.position = "bottom") +
    ggplot2::scale_color_gradientn(
      colours = palette,
      values = scales::rescale(target_damage),
      guide = guide_colorbar(
        title.position = "top",
        barwidth = 10
      )
    ) +
    # ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(x = NULL, y = "Mito. prop", color = "Damage score") +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.background = element_rect(fill = "#F5F5F5", color = NA),
      panel.grid.minor = element_blank(),
      axis.ticks = element_blank(),
      axis.title.y = element_text(face = "bold", size = base_size, vjust = 2),
      axis.title.x = element_text(face = "bold", size = base_size),
      strip.text = element_text(face = "bold", size = base_size * 0.9),
      strip.placement = "outside",
      legend.position = "bottom",
      legend.title = element_text(
        face = "bold", hjust = 0.5, vjust = 2, size = base_size),
      legend.margin = margin(t = 0, b = 0),
      legend.key.height = unit(0.5, "cm")
    )

  return(plot)

}

#' plot_unaltered_counts
#'
#' Function to generate a scatter plot of quality control metrics for unaltered
#' data.
#'
#' This function visualizes the distribution of features and proportions of
#' mitochondrial and ribosomal genes for unaltered cells. It provides a
#' reference plot to assess the original data's quality before any alterations.
#'
#' @param qc_summary A data frame containing the quality control summary for
#' cells.
#' @return A `ggplot2` object representing the scatter plot of quality control
#'  metrics for unaltered cells.
#' @importFrom ggplot2 ggplot aes geom_point facet_wrap scale_y_continuous
#'  labs theme_minimal
#' @importFrom dplyr select rename mutate filter bind_rows
#' @importFrom scales rescale
#' @keywords internal
plot_unaltered_counts <- function(
  qc_summary,
  base_size = 14
){

  # Isolate unaltered counts for visualizing QC metrics
  qc_summary_long_filtered <- qc_summary %>%
    tidyr::pivot_longer(
      cols = c(Original_Features, New_Features,
               Original_MitoProp, New_MitoProp,
               Original_RiboProp, New_RiboProp),
      names_to = c("State", ".value"),
      names_pattern = "(Original|New)_(.*)"
    ) %>%
    dplyr::mutate(
      State = factor(
        ifelse(
          .data$State == "Original", "Original", "Altered"
        ),
        levels = c("Original", "Altered")
      )
    ) %>%
    dplyr::rename(
      `Ribo. prop` = RiboProp
    ) %>%
    dplyr::filter(.data$State == "Original") %>%
    tidyr::pivot_longer(
      cols = c(Features, `Ribo. prop`),
      names_to = "X_Variable",
      values_to = "X_Value"
    )

  # Create scatter plot showing QC metric distribution
  plot <- ggplot2::ggplot(qc_summary_long_filtered,
                          aes(x = .data$X_Value, y = .data$MitoProp)) +
    ggplot2::geom_point(alpha = 0.7, size = 0.7, colour = "grey") +
    ggplot2::facet_wrap(~ X_Variable,
                        scales = "free_x",
                        strip.position = "bottom") +
    ggplot2::labs(x = NULL, y = "Mito. prop", color = "Damage score") +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.background = element_rect(fill = "#F5F5F5", color = NA),
      panel.grid.minor = element_blank(),
      axis.ticks = element_blank(),
      axis.title.y = element_text(face = "bold", size = base_size, vjust = 2),
      axis.title.x = element_text(face = "bold", size = base_size),
      strip.text = element_text(face = "bold", size = base_size * 0.9),
      strip.placement = "outside"
    )

  return(plot)
}

#' plot_ribosomal_penalty
#'
#' Function to generate a scatter plot of the simulated data focusing on
#' ribosomal proportion.
#'
#' This function visualizes the distribution of mitochondrial and ribosomal
#' proportion for cells with altered counts.
#'
#' @param qc_summary A data frame containing the quality control summary for
#' cells.
#' @param palette A character vector specifying the color gradient used for
#' coloring the damage levels.
#' @param target_damage Numeric vector specifying the target damage levels for
#' color scaling.
#'
#' @return A ggplot2 plot
#' @importFrom ggplot2 ggplot aes geom_point facet_wrap
#' @importFrom ggplot2 scale_y_continuous labs theme_minimal
#' @importFrom dplyr select rename mutate filter bind_rows
#' @importFrom tidyr pivot_longer
#' @importFrom scales rescale
#' @keywords internal
plot_ribosomal_penalty <- function(
  qc_summary,
  palette = c("grey", "#7023FD", "#E60006"),
  target_damage = c(0.5, 1),
  base_size = 14
){
  # Isolate altered counts for visualizing QC metrics
  altered_data <- qc_summary %>%
    tidyr::pivot_longer(
      cols = c(Original_MitoProp, New_MitoProp,
               Original_RiboProp, New_RiboProp),
      names_to = c("State", ".value"),
      names_pattern = "(Original|New)_(.*)"
    ) %>%
    dplyr::mutate(
      State = factor(ifelse(
        .data$State == "Original", "Original", "Altered"),
        levels = c("Original", "Altered")
      )
    ) %>%
    dplyr::rename(
      `Ribo. prop` = RiboProp
    ) %>%
    dplyr::filter(.data$State == "Altered")


  # Create scatter plot showing QC metric distribution
  plot <- ggplot2::ggplot(altered_data,
                                  aes(x = .data$`Ribo. prop`,
                                      y = .data$MitoProp,
                                      color = .data$Damaged_Level)) +
    ggplot2::geom_point(alpha = 0.7, size = 0.7) +
    ggplot2::scale_color_gradientn(
      colours = palette,
      values = scales::rescale(target_damage),
      guide = guide_colorbar(
        title.position = "top",
        barwidth = 10
      )
    ) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(x = "Ribo. prop", y = "Mito. prop", color = "Damage score") +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.background = element_rect(fill = "#F5F5F5", color = NA),
      panel.grid.minor = element_blank(),
      axis.ticks = element_blank(),
      axis.title.y = element_text(face = "bold", size = base_size, vjust = 2),
      axis.title.x = element_text(face = "bold", size = base_size),
      legend.position = "bottom",
      legend.title = element_text(
        face = "bold", hjust = 0.5, vjust = 2, size = base_size),
      legend.margin = margin(t = 0, b = 0),
      legend.key.height = unit(0.5, "cm")
    )

  return(plot)

}
