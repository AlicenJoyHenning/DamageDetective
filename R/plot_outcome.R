#' plot_outcome
#'
#' Helper function to generate a scatter plot of summary statistics.
#'
#' Provides customizable visualization of the distribution of quality control
#' metrics of cells. Used by `detect_damage`, `select_penalty`, and
#' `simulate_counts` to visualise the output of damage simulation.
#'
#' In all instances of use in this package, `plot_outcome` plots the y axis as
#' the proportion of counts originating from the mitochondrial genome. This
#' provides a rough estimate of damage where we anticipate that most cells
#' cluster near a common mitochondrial proportion and fewer, likely damaged,
#' outliers will lie at higher points. The x axis is either the number of
#' non-zero genes expressed or the proportion of counts coming from ribosomal
#' genes.
#'
#' Based on observations of true single-cell data, we anticipate these plots
#' to follow the same shape—an exponential decay—where cells with high
#' mitochondrial proportions have fewer expressed features and low
#' ribosomal proportions. Any egregious deviation from this shape indicates
#' that the simulated data poorly represents true data, especially when
#' compared directly to unaltered cells.
#'
#' @param data Data frame containing summary statistics
#' @param x Column name of the input data frame containing the values to be
#'  plotted on the x axis of the scatter plot.
#' @param y Column name of the input data frame containing the values to be
#'  plotted on the y axis of the scatter plot.
#' @param altered Boolean specifying whether the plot will be focusing on the
#'  data before or after alteration.
#'
#' * Default is FALSE.
#' @param mito_ribo Boolean specifying whether mitochondrial and ribosomal
#'  proportions are being plotted.
#' @param damage_column Column name of the input data frame containing the
#'   values to use for colouring each point in the plot.
#' @param target_damage Numeric vector specifying the range of target damage
#'   levels that
#' @param palette String specifying the three colours that will be used to
#'  create the continuous colour palette for colouring the 'damage_column'.
#'
#'  * Default is a range from purple to red,
#'  c("grey", "#7023FD", "#E60006").
#' @return 'ggplot2' object
#' @import ggplot2
#' @import patchwork
#' @import scales
#' @export
#' @examples
#' set.seed(42)  # For reproducibility
#' library(ggplot2)
#'
#' # Generate example an example data frame
#' df <- data.frame(
#'  Features = runif(1000, 0, 6000),
#'  Mt.percent = runif(1000, 0, 1),
#'  Damaged_Level = runif(1000, 0, 1)
#' )
#'
#' plot <- plot_outcome(
#'  data = df,
#'  x = "Features",
#'  y = "Mt.percent",
#'  damage_column = "Damaged_Level"
#')
plot_outcome <- function(
    data,
    x,
    y,
    altered = FALSE,
    mito_ribo = FALSE,
    damage_column = "Damaged_Level",
    target_damage =  c(0.1, 0.8),
    palette = c("grey", "#7023FD", "#E60006")
) {

  # Altered counts are plotted with damaged cells coloured
  if (altered) {
    p <- ggplot2::ggplot(data, aes(x = .data[[x]], y = .data[[y]], colour = .data[[damage_column]])) +
      ggplot2::scale_color_gradientn(
        colours = palette,
        values = scales::rescale(c(0, target_damage[[1]], 1)),
        limits = c(0, 1)
      ) +
      ggplot2::scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
      ggplot2::geom_point(size = 0.2) +
      ggplot2::theme_classic() +
      ggplot2::theme(
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.position = "bottom",
        legend.justification = "center",
        legend.title = element_text(face = "bold", hjust = 0.5)
      ) +
      ggplot2::guides(color = ggplot2::guide_colorbar(title = "Damage Level", title.position = "top"))

  } else {
    # Unaltered counts are plotted without damaged cells coloured
    p <- ggplot2::ggplot(data, aes(x = .data[[x]], y = .data[[y]])) +
      ggplot2::scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
      ggplot2::geom_point(size = 0.2, color = "gray") +
      ggplot2::theme_classic() +
      ggplot2::theme(
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.position = "none"
      )
  }

  if (mito_ribo){
    suppressMessages(suppressWarnings({
      p <- p + ggplot2::scale_x_continuous(limits = c(0, 1), labels = scales::number_format(accuracy = 0.1)) +
      ggplot2::scale_y_continuous(limits = c(0, 1), labels = scales::number_format(accuracy = 0.1)) +
      ggplot2::labs(x = "Ribosomal proportion", y = "Mitochondrial proportion")
    }))
    } else {
      suppressMessages(suppressWarnings({
        p <- p + ggplot2::labs(x = "Features expressed", y = "Mitochondrial proportion") +
          ggplot2::scale_y_continuous(limits = c(0, 1), labels = scales::number_format(accuracy = 0.1))
        }))
    }

  return(p)

}


