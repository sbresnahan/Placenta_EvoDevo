ex_plot_theme <- function() {
  ggplot2::theme_bw(base_family = "Liberation Sans") +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}
ex_plot_forest <- function(d, path, mediation = TRUE, top = 20L) {
  if (mediation) {est <- "indirect_beta"; lo <- "indirect_low"; hi <- "indirect_high"} else {
    est <- "slope"; lo <- "slope_low"; hi <- "slope_high"
    if (!lo %in% names(d)) d[, `:=`(slope_low = slope - 1.96 * slope_se, slope_high = slope + 1.96 * slope_se)]
  }
  d <- data.table::copy(d[status == "eligible" & is.finite(get(est))])
  if (!nrow(d)) return(invisible(NULL))
  d <- head(d[order(primary_q, primary_p, na.last = TRUE)], top)
  d[, label := paste(feature_id, family_id, exposure)]
  p <- ggplot2::ggplot(d, ggplot2::aes(y = stats::reorder(label, get(est)), x = get(est))) +
    ggplot2::geom_vline(xintercept = 0, color = "grey65") +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = get(lo), xmax = get(hi)), orientation = "y", width = .15) +
    ggplot2::geom_point(color = "#0072B2") + ggplot2::labs(x = "POE beta change", y = NULL) + ex_plot_theme()
  ggplot2::ggsave(path, p, width = 9, height = max(3, .28 * nrow(d) + 1.2), dpi = 200)
  invisible(p)
}
ex_plot_beta_curve <- function(result, path, M = 0) {
  if (is.null(result$fit)) return(invisible(NULL))
  x <- seq(min(result$samples$E), max(result$samples$E), length.out = 100)
  d <- data.table::rbindlist(lapply(x, function(e) {
    v <- ex_contrast(result, ex_weights(result, e, M), result$dosage_null)
    data.table::data.table(exposure = e, beta = v["beta"], low = v["low"], high = v["high"])
  }))
  p <- ggplot2::ggplot(d, ggplot2::aes(exposure, beta)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = low, ymax = high), fill = "#56B4E9", alpha = .25) +
    ggplot2::geom_line(color = "#0072B2") + ggplot2::geom_hline(yintercept = 0, color = "grey65") +
    ggplot2::labs(x = "Exposure", y = "POE beta") + ex_plot_theme()
  ggplot2::ggsave(path, p, width = 6, height = 4, dpi = 200); invisible(p)
}
ex_plot_signed_magnitude <- function(d, path) {
  if (!nrow(d)) return(invisible(NULL))
  p <- ggplot2::ggplot(d[status == "eligible"], ggplot2::aes(delta_signed, delta_magnitude)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey75") +
    ggplot2::geom_vline(xintercept = 0, color = "grey75") +
    ggplot2::geom_point(color = "#0072B2", alpha = .6) +
    ggplot2::labs(x = "Signed beta change", y = "Absolute beta change") + ex_plot_theme()
  ggplot2::ggsave(path, p, width = 5, height = 4, dpi = 200); invisible(p)
}
ex_plot_joint_support <- function(samples, path) {
  if (!nrow(samples) || !all(c("E", "M") %in% names(samples))) return(invisible(NULL))
  p <- ggplot2::ggplot(samples, ggplot2::aes(E, M)) + ggplot2::geom_point(color = "#0072B2", alpha = .65) +
    ggplot2::labs(x = "Maternal PFAS z-score", y = "TE abundance z-score") + ex_plot_theme()
  ggplot2::ggsave(path, p, width = 5, height = 4, dpi = 200); invisible(p)
}
ex_plot_pfas_correlations <- function(metadata, path) {
  ex_need(metadata, EX_PFAS)
  m <- stats::cor(as.matrix(metadata[, ..EX_PFAS]), use = "pairwise.complete.obs")
  d <- data.table::as.data.table(as.table(m)); data.table::setnames(d, c("PFAS1", "PFAS2", "correlation"))
  p <- ggplot2::ggplot(d, ggplot2::aes(PFAS1, PFAS2, fill = correlation)) +
    ggplot2::geom_tile() + ggplot2::scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00",
      limits = c(-1, 1)) + ggplot2::labs(x = NULL, y = NULL, fill = "Correlation") +
    ex_plot_theme() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  ggplot2::ggsave(path, p, width = 6, height = 5, dpi = 200); invisible(p)
}
