
mplot_func <- function(mplot, ax, plot_title, title_more, insamp, colo) {
  # if signal forecasts are computed then revision matrix is longer than data 
  # matrix
  ymin <- min(mplot, na.rm = TRUE)
  ymax <- max(mplot, na.rm = TRUE)
  # color vector for the lines and the text
  if(length(colo) == 0) {
    colo <- rainbow(ncol(mplot))
  }
  # plot the first data column
  plot(mplot[, 1], ylim = range(mplot, na.rm = TRUE), type = "l", axes = FALSE,
       col = colo[1], main = plot_title, xlab = "", ylab = "")
  abline(h = 0)
  # plot additional titles
  if(length(title_more) > 0) {
    mtext(title_more[1], col = colo[1])
  }
  # plot the remaining data columns
  if(ncol(mplot) > 1) {
    for(j in 2 : ncol(mplot)) {
      lines(mplot[, j], col = colo[j])
      # plot the remaining additional titles
      if(length(title_more) > 0) {
        mtext(title_more[j], col = colo[j], line = -(j - 1))
      }
      # plot a vertical line at the end of the in-sample window
      if((insamp[1] < 1.e+90) & (length(insamp) == ncol(mplot) - 1)) {
        abline(v = insamp[j - 1], col = colo[j])
        text(x = insamp[j - 1], y = -2 - (j - 1) * 0.3, col = colo[j],
             "End of In-Sample Window")
      }
    }
  }
  # add the x- and y-axis and draw a box around the current plot
  axis(1, at = 1 : nrow(mplot), labels = ax)
  axis(2)
  box()
}
