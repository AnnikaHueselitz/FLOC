#'Plot detection
#'
#'plots data, the signal and the true changepoint,
#'endpoint is the time of detection
#'
#'@param data vector of data points
#'@param tau_0 true changepoint; integer
#'@param alpha_minus intercept at tau_0 of the first segment; numeric
#'@param beta_minus slope on the first segment; numeric
#'@param alpha_plus intercept at tau_0 of the second segment; numeric
#'@param beta_plus slope on the secodn segment; numeric
#'
#'@return plot
#'
#'@export
#'
#'@import ggplot2
plot_detection <- function(data,
                           tau_0,
                           alpha_minus, beta_minus,
                           alpha_plus, beta_plus) {
  n <- length(data)
  time <- 1:n
  signal <- seg_lin_fun(1:n, tau_0, alpha_minus, beta_minus, alpha_plus, beta_plus)
  plot_data <- data.frame(data, signal, time)

  ggplot(plot_data) +

    ylab(NULL) +
    xlab(NULL) +

    #plot data
    geom_line(aes(x = time, y = data), alpha = 0.5) +

    #plot signal
    geom_line(aes(x = time, y = signal),color = "black") +

    #plot true changepoint
    geom_vline(xintercept = tau_0, color = "red", linetype = "dashed")
}
