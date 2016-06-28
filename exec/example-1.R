library(ggplot2)

ggplot2::theme_set(theme_classic() +
                       ggplot2::theme(axis.line.x = element_line(colour = 'black', size = 0.5, linetype = 'solid'),
                                      axis.line.y = element_line(colour = 'black', size = 0.5, linetype = 'solid')) +
                       theme(legend.position = "none", text = element_text(size = 16)))


data <- generate_data(T = 2, amplitude = 20, sensor_sd = 1.7, as_df = TRUE)


f <- function(x, t, Time, A, sd, N) {
    rnorm(n = N, mean =  x[, t-1] + A * sin(2*pi*(t-1)/Time), sd = sd)
}

params <- list(sd_x = 1.9, sd_y = 4.0, A = 2.0, fun_x = f)

out <- particle_filter(data = data, N = 100, Time = length(data$observations),
                       x_init = 0, sdx_init = 0.5,
                       params, resample = TRUE, rs_thresh = 0.5)
plot_filtering_estimates(out, data = data)
