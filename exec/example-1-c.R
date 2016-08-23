
motiondata <- generate_data(T = 2, amplitude = 20, sensor_sd = 1.7, as_df = TRUE)

data_missing <- motiondata

# data_missing$observations[10:length(data_missing$observations)] <- NA
data_missing$observations[10:15] <- NA

f_0 <- function(x, t, Time, A, sd, N) {
    rnorm(n = N, mean =  x[, t-1], sd = sd)
}

f_control <- function(x, t, Time, A, sd, N) {
    rnorm(n = N, mean =  x[, t-1] + dt*A*sin(2*pi*(t-1)/Time), sd = sd)
}

dt <- 0.1

params <- list(sd_y = 1.8, A = 20.0,
               fun = f_control,
               sd_x = 0.2, N = 1000,
               x_init = 0, sd_x_init = 0.5)

out <- particle_filter(data = data_missing,
                       params = params,
                       resample_particles = TRUE,
                       rs_thresh = 0.5)

plot_filtering_estimates(out, data = data_missing, predict = TRUE)
