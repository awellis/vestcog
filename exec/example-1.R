library(ggplot2)

motiondata <-  generate_data(T = 2, amplitude = 20, sensor_sd = 2.0, as_df = TRUE)
plot_trajectories(motiondata)

f_0 <- function(x, t, Time, A, sd, N) {
    rnorm(n = N, mean =  x[, t-1], sd = sd)
}

f_control <- function(x, t, Time, A, sd, N) {
    rnorm(n = N, mean =  x[, t-1] + dt*A*sin(2*pi*(t-1)/Time), sd = sd)
}

dt <- 0.1

params <- list(sd_y = 2.0, A = 20.0,
               fun = f_0,
               sd_x = 0.2, N = 3000,
               x_init = 0, sd_x_init = 0.5)

out <- particle_filter(data = motiondata,
                       params = params,
                       resample_particles = TRUE,
                       rs_thresh = 0.9)

p <- plot_filtering_estimates(out, data = motiondata, predict = FALSE)

print(p)

print(out$logliksum)
# print(out$ll)

# out$loglik %>% plot()
