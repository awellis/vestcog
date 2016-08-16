library(ggplot2)

# ggplot2::theme_set(
#     theme_classic() +
#         ggplot2::theme(
#             axis.line.x = element_line(
#                 colour = 'black',
#                 size = 0.5,
#                 linetype = 'solid'
#             ),
#             axis.line.y = element_line(
#                 colour = 'black',
#                 size = 0.5,
#                 linetype = 'solid'
#             )
#         ) +
#         theme(legend.position = "none", text = element_text(size = 28))
# )


motion_data <- generate_data(T = 2, amplitude = 20, sensor_sd = 1.7, as_df = TRUE)
plot_trajectories(motion_data, facet = FALSE)

fun_control <- function(A, t, Time, N) {
    rnorm(n = N,
          mean = A * sin(2 * pi * (t - 1) / Time),
          sd = 0.1)
}

fun_process <- function(x, t, Time, A, sd, N) {
    rnorm(
        n = N,
        mean =  x[, t - 1] + fun_control(A, t, Time, N),
        sd = sd
    )
}

# fun_process <- function(x, t, Time, A, sd, N) {
#     rnorm(n = N, mean =  x[, t-1] + A*sin(2*pi*(t-1)/Time), sd = sd)
# }

# f_0 <- function(x, t, Time, A, sd, N) {
#     rnorm(n = N, mean =  x[, t-1], sd = sd)
# }
#
# f_2 <- function(x, t, Time, A, sd_active, sd_passive, N) {
#     passive <- rnorm(n = N, mean =  passive[, t-1], sd = sd_passive)
#     active <- rnorm(n = N, mean =  active[, t-1], sd = sd_active)
#     velocity <- passive + passive
# }


# f2 <- function(x, t, Time, A, sd, N) {
#     # A * sin(2*pi*(t-1)/Time)
#     rnorm(n = N, mean =  A * sin(2*pi*(t-1)/Time), sd = sd)
# }


params <- list(sd_y = 4.0, A = 2.0,
               fun = fun_process,
               sd_x = 1.1, N = 1000,
               x_init = 0, sd_x_init = 0.5)

out <- particle_filter(data = motion_data,
                       params = params,
                       resample_particles = FALSE,
                       rs_thresh = 0.5)

plot_filtering_estimates(out, data = motion_data)

print(out$logliksum)
# print(out$ll)

# out$loglik %>% plot()
