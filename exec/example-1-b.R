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
#         theme(legend.position = "none", text = element_text(size = 16))
# )

set.seed(8529573L)
data_offline <- generate_data(T = 2, amplitude = 0, sensor_sd = 1.7, as_df = TRUE)

plot_trajectories(data = data_offline)


params <- list(sd_y = 4.0, A = -10.0,
               fun = fun_process,
               sd_x = 1.1, N = 1000,
               x_init = 0, sd_x_init = 0.5)


out <- particle_filter(data = data_offline, params,
                       resample_particles = FALSE,
                       rs_thresh = 0.5)


plot_filtering_estimates(out, data = data_offline)

print(out$logliksum)
print(out$ll)