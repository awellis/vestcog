
# Rbiips mode 3: ----
library(ggplot2)
library(dplyr)
library(tidyr)
library(Rbiips)

library(viridis)
color_palette <- viridis::plasma(n = 9)

sd2precision <- function(sd) {
    prec <- 1/(sd^2)
    prec
}

neff <- function(weights) {
    1/sum(weights^2)
}

Mode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
}

generate_data <- function(T = 2, dt = 0.1, amplitude = 20, sensor_sd = 1.7) {
    nsteps <- T/dt
    t <- seq(from = 0, to = T, length.out = nsteps)

    # the following generates a motion profile with single-cycle sinusoidal
    # acceleration
    time <- t
    position <- amplitude*T/(2*pi) * (t-(T/(2*pi)) * sin(2*pi*t/T))
    velocity <- amplitude * T/(2 * pi) * (1-cos(2 * pi * t/T))
    acceleration <- amplitude * sin(2 * pi * t/T)
    trajectory <- rbind(position, velocity, acceleration)

    observations <- rnorm(ncol(trajectory), trajectory[2,], sensor_sd)
    out <- rbind(time, trajectory, observations)
}

plot_trajectories <- function(data) {
    set.seed(44234)
    data <- data %>% gather(key = "key", value = "value", -time)

    g <- ggplot(data = data, aes(x = time, y = value, color = key)) +
        geom_line(data = filter(data,
                                key %in% c("acceleration",
                                           "velocity",
                                           "position")),
                  size = 2) +
        geom_point(data = filter(data, key == "observations"),
                   size = 3, shape = 15, color = "black") +
        xlab("Time") + ylab("Angular velocity [deg]") +
        # labs(title = "Natural head motion") +
        # scale_colour_brewer(palette = "Set1")
        scale_colour_manual(values = sample(color_palette)) +
        theme(legend.title = element_blank())
    print(g)
}

set.seed(8529573L)
trajectory_data <- generate_data(T = 2, amplitude = 20, sensor_sd = 3.0) %>%
    t() %>%
    data.frame()

plot_trajectories(data = trajectory_data)


## plot filtering estimates ----
plot_filtering_estimates <- function(df) {
    p <- ggplot(data = df, aes(x = t)) +
        geom_ribbon(aes(ymin = pos_lower, ymax = pos_upper), alpha = 0.3,
                    fill = color_palette[6]) +
        geom_ribbon(aes(ymin = vel_lower, ymax = vel_upper), alpha = 0.3,
                    fill = color_palette[2]) +
        geom_line(aes(y = velocity), colour = color_palette[7], alpha = 0.9,
                  linetype = 2, size = 1.2) +

        geom_line(aes(y = vel_mean), colour = color_palette[2], size = 1.6) +
        geom_line(aes(y = pos_mean), colour = color_palette[6], size = 1.6) +

        geom_point(aes(y = observations), colour = color_palette[1],
                   size = 2, shape = 15, alpha = 0.6) +
        theme_bw() +
        theme(legend.title = element_blank()) +
        ylab("") + xlab("")
    # ggtitle("Velocity storage")

    # geom_vline(xintercept=as.numeric(as.Date("1959-12-01")), linetype=2) +
    # ggtitle(paste0("ARIMA -- Holdout MAPE = ", round(100*MAPE,2), "%")) +
    # theme(axis.text.x=element_text(angle = -90, hjust = 0))

    print(p)
}



dt <- 0.1

t_max <- length(trajectory_data$observations)

# starting position
mean_x_init <- c(0, 0)

prec_x_init <- diag(c(
    sd2precision(1),
    sd2precision(1)), nrow = 2)

prec_x <- diag(c(
    sd2precision(0.5),
    sd2precision(1)), nrow = 2)

prec_y <- sd2precision(3.7)


# A: process model: implements Newton's laws of motion:
# x[1, t] = 1 * x[1, t-1] + dt * x[2, t-1]
# x[2, t] = 0 * x[1, t-1] + 1 * x[2, t-1]

A <- matrix(c(1, dt,
              0, 1), nrow = 2, byrow = TRUE)

# C: measurement model: y measures velocity
# y[t] = 0 * x[1, t] + 1 * x[2, t]
C <- matrix(c(0, 1), nrow = 1, byrow = TRUE)

# B: control input
B <- matrix(c(dt^2/2, dt), nrow = 2, byrow = TRUE)

## specify data ----
data = list(t_max = t_max,
            y = trajectory_data$observations,
            a = trajectory_data$acceleration/20,
            mean_x_init = mean_x_init,
            prec_x_init = prec_x_init,
            prec_x = prec_x,
            prec_y = prec_y,
            A = A, B = B, C = C)



## Rbiips model -----

model_string <- "
var x[2, t_max], y[t_max]

model {
x[, 1] ~ dmnorm(mean_x_init, prec_x_init)
y[1] ~ dnorm(C %*% x[, 1], prec_y)
amplitude <- 20

for (t in 2:t_max) {
x[, t] ~ dmnorm(A %*% x[, t-1] + amplitude * B %*% a[t-1], prec_x)
y[t] ~ dnorm(C %*% x[, t], prec_y)
}
}"

model = biips_model(file = textConnection(model_string), data = data, sample_data = FALSE)
data = model$data()

n_particles <- 1000 # Number of particles

variables <- c('x') # Variables to be monitored
mn_type <- 'f'
rs_type <- c('systematic', 'multinomial', 'stratified')
rs_thres <- 0.5 # Optional parameters

# Run SMC ----
out_smc <- biips_smc_samples(model, variables, n_particles,
                             type = mn_type, rs_type = rs_type[2],
                             rs_thres = rs_thres)

diag_smc <- biips_diagnosis(out_smc)

summ_smc <- biips_summary(out_smc, probs=c(.025, .975))


# #### Plot Filtering estimates ----
pos_f_mean <- summ_smc$x$f$mean[1, ]
pos_f_lower <- summ_smc$x$f$quant$`0.025`[1, ]
pos_f_upper <- summ_smc$x$f$quant$`0.975`[1, ]

vel_f_mean <- summ_smc$x$f$mean[2, ]
vel_f_lower <- summ_smc$x$f$quant$`0.025`[2, ]
vel_f_upper <- summ_smc$x$f$quant$`0.975`[2, ]

df <- data.frame(t = 1:t_max,
                 pos_mean = pos_f_mean,
                 pos_lower = pos_f_lower,
                 pos_upper = pos_f_upper,
                 vel_mean = vel_f_mean,
                 vel_lower = vel_f_lower,
                 vel_upper = vel_f_upper,
                 velocity = trajectory_data$velocity,
                 observations = trajectory_data$observations)

plot_filtering_estimates(df)

