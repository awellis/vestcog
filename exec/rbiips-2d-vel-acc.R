## Rbiips example with 2-d state space:
## velocity and acceleration

library(Rbiips)
library(ggplot2)
# library(dplyr)

# plot settings -----------------------------------------------------------
cb_palette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
                "#0072B2", "#D55E00", "#CC79A7", "grey80")
ggplot2::theme_set(theme_bw())

set.seed(8573L)


# generate trajectory and noisy measurements ----
dt <- 0.1
amplitude <- 20
T <- 2
nsteps <- T/dt
t <- seq(from = 0, to = T, length.out = nsteps)

# the following generates a motion profile with single-cycle sinusiodal
# acceleration

position <- amplitude*T/(2*pi) * (t-(T/(2*pi)) * sin(2*pi*t/T))
velocity <- amplitude * T/(2 * pi) * (1-cos(2 * pi * t/T))
acceleration <- amplitude * sin(2 * pi * t/T)
trajectory <- rbind(position, velocity, acceleration)

sensor_sd <- 1.7
observations <- rnorm(ncol(trajectory), trajectory[2,], sensor_sd)

# observations[10:15] <- NA


# converts standard deviation to precision ----
sd2precision <- function(sd) {
    prec <- 1/(sd^2)
    prec
}



## 2 states ----
## x[1,] is velocity
## x[2,] is acceleration
t_max <- length(observations)

# starting position
mean_x_init <- c(0, 0)

prec_x_init <- diag(c(
    sd2precision(1),
    sd2precision(1)), nrow = 2)

prec_x <- diag(c(
    sd2precision(1.5),
    sd2precision(1.5)), nrow = 2)

prec_y <- sd2precision(1.7)


# A: process model: implements Newton's laws of motion:
# x[1, t] = 1 * x[1, t-1] + dt * x[2, t-1]
# x[2, t] = 0 * x[1, t-1] + 1 * x[2, t-1]

A <- matrix(c(1, dt,
              0, 1), nrow = 2, byrow = TRUE)

# C: measurement model: y measures velocity
# y[t] = 1 * x[1, t] + 0 * x[2, t]
C <- matrix(c(1, 0), nrow = 1, byrow = TRUE)

## run in offline mode ----
## this doesn't work yet, we need to apply a process model, because there is nothin to learn from the data
# observations[2:length(observations)] <- NA
# prec_x <- diag(c(
#     sd2precision(0.1),
#     sd2precision(0.1)), nrow = 2)
# prec_y <- sd2precision(10.7)

## specify data ----
data = list(t_max = t_max,
            y = observations,
            mean_x_init = mean_x_init,
            prec_x_init = prec_x_init,
            prec_x = prec_x,
            prec_y = prec_y,
            A = A, C = C)




model_string <- "
var x[2, t_max], y[t_max]

model {
x[, 1] ~ dmnorm(mean_x_init, prec_x_init)
y[1] ~ dnorm(C %*% x[, 1], prec_y)

for (t in 2:t_max) {
x[, t] ~ dmnorm(A %*% x[, t-1], prec_x)
y[t] ~ dnorm(C %*% x[, t], prec_y)
}
}"


model = biips_model(file = textConnection(model_string), data = data, sample_data = FALSE)
data = model$data()

n_particles <- 1000 # Number of particles

variables <- c('x') # Variables to be monitored
mn_type <- 'f'
# rs_type <- c('systematic', 'multinomial', 'stratified')
rs_thres <- 0.5 # Optional parameters

# Run SMC ----
out_smc <- biips_smc_samples(model, variables, n_particles,
                             type = mn_type, rs_type = 'multinomial',
                             rs_thres = rs_thres)

diag_smc <- biips_diagnosis(out_smc)

summ_smc <- biips_summary(out_smc, probs=c(.025, .975))


# #### Plot Filtering estimates ----
vel_f_mean <- summ_smc$x$f$mean[1, ]
vel_f_lower <- summ_smc$x$f$quant$`0.025`[1, ]
vel_f_upper <- summ_smc$x$f$quant$`0.975`[1, ]

acc_f_mean <- summ_smc$x$f$mean[2, ]
acc_f_lower <- summ_smc$x$f$quant$`0.025`[2, ]
acc_f_upper <- summ_smc$x$f$quant$`0.975`[2, ]

df <- data.frame(t = 1:t_max,
                 vel_mean = vel_f_mean,
                 vel_lower = vel_f_lower,
                 vel_upper = vel_f_upper,
                 acc_mean = acc_f_mean,
                 acc_lower = acc_f_lower,
                 acc_upper = acc_f_upper,
                 velocity = trajectory[2,],
                 observations = observations)

## plot filtering estimates ----

plot_filtering_estimates <- function(df) {
    p <- ggplot(data = df, aes(x = t)) +
        geom_ribbon(aes(ymin = vel_lower, ymax = vel_upper), alpha = 0.3,
                    fill = cb_palette[6]) +
        # geom_ribbon(aes(ymin = acc_lower, ymax = acc_upper), alpha = 0.3,
        # fill = cb_palette[5]) +
        geom_line(aes(y = velocity), colour = cb_palette[7], alpha = 0.9,
                  linetype = 2, size = 1.2) +

        geom_line(aes(y = acc_mean), colour = cb_palette[5], size = 1.6) +
        geom_line(aes(y = vel_mean), colour = cb_palette[6], size = 1.6) +

        geom_point(aes(y = observations), colour = cb_palette[1],
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

plot_filtering_estimates(df)


## session info ----
# devtools::session_info()
