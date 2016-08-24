# Rbiips mode 4: simple model with control ----
library(ggplot2)
library(dplyr)
library(tidyr)
library(vestcog)
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


trdata <- vestcog::generate_data(T = 2, amplitude = 20, sensor_sd = 1.7,
                                 as_df = TRUE, seed = TRUE)

dt <- 0.1
t_max <- length(trdata$observations)
mean_w_init <- 0
prec_w_init <- sd2precision(1)
mean_p_init <- 0
prec_p_init <- sd2precision(1)

prec_w <- sd2precision(0.5)
prec_p <- sd2precision(1.0)
prec_y <- sd2precision(1.7)


data = list(t_max = t_max,
            dt = dt,
            y = trdata$observations,
            #a = trajectory_data$acceleration/20,
            A = 20,
            mean_w_init = mean_w_init,
            prec_w_init = prec_w_init,
            mean_p_init = mean_p_init,
            prec_p_init = prec_p_init,
            prec_w = prec_w,
            prec_p = prec_p,
            prec_y = prec_y,
            two_pi = 2*pi)

model_string <- "
model {
w[1] ~ dnorm(mean_w_init, prec_w_init)
p[1] ~ dnorm(mean_p_init, prec_p_init)
y[1] ~ dnorm(w[1], prec_y)

for (t in 2:t_max) {
# p: position
# w: velocity
# a: control input
a[t-1] <- A * sin(two_pi*(t-1)/t_max)
p[t] ~ dnorm(p[t-1] + dt * w[t-1] + dt^2/2 * a[t-1], prec_p)
w[t] ~ dnorm(w[t-1] + dt * a[t-1], prec_w)

y[t] ~ dnorm(w[t], prec_y)
}
}"


model = biips_model(file = textConnection(model_string), data = data,
                    sample_data = FALSE)
n_particles <- 1000 # Number of particles

variables <- c('w', 'p') # Variables to be monitored
mn_type <- 'f'
rs_type <- c('systematic', 'multinomial', 'stratified')
rs_thres <- 0.5 # n_particles # Optional parameters

# Run SMC ----
out_smc <- biips_smc_samples(model, variables, n_particles,
                             type = mn_type, rs_type = rs_type[2],
                             rs_thres = rs_thres)
summ_smc <- biips_summary(out_smc, probs=c(.025, .975))

vel_f_mean <- summ_smc$w$f$mean
vel_f_lower <- summ_smc$w$f$quant$`0.025`
vel_f_upper <- summ_smc$w$f$quant$`0.975`

pos_f_mean <- summ_smc$p$f$mean
pos_f_lower <- summ_smc$p$f$quant$`0.025`
pos_f_upper <- summ_smc$p$f$quant$`0.975`

df <- data.frame(t = 1:t_max,
                 vel_mean = vel_f_mean,
                 vel_lower = vel_f_lower,
                 vel_upper = vel_f_upper,
                 pos_mean = pos_f_mean,
                 pos_lower = pos_f_lower,
                 pos_upper = pos_f_upper,
                 velocity = trdata$velocity,
                 observations = trdata$observations)

plot_filtering_estimates <- function(df) {
    p <- ggplot(data = df, aes(x = t)) +

        geom_ribbon(aes(ymin = vel_lower, ymax = vel_upper), alpha = 0.3,
                    fill = color_palette[2]) +
        geom_line(aes(y = velocity), colour = color_palette[7], alpha = 0.9,
                  linetype = 2, size = 1.2) +

        geom_line(aes(y = vel_mean), colour = color_palette[2], size = 1.6) +

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


plot_filtering_estimates <- function(df, pos = TRUE) {
    p <- ggplot(data = df, aes(x = t)) +
        geom_ribbon(aes(ymin = vel_lower, ymax = vel_upper), alpha = 0.3,
                    fill = color_palette[2]) +
        geom_line(aes(y = velocity), colour = color_palette[7], alpha = 0.9,
                  linetype = 2, size = 1.2) +

        geom_line(aes(y = vel_mean), colour = color_palette[2], size = 1.6) +

        geom_point(aes(y = observations), colour = color_palette[1],
                   size = 2, shape = 15, alpha = 0.6) +
        theme_bw() +
        theme(legend.title = element_blank()) +
        ylab("") + xlab("")
    # ggtitle("Velocity storage")

    # geom_vline(xintercept=as.numeric(as.Date("1959-12-01")), linetype=2) +
    # ggtitle(paste0("ARIMA -- Holdout MAPE = ", round(100*MAPE,2), "%")) +
    # theme(axis.text.x=element_text(angle = -90, hjust = 0))

    if (pos) {
        p <- p + geom_ribbon(aes(ymin = pos_lower, ymax = pos_upper), alpha = 0.3,
                             fill = color_palette[6]) +
            geom_line(aes(y = pos_mean), colour = color_palette[6], size = 1.6)
    }

    print(p)
}


plot_filtering_estimates(df, pos = TRUE)
