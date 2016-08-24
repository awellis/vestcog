library(mvtnorm)
library(vestcog)
library(ggplot2)

eye <- function(n) {
    diag(nrow = n, ncol = n)
}

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
plot_trajectories <- function(data) {
    # set.seed(44234)
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

plot_filtering_estimates <- function(object, data, predict = FALSE) {

    color_palette <- viridis::plasma(n = 9)
    ggplot2::theme_set(
        theme_classic() +
            ggplot2::theme(
                axis.line.x = element_line(
                    colour = 'black',
                    size = 0.5,
                    linetype = 'solid'
                ),
                axis.line.y = element_line(
                    colour = 'black',
                    size = 0.5,
                    linetype = 'solid'
                )
            ) +
            theme(legend.position = "none", text = element_text(size = 24))
    )

    df <- with(object, {
        dplyr::data_frame(t = seq(1, Time),
                          vel_mean = vel_means,
                          # median = vel_medians,
                          vel_lower = vel_quantiles[1, ],
                          vel_upper = vel_quantiles[2, ],
                          acc_mean = acc_means,
                          acc_lower = acc_quantiles[1, ],
                          acc_upper = acc_quantiles[2, ],
                          x_true = data$velocity,
                          observations = data$observations,
                          loglik = loglik)
    })

    p <- ggplot2::ggplot(data = df, aes(x = t)) +
        ggplot2::geom_hline(yintercept = 0, linetype = "solid", alpha = 0.4) +
        ggplot2::geom_ribbon(aes(ymin = vel_lower, ymax = vel_upper), alpha = 0.2,
                             fill = "steelblue") +
        ggplot2::geom_line(aes(y = x_true), colour = "orangered1", alpha = 1,
                           linetype = "dashed", size = 2) +
        # geom_line(aes(y = mean), colour = color_palette[6], size = 1.4) +


        ggplot2::geom_line(aes(y = observations), colour = "darkgrey", size = 1.5,
                           linetype = "dotted") +

        ggplot2::geom_point(aes(y = observations), alpha = 1, fill = "white", colour = "white", shape = 21, size = 6) +
        ggplot2::geom_point(aes(y = observations), alpha = 1, fill = "white", colour = "grey40", shape = 21, size = 4) +
        # ggplot2::geom_line(aes(y = acc_mean), colour = "grey40",
        #                    linetype = "solid", size = 2, alpha = 1) +

        ggplot2::geom_line(aes(y = vel_mean), colour = "steelblue",
                           linetype = "solid", size = 2, alpha = 1) +

        ggplot2::scale_x_continuous(limits = c(1, 20), breaks = c(0, 5, 10, 15, 20),
                                    labels = c("0", "0.5", "1", "1.5", "2")) +
        ggplot2::ylab(expression(paste("Latent state: ", omega))) +

        xlab("Time [sec]")

    if (predict) {
        p <- p +
            ggplot2::geom_vline(xintercept = 10, linetype = "dashed",
                                color = color_palette[4], size = 1.5, alpha = 0.6) +
            ggplot2::geom_ribbon(data = filter(df, t > 9),
                                 aes(ymin = vel_lower, ymax = vel_upper), alpha = 0.4,
                                 fill = color_palette[4])
    }
    print(p)
}

# generate_data <- function(T = 2, dt = 0.1, amplitude = 20, sensor_sd = 1.7) {
#     nsteps <- T/dt
#     t <- seq(from = 0, to = T, length.out = nsteps)
#
#     # the following generates a motion profile with single-cycle sinusoidal
#     # acceleration
#     time <- t
#     position <- amplitude*T/(2*pi) * (t-(T/(2*pi)) * sin(2*pi*t/T))
#     velocity <- amplitude * T/(2 * pi) * (1-cos(2 * pi * t/T))
#     acceleration <- amplitude * sin(2 * pi * t/T)
#     trajectory <- rbind(position, velocity, acceleration)
#
#     observations <- rnorm(ncol(trajectory), trajectory[2,], sensor_sd)
#     out <- rbind(time, trajectory, observations)
# }


move <- function(A, t, Time) {
    # A * sin(2*pi*(t-1)/Time)
    A * sin(2 * pi * (t - 1)/Time)
}

# f <- function(x, t, Time, A, sigma, N, dt = 1) {
#     mu <- matrix(rep(0, 2*N), nrow = N)
#     mu[, 1] = x[, t-1, 1] + dt * 0.5 * move(A, t, Time)
#     mu[, 2] = x[, t-1, 2] + dt * move(A, t, Time)
#     # TODO: more elegant way of doing this?
#     t(apply(mu, 1, function(x) {rmvnorm(1, mean = x, sigma = sigma)}))
#     }
#
transfun <- function(x, t, Time, A, sigma, N, dt = 0.1, active = TRUE) {
    mu <- matrix(rep(0, 2*N), nrow = N)
    if (active) {
        mu[, 1] = x[, t-1, 1] + dt * move(A, t, Time)
    } else {
        mu[, 1] = x[, t-1, 1] +  x[, t-1, 2]
    }

    mu[, 2] = x[, t-1, 2]
    # TODO: more elegant way of doing this?
    t(apply(mu, 1, function(x) {rmvnorm(1, mean = x, sigma = sigma)}))
}

# transfun_passive <- function(x, t, Time, A, sigma, N, dt = 0.1) {
#     mu <- matrix(rep(0, 2*N), nrow = N)
#     mu[, 1] = x[, t-1, 1] + dt * x[, t-1, 2]
#     mu[, 2] = x[, t-1, 2]
#     # TODO: more elegant way of doing this?
#     t(apply(mu, 1, function(x) {rmvnorm(1, mean = x, sigma = sigma)}))
# }


run_particle_filter <- function(data, N, x_init = c(0, 0), sdx_init = 0.5*eye(2),
                                params, resample = TRUE, rs_thresh = 0.5,
                                active = TRUE) {


    sdx <- params$sdx
    sdy <- params$sdy
    transfun <- params$transfun
    A <- params$A
    Time = length(data$observations)

    x <- array(rep(0, N*Time), dim = c(N, Time, 2))
    # dimnames(x) <- c("particles", "time", "statevar")
    weights <- matrix(rep(0, N*Time), nrow = N, ncol = Time) #matrix(nrow =  N, ncol = Time)
    loglik <- rep(0, Time)

    x[, 1, ] <- mvtnorm::rmvnorm(N, x_init, sdx_init)

    weights[, 1] <- dnorm(data$observations[1], x[, 1, 1], sdy)
    loglik[1] <- log(sum(weights[, 1]))
    weights[, 1] <- weights[, 1]/sum(weights[, 1])

    idx <- sample(dim(x)[1], replace = TRUE, size = N, prob = weights[, 1])
    x[, 1, ] <- x[idx, 1, ]

    # x[, 1, ] <- sample(x[, 1, ], replace = TRUE, size = N, prob = weights[, 1])


    for (t in seq(2, Time)) {
        x[, t, ] <- transfun(x, t, Time, A, sigma = sdx, N, active = active)

        if (!is.na(data$observations[t])) {
            weights[, t] <- dnorm(data$observations[t], x[, t, 1], sdy)
        } else {
            weights[, t] <- 1/N
        }
        loglik[t] <- log(sum(weights[, t]))

        if (resample) {
            if (neff(weights[, t]) < rs_thresh * N) {
                # x[, t, ] <- sample(x[, 1, 2], replace = TRUE, size = N, prob = weights[, t])
                idx <- sample(dim(x)[1], replace = TRUE, size = N, prob = weights[, t])
                x[, t, ] <- x[idx, t, ]
            }
        }
    }

    vel <- x[, , 1]
    acc <- x[, , 2]
    out <- list(x = x, weights = weights, loglik = loglik,
                vel = vel, acc = acc,
                acc_means = apply(acc, 2, mean),
                vel_means = apply(vel, 2, mean),
                acc_quantiles = apply(acc, 2,
                                      function(x) quantile(x,
                                                           probs = c(0.025,
                                                                     0.975))),
                vel_quantiles = apply(vel, 2,
                                      function(x) quantile(x,
                                                           probs = c(0.025,
                                                                     0.975))),
                Time = Time,
                N = N)
}


data <- vestcog::generate_data(T = 2, amplitude = 20,
                               sensor_sd = 1.7, as_df = TRUE,
                               seed = TRUE)
# data$observations[10:length(data$observations)] <- NA

# plot_trajectories(data = data)

params <- list(sdx = c(0.1, 1.2) * eye(2), sdy = 3.8, A = 20, transfun = transfun)

out <- run_particle_filter(data = data, N = 3000,
                           x_init = c(0, 0), sdx_init = 0.1 * eye(2),
                           params, resample = TRUE, rs_thresh = 0.5,
                           active = TRUE)

plot_filtering_estimates(out, data = data, predict = FALSE)
