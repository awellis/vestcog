#' @import ggplot2
#' @import dplyr
#' @import tidyr
# ' @importFrom viridis plasma

#' @export
update <- function(prior, data, sd_y) {
    out <- dnorm(data, prior, sd_y, log = FALSE)
}

#' @export
normalize <- function(weights) {
    out <- weights/sum(weights)
}

#' @export
resample <- function(samples, weights) {
    out <- sample(samples, replace = TRUE, size = length(samples), prob = weights)
}


#' @export
particle_filter <- function(data, params, resample_particles = TRUE, rs_thresh = 0.5) {

    # unpack parameters
    sd_c <- params$sd_c
    # sd_x <- params$sd_x
    sd_x <- params$sd_x
    sd_y <- params$sd_y
    fun <- params$fun

    # fun_c <- params$fun_c
    A <- params$A
    N <- params$N
    Time = length(data$observations)
    x_init <- params$x_init
    sd_x_init <- params$sd_x_init

    # initialize variables
    x <- matrix(rep(0, N*Time), nrow = N, ncol = Time)  #matrix(nrow =  N, ncol = Time)
    weights <- matrix(rep(0, N*Time), nrow = N, ncol = Time) #matrix(nrow =  N, ncol = Time)
    n_eff <- rep(0, Time)
    resampled <- logical(length = Time)
    loglik <- rep(0, Time)
    ll <- 0

    # sample from prior distribution
    x[, 1] <- rnorm(N, x_init, sd_x_init)
    weights[, 1] <- update(data$observations[1], x[, 1], sd_y)
    loglik[1] <- log(sum(weights[, 1])) - log(N)
    ll <- ll + log(sum(weights[, 1])) - log(N)

    weights[, 1] <- normalize(weights[, 1])
    n_eff[1] <- neff(weights[, 1])

    # x[, 1] <- resample(x[, 1], weights[, 1])

    for (t in seq(2, Time)) {

        # predict 1 step ahead using process model as proposal distribution
        x[, t] <- fun(x, t, Time, A, sd = sd_x, N)

        if (!is.na(data$observations[t])) {
            weights[, t] <- update(data$observations[t], x[, t], sd_y)
        } else {
            weights[, t] <- 1/N
        }

        loglik[t] <- log(sum(weights[, t])) - log(N)
        ll <- ll + log(sum(weights[, t])) - log(N)

        weights[, t] <- normalize(weights[, t])
        n_eff[t] <- neff(weights[, t])

        if (resample_particles) {
            if (n_eff[t] < rs_thresh * N) {
                x[, t] <- resample(x[, t], weights[, t])
                resampled[t] <- TRUE
            } else resampled[t] <- FALSE
        } else resampled[t] <- FALSE
    }

    logliksum <- sum(loglik)
    # logliksum <- sum(loglik) - log(N)


    out <- list(x = x, weights = weights, loglik = loglik,
                logliksum = logliksum, ll = ll,
                ess = round(n_eff, digits = 0), resample = resampled,
                x_means = apply(x, 2, mean),
                x_medians = apply(x, 2, median),
                x_quantiles = apply(x, 2,
                                    function(x) quantile(x,
                                                         probs = c(0.025,
                                                                   0.975))),
                Time = Time,
                N = N)
}


#' @export
plot_filtering_estimates <- function(object, data, predict = FALSE,
                                     bw = TRUE) {

    # color_palette <- viridis::plasma(n = 9)
    color_palette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    if (bw) {
        true_color <- "grey40"
        estimate_color <- "grey60"
        predict_color <- "grey50"
    } else {
        true_color <- color_palette[2]
        estimate_color <- color_palette[6]
        predict_color <- color_palette[8]
    }


    if (object$Time == 20) {
        xend <- 20
        xbreaks <- c(0, 5, 10, 15, 20)
        xlabels = c("0", "0.5", "1", "1.5", "2")
    } else {
        xend <- 80
        xbreaks <- seq(from = 0, to = 80, by = 10)
        xlabels = as.character(seq(from = 0, to = 8, by = 1))
    }

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
                          xmean = x_means,
                          xmedian = x_medians,
                          xlower = x_quantiles[1, ],
                          xupper = x_quantiles[2, ],
                          x_true = data$velocity,
                          observations = data$observations,
                          loglik = loglik,
                          ess = ess,
                          resample = resample)
    })

    p <- ggplot2::ggplot(data = df, aes(x = t)) +
        ggplot2::geom_hline(yintercept = 0, linetype = "solid", alpha = 0.4) +
        ggplot2::geom_ribbon(aes(ymin = xlower, ymax = xupper), alpha = 0.1,
                             fill = estimate_color) +

        ggplot2::geom_line(aes(y = x_true), colour = true_color,
                           linetype = "dashed", size = 2) +
        # geom_line(aes(y = mean), colour = color_palette[6], size = 1.4) +


        # ggplot2::geom_line(aes(y = observations), colour = "darkgrey", size = 1.5,
        #                    linetype = "dotted") +

        ggplot2::geom_point(aes(y = observations), alpha = 1.0, fill = "white",
                            colour = "black", shape = 21, size = 6) +
        # ggplot2::geom_point(aes(y = observations), alpha = 0.8, fill = "white", colour = "grey40", shape = 21, size = 4) +

        ggplot2::geom_line(aes(y = xmean), colour = "black",
                           linetype = "solid", size = 2, alpha = 1) +

        ggplot2::scale_x_continuous(limits = c(1, xend),
                                    breaks = xbreaks,
                                    labels = xlabels) +
        # ggplot2::ylab(expression(paste("Latent state: ", omega))) +
        ggplot2::ylab(expression(paste("Angular velocity estimate [deg/s]"))) +
        ggplot2::xlab("Time [sec]")

    if (predict) {
        p <- p +
            ggplot2::geom_vline(xintercept = c(10, 15), linetype = "dashed",
                                color = predict_color, size = 1.5, alpha = 0.6) +
            ggplot2::geom_ribbon(data = dplyr::filter(df, t > 9, t < 16),
                                 aes(ymin = xlower, ymax = xupper), alpha = 0.4,
                                 fill = predict_color)
    }
    print(p)
}


# f <- function(x, t, Time, A, sd, N) {
#     rnorm(n = N, mean =  x[, t-1] + A * sin(2*pi*(t-1)/Time), sd = sd)
# }
#
# f1 <- function(x, t, Time, A, sd, N) {
#     rnorm(n = N, mean =  x[, t-1], sd = sd)
# }
#
# f2 <- function(x, t, Time, A, sd, N) {
#     # A * sin(2*pi*(t-1)/Time)
#     rnorm(n = N, mean =  A * sin(2*pi*(t-1)/Time), sd = sd)
# }



