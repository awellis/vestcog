#' @import ggplot2
#' @import dplyr
#' @importFrom viridis plasma

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
particle_filter <- function(data, Time, params, resample_particles = TRUE, rs_thresh = 0.5) {

    # unpack parameters
    sd_c <- params$sd_c
    # sd_x <- params$sd_x
    sd_x <- params$sd_x
    sd_y <- params$sd_y
    fun <- params$fun

    # fun_c <- params$fun_c
    A <- params$A
    N <- params$N
    x_init <- params$x_init
    sd_x_init <- params$sd_x_init

    # initialize variables
    x <- matrix(rep(0, N*Time), nrow = N, ncol = Time)  #matrix(nrow =  N, ncol = Time)
    weights <- matrix(rep(0, N*Time), nrow = N, ncol = Time) #matrix(nrow =  N, ncol = Time)
    loglik <- rep(0, Time)
    ll <- 0

    # sample from prior distribution
    x[, 1] <- rnorm(N, x_init, sd_x_init)
    weights[, 1] <- update(data$observations[1], x[, 1], sd_y)
    loglik[1] <- log(sum(weights[, 1])) - log(N)
    ll <- ll + log(sum(weights[, 1])) - log(N)

    weights[, 1] <- normalize(weights[, 1])
    # x[, 1] <- resample(x[, 1], weights[, 1])

    for (t in seq(2, Time)) {

        if (resample_particles) {
            if (neff(weights[, t]) < rs_thresh * N) {
                x[, t] <- resample(x[, t], weights[, t])
            }
        }

        # predict 1 step ahead using process model as proposal distribution
        x[, t] <- fun(x, t, Time, A, sd = sd_x, N)

        if (!is.na(data$observations[t])) {
            weights[, t] <- update(data$observations[t], x[, t], sd_y)
        } else {
            weights[, t] <- 1/N
        }

        loglik[t] <- log(sum(weights[, t])) - log(N)
        ll <- ll + log(sum(weights[, t])) - log(N)
    }

    logliksum <- sum(loglik)
    # logliksum <- sum(loglik) - log(N)


    out <- list(x = x, weights = weights, loglik = loglik,
                logliksum = logliksum, ll = ll,
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
plot_filtering_estimates <- function(object, data, predict = FALSE) {

    color_palette <- viridis::plasma(n = 9)
    color_palette <- viridis::plasma(n = 9)

    df <- with(object, {
        dplyr::data_frame(t = seq(1, Time),
                          mean = x_means,
                          median = x_medians,
                          lower = x_quantiles[1, ],
                          upper = x_quantiles[2, ],
                          x_true = data$velocity,
                          observations = data$observations,
                          loglik = loglik)
    })

    p <- ggplot2::ggplot(data = df, aes(x = t)) +
        ggplot2::geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1,
                             fill = color_palette[2]) +
        ggplot2::geom_line(aes(y = x_true), colour = color_palette[7], alpha = 0.9,
                           linetype = "dashed", size = 1.2) +
        # geom_line(aes(y = mean), colour = color_palette[6], size = 1.4) +
        ggplot2::geom_line(aes(y = median), colour = color_palette[1], size = 1.4, alpha = 0.6) +

        ggplot2::geom_point(aes(y = observations), colour = "black",
                            size = 3, shape = 15, alpha = 0.5) +
        ggplot2::geom_line(aes(y = observations), colour = "black", size = 1,
                           alpha = 0.2, linetype = "dotted") +
        ggplot2::scale_x_continuous(limits = c(1, 20), breaks = c(0, 5, 10, 15, 20),
                                    labels = c("0", "0.5", "1", "1.5", "2")) +
        ggplot2::ylab(expression(paste("Latent state: ", omega))) + xlab("Time [sec]")

    if (predict) {
        p <- p +
            ggplot2::geom_vline(xintercept = 10, linetype = "dashed",
                                color = color_palette[4], size = 1.5, alpha = 0.6) +
            ggplot2::geom_ribbon(data = filter(df, t > 9),
                                 aes(ymin = lower, ymax = upper), alpha = 0.4,
                                 fill = color_palette[4])
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


smc <- function(y, phi, sigmav, sigmae, nPart, T, x0) {
    #
    # Fully-adapted particle filter for the linear Gaussian SSM
    #
    # Inputs:
    # y:                   observations from the system for t=1,...,T.
    #
    # phi, sigmav, sigmae: the persistence of the state and the
    #                      standard deviations of the state innovations and
    #                      observation noise.
    #
    # nPart:               number of particles (N)
    #
    # T and xo:            the no. observations and initial state.
    #
    # Outputs:
    # xh:                  vector with T elements
    #                      estimates of the filtered state
    #                      for each t=0,1,...,T-1.
    #
    # ll:                  estimate of the log-likelihood at T-1
    #
    #

    #
    # Initialise variables

    xhatf <- matrix(x0, nrow=T, ncol=1)
    p     <- matrix(x0, nrow=nPart, ncol=T+1)
    w     <- matrix(1/nPart, nrow=nPart, ncol=T+1)
    v     <- matrix(1, nrow=nPart, ncol=T+1)
    ll    <- 0

    #
    # Run main loop
    #
    for (tt in 2:T) {


        # Resample ( multinomial )

        nIdx   <- sample(1:nPart, nPart, replace=TRUE, prob = w[,tt-1])


        # Propagate

        Delta  <- ( sigmav^(-2) + sigmae^(-2) )^(-1)
        mup    <- sigmae^(-2) * y[tt] + sigmav^(-2) * phi * p[nIdx,tt-1]
        p[,tt] <- Delta * mup + rnorm(nPart, 0, sqrt(Delta))


        # Compute weights

        v[,tt] <- dnorm(y[tt+1], phi * p[,tt], sqrt( sigmae^2 + sigmav^2), log=TRUE)

        # Rescale log-weights and recover weight
        vmax   <- max(v[,tt])
        v[,tt] <- exp(v[,tt] - vmax)

        # Normalize the weights
        w[,tt] <- v[,tt] / sum(v[,tt])

        # Estimate the state
        xhatf[tt] <- mean(p[,tt])

        # Estimate the log-likelihood
        ll        <- ll + vmax + log(sum(v[,tt])) - log(nPart)

    }

    #
    # Return state estimate and log-likelihood estimate
    #
    list(xh = xhatf, ll=ll, p=p, w=w)
}

