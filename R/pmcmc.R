#' @import ggplot2
#' @import dplyr
#' @importFrom viridis plasma

#' @export
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

#' @export
sd2precision <- function(sd) {
    prec <- 1/(sd^2)
    prec
}

#' @export
neff <- function(weights) {
    1/sum(weights^2)
}

#' @export
Mode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
}

#' @export
particle_filter <- function(data, N, Time, x_init, sdx_init,
                                params, resample = TRUE, rs_thresh = 0.5) {

    sdc <- params$sdc
    sdx <- params$sdx
    sdy <- params$sdy
    fun_x <- params$funx_x
    fun_c <- params$fun_c
    A <- params$A

    x <- matrix(rep(0, N*Time), nrow = N, ncol = Time)  #matrix(nrow =  N, ncol = Time)
    weights <- matrix(rep(0, N*Time), nrow = N, ncol = Time) #matrix(nrow =  N, ncol = Time)
    loglik <- rep(0, Time)

    x[, 1] <- rnorm(N, x_init, sdx_init)
    weights[, 1] <- dnorm(data$observations[1], x[, 1], sdy)
    weights[, 1] <- weights[, 1]/sum(weights[, 1])

    x[, 1] <- sample(x[, 1], replace = TRUE, size = N, prob = weights[, 1])


    for (t in seq(2, Time)) {

        x[, t] <- f(x, t, Time, A, sd = sdx, N)

        if (!is.na(data$observations[t])) {
            weights[, t] <- dnorm(data$observations[t], x[, t], sdy)
        } else {
            weights[, t] <- 1/N
        }


        if (resample) {
            if (neff(weights[, t]) < rs_thresh * N) {
                x[, t] <- sample(x[, t], replace = TRUE, size = N, prob = weights[, t])
            }
        }

        # TODO: do we need to reset weights after resampling?
        # weights[, t] <- 1/N

    }

    out <- list(x = x, weights = weights, loglik = loglik,
                x_means = apply(x, 2, mean),
                x_means = apply(x, 2, mean),
                x_means = apply(x, 2, mean),
                x_means = apply(x, 2, mean),
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


f <- function(x, t, Time, A, sd, N) {
    rnorm(n = N, mean =  x[, t-1] + A * sin(2*pi*(t-1)/Time), sd = sd)
}

f1 <- function(x, t, Time, A, sd, N) {
    rnorm(n = N, mean =  x[, t-1], sd = sd)
}

f2 <- function(x, t, Time, A, sd, N) {
    # A * sin(2*pi*(t-1)/Time)
    rnorm(n = N, mean =  A * sin(2*pi*(t-1)/Time), sd = sd)
}


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


pmcmc <- function(y, initPar, sigmav, sigmae, nPart, T, x0, nIter, stepSize) {

    #
    # Initialise variables
    #
    th     <- matrix(0, nrow=nIter, ncol=1)
    thp    <- matrix(0, nrow=nIter, ncol=1)
    ll     <- matrix(0, nrow=nIter, ncol=1)
    llp    <- matrix(0, nrow=nIter, ncol=1)
    accept <- matrix(0, nrow=nIter, ncol=1)

    # Set the initial parameter and estimate the initial log-likelihood
    th[1]  <- initPar
    ll[1]  <- sm(y, th[1], sigmav, sigmae, nPart, T, x0)$ll

    #
    # Run main loop
    #
    for (kk in 2:nIter) {

        # Propose a new parameter
        thp[kk] <- th[kk-1] + stepSize * rnorm(1)

        # Estimate the log-likelihood (don't run if unstable system)
        if (abs(thp[kk]) < 1.0) {
            llp[kk] <- smc(y, thp[kk], sigmav, sigmae, nPart, T, x0)$ll
        }

        # Compute the acceptance probability
        aprob <- exp(dnorm(thp[kk], log=TRUE) - dnorm(th[kk-1], log=TRUE) + llp[kk] - ll[kk-1])

        # Generate uniform random variable in U[0,1]
        u = runif(1)

        # Accept / reject step
        # Check if | phi | > 1.0, in that case always reject.
        if ((u < aprob) && ( abs( thp[kk] ) < 1.0 )) {
            # Accept the parameter
            th[kk]     <- thp[kk]
            ll[kk]     <- llp[kk]
            accept[kk] <- 1.0
        } else {
            # Reject the parameter
            th[kk]     <- th[kk-1]
            ll[kk]     <- ll[kk-1]
            accept[kk] <- 0.0
        }

        # Write out progress
        if (kk%%100 == 0) {
            cat(sprintf("#####################################################################\n"))
            cat(sprintf(" Iteration: %d of : %d completed.\n \n", kk, nIter))
            cat(sprintf(" Current state of the Markov chain:       %.4f \n", th[kk] ))
            cat(sprintf(" Proposed next state of the Markov chain: %.4f \n", thp[kk] ))
            cat(sprintf(" Current posterior mean:                  %.4f \n", mean(th[0:kk]) ))
            cat(sprintf(" Current acceptance rate:                 %.4f \n", mean(accept[0:kk]) ))
            cat(sprintf("#####################################################################\n"))
        }
    }

    #
    # Return traces of the parameters
    #
    th
}
