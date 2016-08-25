

particle_filter_active <- function(data, Time, params, resample_particles = TRUE, rs_thresh = 0.5) {

    # unpack parameters
    sd_c <- params$sd_c
    sd_x <- params$sd_x
    sd_y <- params$sd_y
    fun <- params$fun
    sd_active <- params$sd_active
    sd_passive <- params$sd_passive
    A <- params$A
    N <- params$N
    x_init <- params$x_init
    sd_x_init <- params$sd_x_init

    # initialize variables
    x <- matrix(rep(0, N*Time), nrow = N, ncol = Time)  #matrix(nrow =  N, ncol = Time)
    weights <- matrix(rep(0, N*Time), nrow = N, ncol = Time) #matrix(nrow =  N, ncol = Time)
    loglik <- rep(0, Time)

    # sample from prior distribution
    x[, 1] <- rnorm(N, x_init, sd_x_init)
    weights[, 1] <- update(data$observations[1], x[, 1], sd_y)
    loglik[1] <- log(sum(weights[, 1]))

    weights[, 1] <- normalize(weights[, 1])
    x[, 1] <- resample(x[, 1], weights[, 1])

    for (t in seq(2, Time)) {
        # predict 1 step ahead using process model as proposal distribution
        x[, t] <- fun(x, t, Time, A, sd = sd_x, N)

        if (!is.na(data$observations[t])) {
            weights[, t] <- update(data$observations[t], x[, t], sd_y)
        } else {
            weights[, t] <- 1/N
        }
        loglik[t] <- log(sum(weights[, t]))

        if (resample_particles) {
            if (neff(weights[, t]) < rs_thresh * N) {
                x[, t] <- resample(x[, t], weights[, t])
            }
        }
    }
    logliksum <- sum(loglik) - log(N)
    out <- list(x = x, weights = weights, loglik = loglik, logliksum = logliksum,
                x_means = apply(x, 2, mean),
                x_medians = apply(x, 2, median),
                x_quantiles = apply(x, 2,
                                    function(x) quantile(x,
                                                         probs = c(0.025,
                                                                   0.975))),
                Time = Time,
                N = N)
}