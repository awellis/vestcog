

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
