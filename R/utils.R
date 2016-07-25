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
generate_data <- function(T = 2, dt = 0.1, amplitude = 20, sensor_sd = 1.7,
                          as_df = TRUE, seed = TRUE) {
    nsteps <- T/dt
    t <- seq(from = 0, to = T, length.out = nsteps)

    # the following generates a motion profile with single-cycle sinusoidal
    # acceleration
    time <- t
    position <- amplitude*T/(2*pi) * (t-(T/(2*pi)) * sin(2*pi*t/T))
    velocity <- amplitude * T/(2 * pi) * (1-cos(2 * pi * t/T))
    acceleration <- amplitude * sin(2 * pi * t/T)
    trajectory <- rbind(position, velocity, acceleration)

    if(seed) {set.seed(4879863)}
    observations <- rnorm(ncol(trajectory), trajectory[2,], sensor_sd)
    out <- rbind(time, trajectory, observations)

    if(as_df) {
        out <- out %>% t() %>% data.frame()
    }

    return(out)
}

#' @export
plot_trajectories <- function(data, velocity_only = FALSE) {

    color_palette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

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
            theme(legend.position = "none", text = element_text(size = 28))
    )

    set.seed(44234)
    data <- data %>% tidyr::gather(key = "key", value = "value", -time)

    if (velocity_only) {
        keylist <- "velocity"
    } else {
        keylist <- c("acceleration", "velocity", "position")
    }

    g <- ggplot(data = data, aes(x = time, y = value, color = key)) +
                                 # linetype = key)) +
        # geom_line(data = filter(data, key %in% keylist),
        #           size = 2) +
        geom_line(data = filter(data, key %in% "position"),
                  size = 2, linetype = "longdash", color = "black", alpha = 0.5) +
        geom_line(data = filter(data, key %in% "velocity"),
                  size = 2, linetype = "solid", color = "black",  alpha = 1.0) +
        geom_line(data = filter(data, key %in% "acceleration"),
                  size = 2, linetype = "dotdash", color = "black",  alpha = 0.5) +

        # geom_point(data = filter(data, key == "observations"),
        #            size = 3, shape = 15, color = "black", alpha = 0.5) +
        # geom_line(data = filter(data, key == "observations"),
        #           colour = "black", size = 1, alpha = 0.2, linetype = "dotted") +

        geom_line(data = filter(data, key == "observations"), colour = "darkgrey", linetype = "solid", size = 1.5, alpha = 1) +
        geom_point(data = filter(data, key == "observations"), alpha = 1., fill = "white", colour = "white", shape = 21, size = 6) +
        geom_point(data = filter(data, key == "observations"), alpha = 1., fill = "white", colour = "grey40", shape = 21, size = 4) +


        # ggplot2::geom_point(aes(data = filter(data, key == "observations"),
        #                         y = observations), alpha = 0.8, fill = "white", colour = "white", shape = 21, size = 6) +
        # ggplot2::geom_point(aes(data = filter(data, key == "observations"),
        #                         y = observations), alpha = 0.8, fill = "white", colour = "grey40", shape = 21, size = 4) +
        # ggplot2::geom_line(aes(data = filter(data, key == "observations"),
        #                        y = median), colour = "steelblue",
        #                    linetype = "solid", size = 2, alpha = 1) +


    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.4) +

        xlab("Time [s]") + ylab("Angular velocity [deg]") +

        # labs(title = "Natural head motion") +
        # scale_colour_brewer(palette = "Set1")
        scale_colour_manual(values =
                                c(color_palette[7], color_palette[2], color_palette[1])) + #sample(color_palette)) +\
        # facet_grid(~ key) +
        theme(legend.title = element_blank())
    print(g)
}