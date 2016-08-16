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
plot_trajectories <- function(motiondata) {
    opar <- par()
    par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0),
        cex.lab = 1.5, font.lab = 2, cex.axis = 1.3, bty = "n", las = 1,
        lwd = 5)

    with(motiondata, {
        plot(observations ~ time, type = 'n', pch = 21, bg = "grey80", cex = 3,
             ylim = c(-20, 20), xlim = c(0, 2),
             col = rgb(red = 0, green = 0, blue = 0, alpha = 0.2),
             ylab = "Position [deg] ", xlab = "Time [s]")

        lines(acceleration ~ time, lty = "dotted")
        lines(velocity ~ time, lty = "solid")
        lines(position ~ time, lty = "dashed")
        points(time, observations,  col = "black", bg = 'white',
               cex = 1.5, pch = 21)
        # points(observations)
    })
    # text(0.5, 9, "Velocity", cex = 1.5, pos = 4)

    legend(0.2, -9, legend = c("Acceleration", "Velocity",
                               "Position", "Observations"),
           pch = c(rep(NA, 3), 21), col = rep("black", 4),
           bg = c(rep(NA, 3), "white"),
           lwd = rep(2.3, 4), lty = c("dotted", "solid", "dashed", NA),
           bty = "n", x.intersp = 0.5, cex = 1.5)
    # par(opar)
}

#' @export
plot_trajectories_2 <- function(motion_data, facet = FALSE) {
    color_palette <- c(
        "#000000",
        "#E69F00",
        "#56B4E9",
        "#009E73",
        "#F0E442",
        "#0072B2",
        "#D55E00",
        "#CC79A7"
    )

    ggplot2::theme_set(
        theme_classic() +
            # ggplot2::theme(
            #     axis.line.x = element_line(
            #         colour = 'black',
            #         size = 0.5,
            #         linetype = 'solid'
            #     ),
            #     axis.line.y = element_line(
            #         colour = 'black',
            #         size = 0.5,
            #         linetype = 'solid'
            #     )
            # ) +
            theme(legend.position = "right", text = element_text(size = 24))
    )

    set.seed(44234)

    data <- motion_data %>% tidyr::gather(

        key = "key",
        value = "value",
        -time,
        -observations,
        factor_key = TRUE
    )

    data <- data %>%
        dplyr::mutate(observations = ifelse(key %in% c("acceleration", "position"),
                                            NA, observations))

    data$key <- ordered(data$key,
                        levels = c("acceleration", "velocity", "position"))


    g <- ggplot(data = data, aes(
            x = time,
            y = value,
            linetype = key
        )) +

        geom_hline(yintercept = 0,
                   linetype = "dashed",
                   alpha = 0.4) +

        geom_line(size = 4, alpha = 0.5) +
        geom_line(
            data = dplyr::filter(data, key == "velocity"),
            alpha = 1.0,
            size = 4,
            linetype = "solid"
        ) +
        geom_point(
            aes(y = observations),
            alpha = 1.,
            fill = "white",
            colour = "white",
            shape = 21,
            size = 12
        ) +
        geom_point(
            aes(y = observations),
            alpha = 1.,
            fill = "white",
            colour = "black",
            shape = 21,
            size = 10
        ) +

        xlab("Time [s]") + ylab("Angular velocity [deg]") +

        scale_linetype_manual(
            values = c("dashed", "solid",
                       "dotted"),
            name = "",
            labels = c("acceleration", "velocity", "position"),
            guide = guide_legend(override.aes = list(
                alpha = c(0.5, 1, 1.0),
                title = NULL
            ))
        )

    if (facet) {
        g <- g + facet_grid(key ~ .) +
             theme(legend.position = "none")
    }
    g
}



#' @export
theme_publication <- function (base_size = 14, base_family = "") {
        theme_bw(base_size = base_size, base_family = base_family) %+replace%
            theme(
                panel.border = element_blank(),
                axis.line = element_line(colour = "black"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                strip.background = element_rect(colour = "black",
                                                size = 0.5),
                legend.key = element_blank(),
                # Tick labels
                axis.text.x = element_text(
                    size = rel(0.86),
                    colour = "black",
                    face = "bold"
                ),
                axis.text.y = element_text(
                    size = rel(0.86),
                    colour = "black",
                    face = "bold"
                ),

                # Axis
                axis.title = element_text(
                    size = rel(1),
                    colour = "black",
                    face = "bold"
                ),
                axis.line.x = element_line(colour = "black", size = 1),
                axis.line.y = element_line(colour = "black", size = 1),
                axis.ticks = element_line(colour = "black", size = 1),

                # Main title
                plot.title = element_text(
                    size = rel(1),
                    colour = "black" ,
                    lineheight = 1.0,
                    face = "bold"
                ),

                legend.position = "bottom",
                legend.title = element_text(
                    size = rel(0.7),
                    face = "bold",
                    colour = "black"
                ),
                legend.text = element_text(
                    size = rel(0.7),
                    face = "plain",
                    colour = "black"
                )
            )
    }
