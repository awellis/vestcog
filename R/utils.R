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
plot_trajectories <- function(motion_data, facet = FALSE) {
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


    g <-
        ggplot(data = data, aes(
            x = time,
            y = value,
            linetype = key
        )) +

        geom_line(size = 2, alpha = 0.5) +
        geom_line(
            data = dplyr::filter(data, key == "velocity"),
            alpha = 1.0,
            size = 2,
            linetype = "solid"
        ) +
        geom_point(
            aes(y = observations),
            alpha = 1.,
            fill = "white",
            colour = "white",
            shape = 21,
            size = 8
        ) +
        geom_point(
            aes(y = observations),
            alpha = 1.,
            fill = "white",
            colour = "black",
            shape = 21,
            size = 6
        ) +

        geom_hline(yintercept = 0,
                   linetype = "dashed",
                   alpha = 0.4) +
        xlab("Time [s]") + ylab("Angular velocity [deg]") +
        # facet_grid(key ~ .) +
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
    print(g)
}

plot_trajectories_old <- function(data, velocity_only = FALSE) {

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
            theme(legend.position = "right", text = element_text(size = 24))
    )

    set.seed(44234)
    data <- data %>% tidyr::gather(key = "key", value = "value", -time)
    data$key <- factor(data$key)

    if (velocity_only) {
        keylist <- "velocity"
    } else {
        keylist <- c("acceleration", "velocity", "position")
    }

   g <- ggplot(data = data, aes(x = time, y = value, linetype = key)) +

        geom_line(size = 2, alpha = 0.5) +
        geom_line(data = dplyr::filter(data, key == "velocity"),
                  alpha = 1.0, size = 2, linetype = "solid") +

        geom_point(data = dplyr::filter(data, key == "observations"), alpha = 1.,
                   fill = "white", colour = "white", shape = 21, size = 8) +
        geom_point(data = dplyr::filter(data, key == "observations"), alpha = 1.,
                   fill = "white", colour = "black", shape = 21, size = 6) +

        geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.4) +

       facet_grid(.~ key) +

        xlab("Time [s]") + ylab("Angular velocity [deg]") +

       facet_grid(key ~ .) +

        # scale_shape_manual(name = "", guide = guide_legend(override.aes = list(
        #     values = c(NULL, 21, NULL, NULL)))) +

        scale_linetype_manual(values = c("twodash", "dashed",
                                         "dotted", "solid"),
                              name = "", labels = c("acceleration", "observations", "position", "velocity"),
                              guide = guide_legend(override.aes = list(
                        alpha = c(0.5, 0.5, 0.5, 1.0)),
                        title = NULL))
   print(g)
}

