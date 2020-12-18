Bayesian longitudinal functional data analysis example
================
John Shamshoian
December 17, 2020

This notebook runs a longitudinal-functional analysis on simulated data.

    library(tidyverse)
    library(MASS)
    library(Rcpp)
    library(splines)
    library(LFBayes)
    library(plotly)

    ### Control parameters

    set.seed(999)
    source("/Users/johnshamshoian/Documents/R_projects/LFBayes/Example/Simulation_funcs.R")
    errorvar <- .025
    SS <- 20
    TT <- 20
    t <- seq(from = 0, to = 1, length.out = TT)
    s <- seq(from = 0, to = 1, length.out = SS)
    n <- 60
    tt <- list()
    tt[[1]] <- 1:(TT*SS)
    tt <- rep(tt, n)
    p1 <- 12
    p2 <- 12
    q1 <- 3
    q2 <- 3
    Bt <- bs(t, df = p1, intercept = TRUE)
    Bs <- bs(s, df = p2, intercept = TRUE)
    Bt1 <- bs(t, df = p1, intercept = TRUE)
    Bs1 <- bs(s, df = p2, intercept = TRUE)

    ### Generate key model quantities

    H <- GenerateH(q1, q2)
    mu1 <- GenerateMu1(s,t)
    Lambda <- Loading.Matern(t, p1, q1, Bt)
    Gamma <- Loading.Brown.Bridge(s, p2, q2)
    Cov <- kronecker(Bs%*%Gamma, Bt%*%Lambda)%*%H%*%
      t(kronecker(Bs%*%Gamma, Bt%*%Lambda)) + errorvar * diag(SS * TT)

    ### Generate data from model

    x <- mvrnorm(n, mu  = as.vector(mu1), Sigma = Cov)
    sx <- sd(x)
    mx <- mean(x)
    x <- (x-mx)/sx
    Smooth_scaled_cov <- (Cov - errorvar * diag(SS * TT)) / sx^2
    mu <- (mu1 - mx)/sx
    y <- lapply(1:n, function(i) x[i,])
    missing <- list()
    for(ii in 1:n){
      missing[[ii]] <- numeric(0)
    }
    X <- cbind(rep(1, n))

    ### Visualize a few trajectories

    nsub <- 6
    small_data <- tibble(id = numeric(),
                         func_time = numeric(),
                         long_time = numeric(),
                         value = numeric())
    small_data <- small_data %>% 
      add_row(id = rep(1:nsub, each = SS * TT),
              func_time = rep(rep(t, SS), nsub),
              long_time = rep(rep(s, each = TT), nsub),
              value = c(t(x[1:nsub,])))

    id.labs <- paste("Subject", 1:nsub)
    names(id.labs) <- "1":nsub
    small_data %>%
      filter(long_time %in% c(s[5], s[10], s[15])) %>%
      ggplot(aes(func_time, value)) +
      geom_line(aes(linetype = factor(long_time))) +
      facet_wrap(. ~ id, labeller = labeller(id = id.labs)) + 
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
      labs(x = "Functional time", y = "Response", linetype = "Longitudinal time") +
      scale_linetype_discrete(labels = round(c(s[5], s[10], s[15]), 2))

![](example_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

    ### MCMC control parameters

    iter <- 500 # Number of iterations
    burnin <- 100 # Burnin iterations
    thin <- 1 # Thinning for each chain
    nchain <- 1 # Number of chains
    q1s <- 3 # Number of latent factors for functional dimension
    q2s <- 3 # Number of latent factors for longitudinal dimension
    alpha <- .05 # Type 1 error
    neig <- 3 # Number of eigenfunctions for inference

    ### Processing

    MCMC <- mcmcWeakChains(y, missing, X, Bs1, Bt1,
                           q1s, q2s, iter, thin, burnin, nchain)
    MCMC_eigen <- eigenLFChains(Bs1, Bt1, MCMC, neig,
                                iter, burnin, nchain, s, t, alpha)

    ### Prepare data for plotting
    postmean <- matrix(MCMC_eigen$postmean, nrow = TT, ncol = SS)
    fig <- plot_ly(z = ~postmean, x = s, y = t) %>% add_surface()
    eigenfunction_tibble <- tibble(time = numeric(), type = character(),
                                   value = numeric(), bound = character(),
                                   number = numeric())
    eigenfunction_tibble <- eigenfunction_tibble %>%
      add_row(value = c(MCMC_eigen$eigvecFuncmean,
                        MCMC_eigen$eigvecFunclower,
                        MCMC_eigen$eigvecFuncupper),
              number = rep(rep(1:neig, each = TT), 3),
              bound = rep(c("mean", "lower", "upper"), each = neig * TT),
              type = "Functional",
              time = rep(t, 3 * neig)) %>%
      add_row(value = c(MCMC_eigen$eigvecLongmean,
                        MCMC_eigen$eigvecLonglower,
                        MCMC_eigen$eigvecLongupper),
              number = rep(rep(1:neig, each = SS), 3),
              bound = rep(c("mean", "lower", "upper"), each = neig * SS),
              type = "Longitudinal",
              time = rep(s, 3 * neig))

    ### Plot mean surface

    axx <- list(
      title = "Longitudinal time"
    )

    axy <- list(
      title = "Functional time"
    )

    axz <- list(
      title = "Posterior mean"
    )

    fig <- plot_ly() %>%
      add_surface(z = ~postmean, x = s, y = t) %>%
      layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz,
                          aspectratio = list(x = .9, y = .8, z = 1))) %>%
      hide_colorbar()
    htmlwidgets::saveWidget(fig, "index.html")
    htmltools::tags$iframe(getwd(), "index.html")

<!--html_preserve-->
<iframe>
/Users/johnshamshoian/Documents/R\_projects/LFBayes/Example index.html
</iframe>
<!--/html_preserve-->

    ### Plot eigenfunctions

    options(repr.plot.width=6, repr.plot.height=3)
    number.labs <- paste("Eigenfunction", 1:neig)
    names(number.labs) <- c("1":neig)
    eigenfunction_tibble %>%
      ggplot(aes(time, value)) +
      geom_line(aes(linetype = bound)) +
      facet_grid(type ~ number, labeller = labeller(number = number.labs)) +
      scale_linetype_manual(values=c("dashed", "solid", "dashed"))+
      theme_bw() +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

![](example_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

    ### Percent variability explained by first few eigenfunctions in functional direction
    100 * rowMeans(MCMC_eigen$eigvalFunc)[p1:(p1-neig+1)]

    ## [1] 63.379161 24.612260  7.079978

    ### Percent variability explained by first few eigenfunctions in longitudinal direction
    100 * rowMeans(MCMC_eigen$eigvalLong)[p2:(p2-neig+1)]

    ## [1] 78.402375 14.563600  2.790165