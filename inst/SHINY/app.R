# app.R
# Shiny MLE Lab: Normal, Exponential, Poisson, Gamma, Weibull
# ------------------------------------------------------------
# Features:
# - Simulate or paste your own data
# - Compute MLEs (optim L-BFGS-B), SEs (observed Fisher info via numDeriv::hessian)
# - Show fitted vs true curves/bars + log-likelihood visualization
# - 95% Wald confidence intervals

library(shiny)
library(ggplot2)
library(numDeriv)

# ---------- Helpers: log-likelihoods, simulation, density/pmf, and defaults -----

# Each dist gets a small spec with:
#  name, param_names, simulate(), density()/pmf(), loglik(), bounds (lower), initials()

dists <- list(

  normal = list(
    label = "Normal(μ, σ)",
    param_names = c("mu","sigma"),
    lower = c(-Inf, 1e-8),
    simulate = function(n, pars) rnorm(n, mean = pars["mu"], sd = pars["sigma"]),
    d_fun = function(x, pars) dnorm(x, mean = pars["mu"], sd = pars["sigma"]),
    loglik = function(par, x) {
      mu <- par[1]; sigma <- par[2]
      if (sigma <= 0) return(-Inf)
      sum(dnorm(x, mu, sigma, log = TRUE))
    },
    initials = function(x) c(mu = mean(x), sigma = sd(x))
  ),

  exponential = list(
    label = "Exponential(λ)  [rate]",
    param_names = c("lambda"),
    lower = c(1e-10),
    simulate = function(n, pars) rexp(n, rate = pars["lambda"]),
    d_fun = function(x, pars) dexp(x, rate = pars["lambda"]),
    loglik = function(par, x) {
      lambda <- par[1]
      if (lambda <= 0 || any(x < 0)) return(-Inf)
      sum(dexp(x, rate = lambda, log = TRUE))
    },
    initials = function(x) c(lambda = 1 / mean(x))
  ),

  poisson = list(
    label = "Poisson(λ)",
    param_names = c("lambda"),
    lower = c(1e-10),
    simulate = function(n, pars) rpois(n, lambda = pars["lambda"]),
    d_fun = function(x, pars) dpois(x, lambda = pars["lambda"]),
    loglik = function(par, x) {
      lambda <- par[1]
      if (lambda <= 0 || any(x < 0) || any(x != floor(x))) return(-Inf)
      sum(dpois(x, lambda = lambda, log = TRUE))
    },
    initials = function(x) c(lambda = mean(x))
  ),

  gamma = list(
    label = "Gamma(α shape, β rate)",
    param_names = c("alpha","beta"),
    lower = c(1e-10, 1e-10),
    simulate = function(n, pars) rgamma(n, shape = pars["alpha"], rate = pars["beta"]),
    d_fun = function(x, pars) dgamma(x, shape = pars["alpha"], rate = pars["beta"]),
    loglik = function(par, x) {
      alpha <- par[1]; beta <- par[2]
      if (alpha <= 0 || beta <= 0 || any(x < 0)) return(-Inf)
      sum(dgamma(x, shape = alpha, rate = beta, log = TRUE))
    },
    initials = function(x) {
      m <- mean(x); v <- var(x)
      # method of moments
      alpha <- m^2 / v
      beta  <- m / v
      c(alpha = max(alpha, 1e-3), beta = max(beta, 1e-3))
    }
  ),

  weibull = list(
    label = "Weibull(k shape, λ scale)",
    param_names = c("shape","scale"),
    lower = c(1e-10, 1e-10),
    simulate = function(n, pars) rweibull(n, shape = pars["shape"], scale = pars["scale"]),
    d_fun = function(x, pars) dweibull(x, shape = pars["shape"], scale = pars["scale"]),
    loglik = function(par, x) {
      shape <- par[1]; scale <- par[2]
      if (shape <= 0 || scale <= 0 || any(x < 0)) return(-Inf)
      sum(dweibull(x, shape = shape, scale = scale, log = TRUE))
    },
    initials = function(x) {
      # crude initials: log-moment heuristic
      lx <- log(pmax(x, 1e-8))
      s  <- sd(lx); m <- mean(x)
      shape <- max(1.2 / (s + 1e-8), 0.5)
      scale <- max(m, 1e-3)
      c(shape = shape, scale = scale)
    }
  )
)

# Utility: optimize MLE and compute SEs via observed Fisher information
mle_fit <- function(dist_key, x) {
  spec <- dists[[dist_key]]
  k <- length(spec$param_names)
  init <- spec$initials(x)
  stopifnot(length(init) == k)
  # bounds
  lower <- spec$lower
  upper <- rep(Inf, k)
  # NEGATIVE log-lik for minimization
  nll <- function(par) -spec$loglik(par, x)
  opt <- tryCatch(
    optim(par = init, fn = nll, method = "L-BFGS-B", lower = lower, upper = upper, control = list(fnscale = 1)),
    error = function(e) NULL
  )
  if (is.null(opt) || opt$convergence != 0) {
    return(list(ok = FALSE, message = "Optimization failed to converge.", par = NA, vcov = NA))
  }
  # Hessian of +loglik (so we negate nll)
  ll <- function(par) spec$loglik(par, x)
  H <- tryCatch(hessian(func = ll, x = opt$par), error = function(e) NULL)
  vc <- NA
  ses <- rep(NA, k)
  if (!is.null(H)) {
    # observed information is -H (since H is second derivative of loglik),
    # Var(theta_hat) approx inv(-H)
    info <- -H
    # add small ridge if needed
    eig <- tryCatch(eigen(info, symmetric = TRUE)$values, error = function(e) NULL)
    if (!is.null(eig) && all(eig > 1e-8)) {
      vc <- solve(info)
      ses <- sqrt(diag(vc))
    }
  }
  list(ok = TRUE, par = setNames(opt$par, spec$param_names), vcov = vc, se = setNames(ses, spec$param_names))
}

# Make a tidy table from MLE output + 95% Wald CI
mle_table <- function(fit) {
  if (!isTRUE(fit$ok)) {
    return(data.frame(Parameter = character(), Estimate = numeric(), SE = numeric(), `95% CI Lower` = numeric(), `95% CI Upper` = numeric()))
  }
  est <- fit$par
  se  <- fit$se
  z   <- 1.96
  lo  <- est - z * se
  hi  <- est + z * se
  data.frame(
    Parameter = names(est),
    Estimate = as.numeric(est),
    SE = as.numeric(se),
    `95% CI Lower` = as.numeric(lo),
    `95% CI Upper` = as.numeric(hi),
    row.names = NULL,
    check.names = FALSE
  )
}

# ---------- Shiny UI ------------------------------------------------------------

ui <- fluidPage(
  titlePanel("Maximum Likelihood Explorer: 5 Univariate Distributions"),
  sidebarLayout(
    sidebarPanel(
      selectInput("dist", "Distribution", choices = setNames(names(dists), sapply(dists, `[[`, "label"))),
      numericInput("n", "Sample size (n)", value = 200, min = 5, step = 1),
      checkboxInput("use_custom", "Use my own data (paste numeric values below)", FALSE),
      conditionalPanel(
        condition = "input.use_custom",
        helpText("Paste a comma/space/newline-separated vector of numeric values (integers for Poisson)."),
        textAreaInput("custom_data", NULL, rows = 4, placeholder = "e.g. 0.2, 0.5, 1.1, 0.7, 1.9")
      ),
      hr(),
      strong("True parameters (used for simulation & 'true' overlay)"),
      uiOutput("param_inputs"),
      numericInput("seed", "Random seed", value = 123, step = 1),
      actionButton("simulate", "Generate / Analyze", class = "btn-primary"),
      width = 4
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Data & Fit",
                 br(),
                 plotOutput("dataPlot", height = 380),
                 br(),
                 tableOutput("mleTable")
        ),
        tabPanel("Log-Likelihood",
                 br(),
                 plotOutput("llPlot", height = 420),
                 helpText("For 1-parameter models: a 1D log-likelihood curve with the MLE marked.",
                          "For 2-parameter models: log-likelihood contours with the MLE (★) and true params (×).")
        ),
        tabPanel("Details",
                 br(),
                 h4("Notes"),
                 tags$ul(
                   tags$li("MLEs are computed with bounded optimization (L-BFGS-B) in natural parameter space."),
                   tags$li("Standard errors and 95% Wald intervals use the observed Fisher information (the negative Hessian of the log-likelihood at the MLE), computed numerically via ", code("numDeriv::hessian"), "."),
                   tags$li("Discrete distributions (Poisson) show bars; continuous distributions show density curves.")
                 ),
                 h4("Tips"),
                 tags$ul(
                   tags$li("Toggle between simulating data and pasting your own to see how the MLE responds."),
                   tags$li("Increase n to watch the likelihood concentrate and the SEs shrink (consistency)."),
                   tags$li("Compare true parameters to MLEs and CIs to build intuition about estimator behavior.")
                 )
        )
      )
    )
  )
)

# ---------- Shiny Server --------------------------------------------------------

server <- function(input, output, session) {

  # Dynamic parameter inputs per distribution
  output$param_inputs <- renderUI({
    key <- req(input$dist)
    spec <- dists[[key]]
    labs <- spec$param_names
    # Provide sensible defaults for the "true" parameters
    defaults <- switch(key,
                       normal     = c(mu = 0, sigma = 1),
                       exponential= c(lambda = 1),
                       poisson    = c(lambda = 3),
                       gamma      = c(alpha = 2, beta = 1),
                       weibull    = c(shape = 1.5, scale = 1),
                       c(1,1)
    )
    tagList(lapply(seq_along(labs), function(i) {
      nm <- labs[i]
      val <- defaults[i]
      minv <- if (grepl("sigma|lambda|alpha|beta|shape|scale", nm)) 1e-6 else -Inf
      numericInput(paste0("par_", nm), nm, value = val, min = minv)
    }))
  })

  # Reactive: data vector x
  data_vec <- eventReactive(input$simulate, {
    key <- req(input$dist)
    spec <- dists[[key]]
    set.seed(input$seed %||% 123)
    # Collect true parameters from UI
    par_vals <- setNames(
      sapply(spec$param_names, function(nm) input[[paste0("par_", nm)]] %||% NA_real_),
      spec$param_names
    )
    validate(need(all(is.finite(par_vals)), "Please supply finite parameter values."))

    if (isTRUE(input$use_custom)) {
      raw <- gsub("[,\\n\\t]", " ", input$custom_data %||% "")
      xs  <- suppressWarnings(as.numeric(strsplit(raw, "\\s+")[[1]]))
      xs  <- xs[is.finite(xs)]
      validate(need(length(xs) >= 5, "Please paste at least 5 numeric values."))
      # For Poisson, enforce integers
      if (key == "poisson") {
        validate(need(all(abs(xs - round(xs)) < 1e-8 & xs >= 0), "Poisson data must be nonnegative integers."))
      }
      list(x = xs, true_par = par_vals)
    } else {
      n <- input$n %||% 200
      validate(need(n >= 5, "Sample size must be at least 5."))
      xs <- spec$simulate(n, par_vals)
      list(x = xs, true_par = par_vals)
    }
  }, ignoreInit = TRUE)

  # Reactive: MLE fit
  fit_obj <- reactive({
    dv <- req(data_vec())
    key <- req(input$dist)
    mle_fit(key, dv$x)
  })

  # Table of estimates, SEs, CIs
  output$mleTable <- renderTable({
    fit <- req(fit_obj())
    tab <- mle_table(fit)
    if (!isTRUE(fit$ok)) {
      return(data.frame(Message = fit$message))
    }
    tab
  }, digits = 5)

  # Data + fitted/true graph
  output$dataPlot <- renderPlot({
    dv <- req(data_vec())
    x  <- dv$x
    key <- req(input$dist)
    spec <- dists[[key]]
    fit <- req(fit_obj())

    # Continuous vs discrete plot
    is_discrete <- key %in% c("poisson")

    if (!is_discrete) {
      df <- data.frame(x = x)
      p <- ggplot(df, aes(x)) +
        geom_histogram(aes(y = ..density..), bins = 30, fill = "grey85", color = "grey40") +
        labs(title = paste0("Data with Fitted & True Density — ", dists[[key]]$label),
             x = "x", y = "Density")
      # Overlay true density
      xs <- seq(min(x) - 0.2 * sd(x), max(x) + 0.2 * sd(x), length.out = 400)
      true_y <- spec$d_fun(xs, dv$true_par)
      p <- p + geom_line(aes(x = xs, y = true_y), linewidth = 0.9, linetype = "dashed")
      # Overlay fitted density if available
      if (isTRUE(fit$ok)) {
        fit_y <- spec$d_fun(xs, fit$par)
        p <- p + geom_line(aes(x = xs, y = fit_y), linewidth = 1.1)
      }
      p
    } else {
      # Poisson bars
      kmin <- max(0, floor(min(x)))
      kmax <- max(ceiling(max(x)), kmin + 10)
      kseq <- kmin:kmax
      df_obs <- as.data.frame(table(factor(x, levels = kseq)))
      names(df_obs) <- c("k", "count")
      df_obs$k <- as.integer(as.character(df_obs$k))
      df_obs$rel <- df_obs$count / sum(df_obs$count)

      p <- ggplot(df_obs, aes(k, rel)) +
        geom_col(width = 0.8, fill = "grey80", color = "grey40") +
        labs(title = paste0("Observed Relative Frequencies vs PMF — ", dists[[key]]$label),
             x = "k", y = "Relative Frequency / PMF")

      # Overlay true PMF
      true_p <- dpois(kseq, lambda = dv$true_par["lambda"])
      p <- p + geom_point(aes(kseq, true_p), shape = 1, size = 3) +
        geom_line(aes(kseq, true_p), linetype = "dashed")

      # Overlay fitted PMF if available
      if (isTRUE(fit$ok)) {
        fit_p <- dpois(kseq, lambda = fit$par["lambda"])
        p <- p + geom_point(aes(kseq, fit_p))
      }
      p
    }
  })

  # Log-likelihood plot
  output$llPlot <- renderPlot({
    dv <- req(data_vec())
    x  <- dv$x
    key <- req(input$dist)
    spec <- dists[[key]]
    fit <- req(fit_obj())
    true_par <- dv$true_par

    # 1 parameter -> 1D curve
    if (length(spec$param_names) == 1) {
      nm <- spec$param_names[1]
      hat <- if (isTRUE(fit$ok)) fit$par[[nm]] else NA_real_
      # Build a grid around MLE / true
      center <- if (is.finite(hat)) hat else true_par[[nm]]
      span <- if (key %in% c("poisson","exponential")) 4 else 2
      lo <- max(spec$lower[1],
                if (center > 0) center / (1 + span) else 1e-6)
      hi <- center * (1 + span)
      if (!is.finite(lo) || lo <= 0) lo <- 1e-6
      if (!is.finite(hi) || hi <= lo) hi <- lo * 20

      grid <- seq(lo, hi, length.out = 300)
      llv  <- sapply(grid, function(g) spec$loglik(g, x))

      df <- data.frame(theta = grid, ll = llv)
      p <- ggplot(df, aes(theta, ll)) +
        geom_line(linewidth = 1) +
        labs(x = nm, y = "Log-likelihood", title = paste0("Log-likelihood — ", spec$label))
      # MLE and true markers
      if (isTRUE(fit$ok)) {
        p <- p + geom_vline(xintercept = fit$par[[nm]], linetype = "solid") +
          annotate("text", x = fit$par[[nm]], y = max(llv), vjust = -0.5, label = "MLE", fontface = 2)
      }
      p <- p + geom_vline(xintercept = true_par[[nm]], linetype = "dashed") +
        annotate("text", x = true_par[[nm]], y = min(llv), vjust = 1.5, label = "True", fontface = 1)
      p
    } else {
      # 2-parameter -> contour
      nms <- spec$param_names
      # Grid around MLE or true
      center <- if (isTRUE(fit$ok)) fit$par else true_par
      span   <- c(0.6, 0.6)  # multiplicative span
      g1 <- seq(max(spec$lower[1], center[1] * (1 - span[1])), center[1] * (1 + span[1]), length.out = 80)
      g2 <- seq(max(spec$lower[2], center[2] * (1 - span[2])), center[2] * (1 + span[2]), length.out = 80)
      g1 <- g1[is.finite(g1) & g1 > spec$lower[1]]
      g2 <- g2[is.finite(g2) & g2 > spec$lower[2]]
      if (length(g1) < 20 || length(g2) < 20) {
        # Fallback: absolute ranges
        g1 <- seq(max(spec$lower[1], 1e-3), max(spec$lower[1], 1e-3) * 50, length.out = 80)
        g2 <- seq(max(spec$lower[2], 1e-3), max(spec$lower[2], 1e-3) * 50, length.out = 80)
      }
      M <- outer(g1, g2, Vectorize(function(a, b) spec$loglik(c(a, b), x)))
      df <- expand.grid(p1 = g1, p2 = g2)
      df$ll <- as.vector(M)
      p <- ggplot(df, aes(p1, p2, z = ll)) +
        geom_contour(aes(z = ll), bins = 15) +
        labs(x = nms[1], y = nms[2], title = paste0("Log-likelihood Contours — ", spec$label))
      # MLE and true markers
      if (isTRUE(fit$ok)) {
        p <- p + geom_point(aes(x = fit$par[1], y = fit$par[2]), shape = 8, size = 3) +
          annotate("text", x = fit$par[1], y = fit$par[2], label = "★ MLE", vjust = -1)
      }
      p <- p + geom_point(aes(x = true_par[1], y = true_par[2]), shape = 4, size = 3) +
        annotate("text", x = true_par[1], y = true_par[2], label = "× True", vjust = 1.5)
      p
    }
  })
}

shinyApp(ui, server)

