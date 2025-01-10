library(WienR)

ezddm <- function(propCorrect, rtCorrectVariance_seconds, rtCorrectMean_seconds, nTrials, s = 1) {
    # This function is adapted from the function of the same name in the `fddm` package
    s2 <- s^2
    v <- as.numeric(NA)
    a <- as.numeric(NA)
    Ter <- as.numeric(NA)
    propCorrect <- (2 * nTrials * propCorrect) / (2 * nTrials + 1) + 1 / (4 * nTrials * nTrials)
    L <- stats::qlogis(propCorrect)
    x <- L * (L * propCorrect^2 - L * propCorrect + propCorrect - 0.5) / rtCorrectVariance_seconds
    v <- sign(propCorrect - 0.5) * s * x^(1/4)
    a <- s2 * stats::qlogis(propCorrect)/v
    y <- -v * a/s2
    MDT <- (a/(2 * v)) * (1 - exp(y))/(1 + exp(y))
    Ter <- rtCorrectMean_seconds - MDT
    return(c("a" = a, "v" = v, "Ter" = Ter))
}

fit_wienr <- function(rt, response, fit_sv = FALSE, fit_sw = FALSE, fit_st0 = FALSE, optim_control = list(), optim_method = "L-BFGS-B", init_par = NULL, run_cleanup_fit = FALSE, ...) {
    if (is.character(response)) response <- as.numeric(as.factor(response))
    if (is.factor(response)) response <- as.numeric(response)
    
    ez_init <- ezddm(mean(response == 2), var(rt[response == 2]), mean(rt[response == 2]), length(response))
    
    if (is.null(init_par)) {
        init_par <- c("a" = unname(ez_init["a"]), "v" = unname(ez_init["v"]), "w" = 0.5, "t0" = unname(min(0.99 * min(rt), ez_init["Ter"])))
        
        if (fit_sv) {
            init_par <- c(init_par, "sv" = 0.5)
        }
        
        if (fit_sw) {
            init_par <- c(init_par, "sw" = 0.1)
        }
        
        if (fit_st0) {
            init_par <- c(init_par, "st0" = 0.1)
        }
    }
    
    if (optim_method %in% c("L-BFGS-B", "nlminb")) {
        lower <- c("a" = 0, "v" = -Inf, "w" = 0, "t0" = 0)
        upper <- c("a" = Inf, "v" = Inf, "w" = 1, "t0" = min(rt))
        
        if (fit_sv) {
            lower <- c(lower, "sv" = 0)
            upper <- c(upper, "sv" = Inf)
        }
        
        if (fit_sw) {
            lower <- c(lower, "sw" = 0)
            upper <- c(upper, "sw" = 0.5)
        }
        
        if (fit_st0) {
            lower <- c(lower, "st0" = 0)
            upper <- c(upper, "st0" = max(rt))
        }
    } else {
        lower <- rep(-Inf, length(init_par))
        upper <- rep(Inf, length(init_par))
    }
    
    neg_log_likelihood <- function(par, rt, response, ...) {
        a <- par["a"]
        v <- par["v"]
        w <- par["w"]
        t0 <- par["t0"]
        
        sv <- ifelse(is.na(par["sv"]), 0, par["sv"])
        sw <- ifelse(is.na(par["sw"]), 0, par["sw"])
        st0 <- ifelse(is.na(par["st0"]), 0, par["st0"])
        
        eval_pdf <- try(WienerPDF(t = rt, response = response, a = a, v = v, w = w, t0 = t0, sv = sv, sw = sw, st0 = st0, ...))
        
        if (any(class(eval_pdf) == "try-error")) return(Inf)
        
        return(-sum(eval_pdf$logvalue))
    }
    
    gradient <- function(par, rt, response, ...) {
        a <- par["a"]
        v <- par["v"]
        w <- par["w"]
        t0 <- par["t0"]
        
        sv <- ifelse(is.na(par["sv"]), 0, par["sv"])
        sw <- ifelse(is.na(par["sw"]), 0, par["sw"])
        st0 <- ifelse(is.na(par["st0"]), 0, par["st0"])
        
        eval_grad <- try(gradWienerPDF(t = rt, response = response, a = a, v = v, w = w, t0 = t0, sv = sv, sw = sw, st0 = st0, ...))
        
        if (any(class(eval_grad) == "try-error")) return(rep(NaN, length(par)))
        
        eval_pdf <- try(WienerPDF(t = rt, response = response, a = a, v = v, w = w, t0 = t0, sv = sv, sw = sw, st0 = st0, ...))
        if (any(class(eval_pdf) == "try-error")) return(rep(NaN, length(par)))
        
        # Derivative of log(f(x)) is f'(x) / f(x)
        
        grad <- c()
        
        for (par_name in names(par)) {
            grad <- c(grad, sum(eval_grad$deriv[,paste0("d", par_name)] / eval_pdf$value))
        }
        
        names(grad) <- names(par)
        
        return(-grad)
    }
    
    if (optim_method == "ucminf") {
        require(ucminf)
        
        fit <- ucminf(
            par = init_par,
            fn = neg_log_likelihood,
            gr = gradient,
            control = optim_control,
            rt = rt,
            response = response
        )
    } else if (optim_method == "nlminb") {
        fit <- nlminb(
            start = init_par,
            objective = neg_log_likelihood,
            gradient = gradient,
            control = optim_control,
            lower = lower,
            upper = upper,
            rt = rt,
            response = response
        )
    } else {
        fit <- optim(
            par = init_par,
            fn = neg_log_likelihood,
            gr = gradient,
            method = optim_method,
            lower = lower,
            upper = upper,
            control = optim_control,
            rt = rt,
            response = response
        )
    }
    
    if (run_cleanup_fit) {
        fit <- nlminb(
            start = fit$par,
            objective = neg_log_likelihood,
            gradient = gradient,
            control = list(eval.max = 10000, iter.max = 10000),
            lower = lower,
            upper = upper,
            rt = rt,
            response = response
        )
    }
    
    return(fit)
}