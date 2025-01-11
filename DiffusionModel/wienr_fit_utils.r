library(WienR)
library(ucminf)

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

fit_wienr <- function(rt, response, fit_sv = FALSE, fit_sw = FALSE, fit_st0 = FALSE, optim_control = list(), init_par = NULL, drift_index = NULL, bound_index = NULL, resid_index = NULL, ...) {
    if (is.character(response)) response <- as.numeric(as.factor(response))
    if (is.factor(response)) response <- as.numeric(response)
    
    par_names <- c()
    
    if (is.null(drift_index)) {
        drift_index <- rep(1, length(rt))
        n_drift <- 1
    } else {
        n_drift <- max(drift_index)
    }
    
    if (is.null(bound_index)) {
        bound_index <- rep(1, length(rt))
        n_bound <- 1
    } else {
        n_bound <- max(bound_index)
    }
    
    if (is.null(resid_index)) {
        resid_index <- rep(1, length(rt))
        n_resid <- 1
    } else {
        n_resid <- max(resid_index)
    }
    
    par_names <- c(paste0("a[", 1:n_bound, "]"), paste0("v[", 1:n_drift, "]"), paste0("w[", 1:n_bound, "]"), paste0("t0[", 1:n_resid, "]"))
    
    if (fit_sv) par_names <- c(par_names, paste0("sv[", 1:n_drift, "]"))
    if (fit_sw) par_names <- c(par_names, paste0("sw[", 1:n_bound, "]"))
    if (fit_st0) par_names <- c(par_names, paste0("st0[", 1:n_resid, "]"))
    
    init_to_use <- rep(NA, length(par_names))
    names(init_to_use) <- par_names
    
    ez_init <- ezddm(mean(response == 2), var(rt[response == 2]), mean(rt[response == 2]), length(response))
    
    init_to_use[startsWith(par_names, "a[")] <- unname(ez_init["a"])
    init_to_use[startsWith(par_names, "v[")] <- unname(ez_init["v"])
    init_to_use[startsWith(par_names, "w[")] <- 0.5
    init_to_use[startsWith(par_names, "t0[")] <- unname(min(0.99 * min(rt), ez_init["Ter"]))
    
    if (fit_sv) init_to_use[startsWith(par_names, "sv[")] <- 0
    if (fit_sw) init_to_use[startsWith(par_names, "sw[")] <- 0
    if (fit_st0) init_to_use[startsWith(par_names, "st0[")] <- 0
    
    if (!is.null(init_par)) {
        overlap <- intersect(names(init_to_use), names(init_par))
        init_to_use[overlap] <- init_par[overlap]
    }
    
    init_par <- init_to_use
    
    lower <- rep(-Inf, length(par_names))
    upper <- rep(Inf, length(par_names))
    
    lower[startsWith(par_names, "a[")] <- 0
    lower[startsWith(par_names, "w[")] <- 0
    lower[startsWith(par_names, "t0[")] <- 0
    lower[startsWith(par_names, "sv[")] <- 0
    lower[startsWith(par_names, "sw[")] <- 0
    lower[startsWith(par_names, "st0[")] <- 0
    
    upper[startsWith(par_names, "w[")] <- 1
    upper[startsWith(par_names, "t0[")] <- unname(tapply(rt, INDEX = resid_index, FUN = min))
    upper[startsWith(par_names, "sw[")] <- 1
    
    neg_log_likelihood <- function(par, rt, response, bound_index, drift_index, resid_index, ...) {
        a <- par[paste0("a[", bound_index, "]")]
        v <- par[paste0("v[", drift_index, "]")]
        w <- par[paste0("w[", bound_index, "]")]
        t0 <- par[paste0("t0[", resid_index, "]")]
        
        if (is.na(par["sv[1]"])) {
            sv <- 0
        } else {
            sv <- par[paste0("sv[", drift_index, "]")]
        }
        
        if (is.na(par["sw[1]"])) {
            sw <- 0
        } else {
            sw <- par[paste0("sw[", bound_index, "]")]
        }
        
        if (is.na(par["st0[1]"])) {
            st0 <- 0
        } else {
            st0 <- par[paste0("st0[", resid_index, "]")]
        }
        
        eval_pdf <- try(WienerPDF(t = rt, response = response, a = a, v = v, w = w, t0 = t0, sv = sv, sw = sw, st0 = st0, ...), silent = TRUE)
        
        if (any(class(eval_pdf) == "try-error")) return(Inf)
        
        return(-sum(eval_pdf$logvalue))
    }
    
    gradient <- function(par, rt, response, bound_index, drift_index, resid_index, ...) {
        a <- par[paste0("a[", bound_index, "]")]
        v <- par[paste0("v[", drift_index, "]")]
        w <- par[paste0("w[", bound_index, "]")]
        t0 <- par[paste0("t0[", resid_index, "]")]
        
        if (is.na(par["sv[1]"])) {
            sv <- 0
            use_sv <- FALSE
        } else {
            sv <- par[paste0("sv[", drift_index, "]")]
            use_sv <- TRUE
        }
        
        if (is.na(par["sw[1]"])) {
            sw <- 0
            use_sw <- FALSE
        } else {
            sw <- par[paste0("sw[", bound_index, "]")]
            use_sw <- TRUE
        }
        
        if (is.na(par["st0[1]"])) {
            st0 <- 0
            use_st0 <- FALSE
        } else {
            st0 <- par[paste0("st0[", resid_index, "]")]
            use_st0 <- TRUE
        }
        
        eval_grad <- try(gradWienerPDF(t = rt, response = response, a = a, v = v, w = w, t0 = t0, sv = sv, sw = sw, st0 = st0, ...), silent = TRUE)
        
        if (any(class(eval_grad) == "try-error")) return(rep(NaN, length(par)))
        
        eval_pdf <- try(WienerPDF(t = rt, response = response, a = a, v = v, w = w, t0 = t0, sv = sv, sw = sw, st0 = st0, ...), silent = TRUE)
        if (any(class(eval_pdf) == "try-error")) return(rep(NaN, length(par)))
        
        # Derivative of log(f(x)) is f'(x) / f(x)
        
        grad <- rep(NaN, length(par))
        names(grad) <- names(par)
        
        for (i in 1:max(bound_index)) {
            grad[paste0("a[", i, "]")] <- sum(eval_grad$deriv[bound_index == i, "da"] / eval_pdf$value[bound_index == i])
            grad[paste0("w[", i, "]")] <- sum(eval_grad$deriv[bound_index == i, "dw"] / eval_pdf$value[bound_index == i])
            if (use_sw) grad[paste0("sw[", i, "]")] <- sum(eval_grad$deriv[bound_index == i, "dsw"] / eval_pdf$value[bound_index == i])
        }
        
        for (i in 1:max(drift_index)) {
            grad[paste0("v[", i, "]")] <- sum(eval_grad$deriv[drift_index == i, "dv"] / eval_pdf$value[drift_index == i])
            if (use_sv) grad[paste0("sv[", i, "]")] <- sum(eval_grad$deriv[drift_index == i, "dsv"] / eval_pdf$value[drift_index == i])
        }
        
        for (i in 1:max(resid_index)) {
            grad[paste0("t0[", i, "]")] <- sum(eval_grad$deriv[resid_index == i, "dt0"] / eval_pdf$value[resid_index == i])
            if (use_st0) grad[paste0("st0[", i, "]")] <- sum(eval_grad$deriv[resid_index == i, "dst0"] / eval_pdf$value[resid_index == i])
        }
        
        return(-grad)
    }
    
    fit1 <- optim(
        par = init_par,
        fn = neg_log_likelihood,
        gr = gradient,
        method = "Nelder-Mead",
        control = c(optim_control, maxit = 10000),
        rt = rt,
        response = response,
        bound_index = bound_index,
        drift_index = drift_index,
        resid_index = resid_index
    )
    
    fit2 <- ucminf(
        par = fit1$par,
        fn = neg_log_likelihood,
        gr = gradient,
        control = optim_control,
        rt = rt,
        response = response,
        bound_index = bound_index,
        drift_index = drift_index,
        resid_index = resid_index
    )
    
    return(fit2)
}
