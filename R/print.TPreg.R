print.TPreg <-
function(x,...){
  if(x$transition == "all"){
    cat("Call:\n")
    print(x$call)
    cat("\nTransitions:\n")
    print(c("11", "12", "13", "23"))
    cat("\n(s,t):\n")
    print(c(x$s, x$t))
    t11 <- length(x$co11$time)
    coeft <- x$co11$coefficients[t11, ]
    sdt <- x$co11$SD[t11, ]
    up.l <- as.numeric(coeft + 1.96*sdt)
    lw.l <- as.numeric(coeft - 1.96*sdt)
    p.v <- x$co11$p.value[t11, ]
    TAB11 <- cbind(Estimate = as.numeric(coeft), St.Err = as.numeric(sdt), LW.L = lw.l, UP.L = up.l, p.value = p.v)
    rownames(TAB11) <- names(coeft)
    colnames(TAB11) <- c("Estimate", "St.Err", "LW.L", "UP.L", "P.value")
    cat("\nCoefficients 11:\n")
    print(TAB11)
    t12 <- length(x$co12$time)
    coeft <- x$co12$coefficients[t12, ]
    sdt <- x$co12$SD[t12, ]
    up.l <- as.numeric(coeft + 1.96*sdt)
    lw.l <- as.numeric(coeft - 1.96*sdt)
    p.v <- x$co12$p.value[t12, ]
    TAB12 <- cbind(Estimate = as.numeric(coeft), St.Err = as.numeric(sdt), LW.L = lw.l, UP.L = up.l, p.value = p.v)
    rownames(TAB12) <- names(coeft)
    colnames(TAB12) <- c("Estimate", "St.Err", "LW.L", "UP.L", "P.value")
    cat("\nCoefficients 12:\n")
    print(TAB12)
    t13 <- length(x$co13$time)
    coeft <- x$co13$coefficients[t13, ]
    sdt <- x$co13$SD[t13, ]
    up.l <- as.numeric(coeft + 1.96*sdt)
    lw.l <- as.numeric(coeft - 1.96*sdt)
    p.v <- x$co13$p.value[t13, ]
    TAB13 <- cbind(Estimate = as.numeric(coeft), St.Err = as.numeric(sdt), LW.L = lw.l, UP.L = up.l, p.value = p.v)
    rownames(TAB13) <- names(coeft)
    colnames(TAB13) <- c("Estimate", "St.Err", "LW.L", "UP.L", "P.value")
    cat("\nCoefficients 13:\n")
    print(TAB13)
    t23 <- length(x$co23$time)
    coeft <- x$co13$coefficients[t23, ]
    sdt <- x$co13$SD[t23, ]
    zval <- coeft/sdt
    up.l <- as.numeric(coeft + 1.96*sdt)
    lw.l <- as.numeric(coeft - 1.96*sdt)
    p.v <- x$co23$p.value[t23, ]
    NMO <- x$n.misobs
    TAB23 <- cbind(Estimate = as.numeric(coeft), St.Err = as.numeric(sdt), LW.L = lw.l, UP.L = up.l, p.value = p.v)
    rownames(TAB23) <- names(coeft)
    colnames(TAB23) <- c("Estimate", "St.Err", "LW.L", "UP.L", "P.value")
    cat("\nCoefficients 23:\n")
    print(TAB23)
    cat("\n\n")
    print(paste(NMO, "observations deleted due to missingness from 'data'"))
  }
  else{
    t <- length(x$co$time)
    coeft <- x$co$coefficients[t, ]
    sdt <- x$co$SD[t, ]
    zval <- coeft/sdt
    up.l <- as.numeric(coeft + 1.96*sdt)
    lw.l <- as.numeric(coeft - 1.96*sdt)
    p.v <- x$co$p.value[t, ]
    NMO <- x$n.misobs
    TAB <- cbind(Estimate = as.numeric(coeft), St.Err = as.numeric(sdt), LW.L = lw.l, UP.L = up.l, p.value = p.v)
    rownames(TAB) <- names(coeft)
    colnames(TAB) <- c("Estimate", "St.Err", "LW.L", "UP.L", "P.value")
    cat("Call:\n")
    print(x$call)
    cat("\nTransition:\n")
    print(x$transition)
    cat("\n(s,t):\n")
    print(c(x$s, x$t))
    cat("\nCoefficients:\n")
    print(TAB)
    cat("\n\n")
    print(paste(NMO, "observations deleted due to missingness from 'data'"))
  }
}
