TPreg <-
function(formula, data, link, s = 0, t = NULL, R = 199, by = NULL, trans)
{
  if (missing(data)) 
    stop("Argument 'data' is missing with no default")
  if (!is.data.frame(data)) 
    stop("Argument 'data' must be a data.frame")
  if (sum(! c("id","Zt","Tt","delta1","delta" ) %in% (colnames(data)))>0)
    stop("data should  contain  id, Zt, Tt, delta1, delta variables")
  if (ncol(data)<=4)
    stop("'data' must have covariables")
  if(sum(is.na(data))!=0){
    miscolumn <- sapply(data, function(x) sum(is.na(x)))
    miscolnam <- names(miscolumn[miscolumn!=0])
    if(sum(miscolumn)!=0) 
      warning(sapply(miscolnam, function(x)paste(x," variable in 'data' has missing value(s)", ", ",sep="")))
  }
  
  formula <- formula
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m  <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())   
  mt <- attr(mf, "terms")
  if (is.empty.model(mt)) {
    stop("'formula' must match with 'data'")
  }
  else{
    X<- model.matrix(mt, mf, contrasts)
    ind = match(colnames(X) , colnames(data))
  ind = ind[!is.na(ind)]
  covname= colnames(data)[ind]
}

  ordata = data[, c("id", "Zt", "Tt", "delta1", "delta", covname )]
  L.or <- nrow(ordata)
  comdata <- ordata[complete.cases(ordata),]
  L.com <- nrow(comdata)
  n.misobs <- L.or - L.com
  if(is.null(by)){
    by <- floor((max(comdata$Zt) - min(comdata$Zt))/quantile(comdata$Zt,0.01))
  }
  
  if(is.null(t)){
    t = max(comdata$Zt[comdata$delta1 == 1], na.rm = T)
  }
  if(t <= s || s<0){
    stop("argument 's' must be smaller than 't' and larger than 0")
  }
  else{
    if (!(link %in% c("logit", "probit", "cauchit"))) 
      stop( paste("binomial family do not have", "'", link, "'",  "link"))
    
    if(sum(comdata$delta1 < comdata$delta) != 0){ 
      stop("'delta' must be 0 when 'delta1' is 0")
    }
    if (!(trans %in% c("11", "12", "13","23","all")))
      stop(paste(trans, "is not a valid transition for a progressive illness-death model")) 
    
    if(s == 0 & (trans == "23" || trans == "all" ))
      stop("for the transition '23' argument 's' must be larger than 0")
    if (trans == "23" || trans == "all" ){
      X2 <- X[comdata$Zt<=s & comdata$Tt>s ,]
      data2 <- comdata[comdata$Zt <= s & comdata$Tt > s,] 
      Sfit23 <- summary( survfit(Surv(data2$Tt, data2$delta == 0)~ +1))
      Shat23 <- rbind(c(0,1),data.frame(time = Sfit23$time, surv = Sfit23$surv))
      Shat.function23 <- function(x){
        Shatx23 <- if(length(Shat23$time) > 1)  tail(subset(Shat23, time <= x)$surv,1) 
        else 1
        return(Shatx23)
      }
    }
    if (trans=="11" ||trans=="12" ||trans=="13" || trans=="all"){
      X <- X[comdata$Zt > s, ]  
      data1 <- comdata[comdata$Zt > s,] 
      Sfit1 <- summary( survfit(Surv(data1$Zt, data1$delta1 == 0)~ +1))
      Shat1 <- rbind(c(0,1), data.frame(time = Sfit1$time, surv = Sfit1$surv))
      Shat.function1 <- function(x){
        Shatx1 <- if(length(Shat1$time) > 1)  tail(subset(Shat1,time <= x)$surv, 1) 
        else 1
        return(Shatx1)
      }
      Sfit <- summary( survfit(Surv(data1$Tt, data1$delta == 0)~ +1))
      Shat <- rbind(c(0,1), data.frame(time = Sfit$time, surv = Sfit$surv))
      Shat.function <- function(x){
        Shatx <- if(length(Shat$time) > 1)  tail(subset(Shat, time <=x )$surv, 1) 
        else 1
        return(Shatx)
      }
    }
    co <- vector("list", 4)
    names(co) <- c("co11", "co12", "co13", "co23")
    if(trans == "11" || trans == "all" ){
      vec.t11 <- data1$Zt
      M11 <- max(vec.t11)
      vec.t11 <- vec.t11[order(vec.t11[vec.t11 > s])]
      vec.t11<- vec.t11[vec.t11 <= t]
      vec.t11<- vec.t11[seq(1, length(vec.t11), by)]
      vec.t11 <- c(vec.t11, t)
      vec.t11 <- unique(vec.t11)
      L.t11 <- length(vec.t11)
      if(vec.t11[L.t11] >= M11) 
        stop("for the tansition '11' the effects can not be estimated for the given 't'(large 't' returns all responses equal to 0) ")
      iii <-NULL
      registerDoParallel(cores = 4)
      eta.list <- lapply( vec.t11, function(x){ 
        jumptime <- x
        res <- (data1$Zt > jumptime)       
        delta1_t <- ifelse(data1$Zt <= jumptime , data1$delta1, 1)
        hatG1 <- sapply(pmin(data1$Zt, jumptime), Shat.function1) 
        wei <- delta1_t/hatG1
        X <- X
        family <- binomial(link)   
        fit <- mod.glm.fit( X, res, family = family, weights = wei)
        eta <- coef(fit)
        data1 <- data1
        r <- foreach(1:R, .combine=rbind,.export=c("iii","mod.glm.fit")) %dopar% {
          X <- X
          iboot <- sample(1:nrow(data1), replace=TRUE)
          iii <- rbind(iii, iboot)
          boot.data <- data1[iboot, ]
          result <- mod.glm.fit(X[iboot, ], res[iboot], family = family, weights = wei[iboot])
          coefficients(result)
        }
        boot.eta <- r
        boot.sd <- apply(boot.eta, 2, sd, na.rm = TRUE)
        return(c(eta, boot.sd))
      })
      
      eta <- data.frame(do.call("rbind", eta.list))
      coef <- eta[ , 1:(length(eta[1,])/2)]
      sd <- eta[ , (length(eta[1, ])/2 + 1):(length(eta[1, ]))]
      colnames(sd) <- colnames(coef)
      CO <- list(transition = "11", time = vec.t11, coefficients = coef, SD = sd, LWL = coef - 1.96*sd, UPL = coef + 1.96*sd, p.value = 2*pnorm(-abs(as.matrix(coef/sd))))
      if(trans == "all"){
        co$co11 = CO
      }
      else {
        co <- list("co" = CO, call = match.call(),transition = trans, s = s, t = t, n.misobs = n.misobs)
        class(co) = "TPreg" 
        return(co)
      }
    }
    
    if(trans == "12" || trans == "all"){
      index <- data1$Zt < data1$Tt 
      vec.t12 <- c(data1$Zt[index], data1$Tt[index])
      M12 <- max(vec.t12)
      vec.t12 <- vec.t12[order(vec.t12[vec.t12 > s])]
      vec.t12 <- vec.t12[vec.t12 <= t]
      vec.t12 <- vec.t12[seq(1, length(vec.t12), by)]
      vec.t12 <- c(vec.t12, t)
      vec.t12 <- unique(vec.t12)
      L.t12 <- length(vec.t12)
      if(vec.t12[L.t12] >= M12) stop(" for the tansition '12' the effects can not be estimated for the given 't', (large 't' returns all responses equal to 0)")
      iii <- NULL
      eta.list <- lapply( vec.t12, function(x){ 
        jumptime <- x
        res <- (data1$Zt <= jumptime & jumptime < data1$Tt)
        delta_t <- ifelse(data1$Tt <= jumptime , data1$delta, 1)
        hatG <- sapply(pmin(data1$Tt, jumptime), Shat.function) 
        wei <- delta_t/hatG
        family <- binomial(link)   
        fit <- mod.glm.fit( X,res,family=family,weights=wei)
        eta <- coef(fit)
        data1 <- data1
        formula1 <- formula
        X <- X
        r <- foreach(1:R, .combine=rbind,.export=c("iii","mod.glm.fit")) %dopar% {
          iboot <- sample(1:nrow(data1), replace=TRUE)
          iii <- rbind(iii, iboot)
          boot.data <- data1[iboot, ]
          result <- mod.glm.fit(X[iboot,],res[iboot],family=family,weights=wei[iboot])
          coefficients(result)
        }
        boot.eta <- r
        boot.sd <- apply(boot.eta, 2, sd, na.rm = TRUE)
        return(c(eta, boot.sd))
      })
      eta <- data.frame(do.call("rbind", eta.list))
      coef <- eta[ , 1:(length(eta[1, ])/2)]
      sd <- eta[ ,(length(eta[1, ])/2+1):(length(eta[1, ]))]
      colnames(sd) <- colnames(coef)
      CO = list( transition = "12", time = vec.t12, coefficients = coef, SD = sd, LWL = coef - 1.96*sd,UPL=coef+1.96*sd, p.value = 2*pnorm(-abs(as.matrix(coef/sd))))
      if(trans == "all"){
        co$co12 = CO
      }
      else {
        co <- list("co" = CO, call = match.call(), transition = trans, s = s, t = t, n.misobs = n.misobs)
        class(co) = "TPreg" 
        return(co)
      }
    }
    
    if(trans == "13" || trans == "all"){
      vec.t13 <- data1$Tt
      M13 <- max(vec.t13)
      vec.t13 <- vec.t13[order(vec.t13[vec.t13 > s])]
      vec.t13 <- vec.t13[vec.t13 <= t]
      vec.t13 <- vec.t13[seq(1, length(vec.t13), by)]
      vec.t13 <- c(vec.t13, t)
      vec.t13 <- unique(vec.t13)
      L.t13 <- length(vec.t13)
      if(vec.t13[L.t13] >= M13) stop(" for the tansition '13' the effects can not be estimated for the given 't', (large 't' returns all responses equal to 1) ")
      iii <- NULL
      eta.list <- lapply( vec.t13, function(x){ 
        jumptime <- x
        res <- (data1$Tt <= jumptime)
        delta_t <- ifelse(data1$Tt <= jumptime , data1$delta, 1)
        hatG <- sapply(pmin(data1$Tt, jumptime), Shat.function) 
        wei <- delta_t/hatG
        family <- binomial(link)   
        fit <- mod.glm.fit( X, res, family = family, weights=wei)
        eta <- coef(fit)
        data1 <- data1
        formula1 <- formula
        X <- X
        r <- foreach(1:R, .combine=rbind, .export = c("iii","mod.glm.fit")) %dopar% {
          iboot <- sample(1:nrow(data1), replace=TRUE)
          
          iii <- rbind(iii, iboot)
          boot.data <- data1[iboot, ]
          result <- mod.glm.fit(X[iboot, ], res[iboot], family = family, weights = wei[iboot])
          coefficients(result)
        }
        boot.eta <- r
        boot.sd <- apply(boot.eta, 2, sd, na.rm = TRUE)
        return(c(eta, boot.sd))
      })
      eta <- data.frame(do.call("rbind", eta.list))
      coef <- eta[ ,1:(length(eta[1, ])/2)]
      sd <- eta[ ,(length(eta[1, ])/2 + 1):(length(eta[1, ]))]
      colnames(sd) <- colnames(coef)
      CO <- list(transition = "13", time = vec.t13, coefficients = coef, SD = sd, LWL = coef - 1.96*sd, UPL = coef + 1.96*sd, p.value = 2*pnorm(-abs(as.matrix(coef/sd))))
      if(trans == "all"){
        co$co13 = CO
      }
      else {
        co <- list("co" = CO, call = match.call(), transition = trans, s = s, t = t, n.misobs = n.misobs)
        class(co) = "TPreg" 
        return(co)
      }
    }
    if(trans == "23" || trans == "all"){
      vec.t23 <- data2$Tt
      M23 <- max(vec.t23)
      vec.t23 <- vec.t23[order(vec.t23[vec.t23 > s])]
      vec.t23 <- vec.t23[vec.t23 <= t]
      vec.t23 <- vec.t23[seq(1, length(vec.t23), by)]
      vec.t23 <- c(vec.t23,t)
      vec.t23 <- unique(vec.t23)
      L.t23 <- length(vec.t23)
      if(vec.t23[L.t23] >= M23) stop(" for the tansition '23' the effects can not be estimated for the given 't'(large 't' returns all responses equal to 1)")
      iii <- NULL
      eta.list <- lapply( vec.t23, function(x){ 
        jumptime <- x
        res <-(data2$Tt <= jumptime)
        delta_t <- ifelse(data2$Tt <= jumptime, data2$delta, 1)
        hatG <- sapply(pmin(data2$Tt, jumptime), Shat.function23) 
        wei <- delta_t/hatG
        family <- binomial(link)   
        fit <- mod.glm.fit( X2, res, family = family, weights=wei)
        eta <- coef(fit)
        data2 <- data2
        formula1 <- formula
        X2 <- X2
        r <- foreach(1:R, .combine = rbind, .export=c("iii","mod.glm.fit")) %dopar% {
          iboot <- sample(1:nrow(data2), replace=TRUE)
          iii <- rbind(iii, iboot)
          boot.data <- data2[iboot, ]
          result <- mod.glm.fit(X2[iboot, ], res[iboot], family = family, weights = wei[iboot])
          coefficients(result)
        }
        boot.eta <- r
        boot.sd <- apply(boot.eta,2, sd, na.rm = TRUE)
        return(c(eta, boot.sd))
      })
      eta <- data.frame(do.call("rbind", eta.list))
      coef <- eta[ , 1:(length(eta[1, ])/2)]
      sd <- eta[ , (length(eta[1, ])/2 + 1):(length(eta[1, ]))]
      colnames(sd) <- colnames(coef)
      CO <- list(transition = "23", time = vec.t23, coefficients = coef, SD = sd, LWL = coef - 1.96*sd, UPL = coef + 1.96*sd, p.value = 2*pnorm(-abs(as.matrix(coef/sd))))
      if(trans == "all"){
        co$co23=CO
      } 
      else {
        co <- list("co" = CO, call = match.call(), transition = trans, s = s, t = t, n.misobs=n.misobs)
        class(co)="TPreg" 
        return(co)
      }
    }
    if(trans == "all"){
      co$call = match.call()
      co$transition = "all"
      co$s = s
      co$t = t
      co$n.misobs = n.misobs
      class(co) = "TPreg"
      return(co)
    }
  }
}
