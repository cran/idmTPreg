iddata <-
function(Stime, Sind , Iltime, Ilind, covar,...){
  
  if (missing(Stime)) 
    stop("Must have a Stime argument")
  if (!is.numeric(Stime)) 
    stop("Survival time variable is not numeric")
  if (missing(Sind))
    stop("Must have an Sind argument")
  
  
  if(sum(!Sind %in% c(0,1))==0){
    delta=Sind
  }
  else if  (is.logical(Sind)) {
    delta <- as.numeric(Sind)
  }
  else  stop("Invalid Sind value")
  
  if(sum(Iltime > Stime)!=0) 
    stop("Illness time can not be larger than Survival time")
  
  if (missing(Iltime)) 
    stop("Must have a Iltime argument")
  if (!is.numeric(Iltime)) 
    stop("Illness  time variable is not numeric")
  
  if (missing(Ilind))
    stop("Must have an Ilind argument")
  
  Zt <- pmin(Iltime, Stime)
  Tt <- Stime
  id <- 1:length(Zt)
  delta1 <- ifelse(delta==0 & Ilind==0,0,1)
  
  iddata <- data.frame(id, Zt, delta1, Tt, delta, covar, ...)
  class(iddata) <- c("iddata","data.frame")
  return(iddata)
}
