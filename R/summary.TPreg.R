summary.TPreg <-
function(object,...){
  cat("Call:\n")
  print(object$call)
  if (object$transition == "all") {
    nlist <- 4
    nam <- c("11", "12", "13", "23")
    cat("Transitions:\n")
    print(nam)
  }
  else {
    nam <- object$transition
    nlist <- 1
  }
  cat("(s,t):\n")
  print(c(object$s, object$t))
  for (i in 1:nlist){
    ilist <- object[[i]]
    cat("\n",paste("Transition ", nam[i], sep=""), " :\n")
    cat("\n     Coefficients:\n")
    print(cbind(time = ilist$time, ilist$coefficients))
    cat("\n     Standard Errors:\n")
    print(cbind(time = ilist$time, ilist$SD))
    cat("\n     Lower limit:\n")
    print(cbind(time = ilist$time, ilist$LWL))
    cat("\n     Upper limit:\n")
    print(cbind(time = ilist$time, ilist$UPL))
    cat("\n     p.value:\n")
    print(data.frame(cbind(time = ilist$time, ilist$p.value)))
    cat("\n\n")
  }
  W.n.misobs <- paste( object$n.misobs, "observation(s) deleted due to missingness from 'data'")
  print(W.n.misobs)  
}
