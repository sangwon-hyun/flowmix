##' Add a growth variable to X.
add_growth <- function(X){
  ## Sunlight variable
  par = X[,"par"]
  ## plot(par, type='l', ylim = c(-2000,5000))

  ## Starting to make growth variable.
  growth0 = par
  cap = 100
  growth0[growth0 >= cap] = cap


  ## Use a "cumsum" loop to make growth variable.
  nighttime = which(growth0!=cap)
  growth = rep(NA, length(par))
  growth[1] = growth0[1]
  for(ii in 1:length(growth0)){
    if(growth0[ii] < 1E-3){
      growth[ii] = par[ii]
    } else {
      if(ii > 1){
        ## print(ii)
        growth[ii] = growth[ii-1] + growth0[ii]
      }
    }
  }
  ## points(growth, type='l', lwd=3, col='green')
  return(cbind(X, growth=growth))
}


##' Add a division varable to X.
add_div <- function(X){
  ## Sunlight variable
  par = X[,"par"]
  ## plot(par, type='l', ylim = c(-2000,5000))

  ## Starting to make growth variable.
  growth0 = par
  cap = 100
  growth0[growth0 >= cap] = cap


  ## Also make a nightime variable
  div0 = rep(0, length(par))
  nighttime = which(par < 100)
  div0[nighttime] = -cap
  ## points(cumsum(div0), type = 'l', lwd=3, col='yellow')

  div = rep(NA, length(par))
  div[1] = div0[1]
  for(ii in 1:length(div0)){
    if(div0[ii] > -1E-3){
      div[ii] = 0
    } else {
      if(ii > 1){
        div[ii] = div[ii-1] + div0[ii]
      }
    }
  }
  ## points(div, type='l', lwd=3, col='yellow')
  return(cbind(X, div=div))
}


##' Adding *hourly* sine and cosine curve to the data matrix X. Assumes that the
##' X matrix has "mm", "dd", and "hh" columns.
##'
##' @param X Data matrix.
##' @param time A vector of time, the same length as X.
##'
##' @return X with augmented columns.
add_sine_to_X <- function(X, time){

  ## Get year, month, day and hour.
  strings = time
  yyyy = as.numeric(sapply(strings, substr, 1, 4))
  mm = as.numeric(sapply(strings, substr, 6, 7))
  dd = as.numeric(sapply(strings, substr, 9, 10))
  hh = as.numeric(sapply(strings, substr, 12, 13))

  ## Order the time points
  days = unique(mm * 33 + dd)
  day_membership = lapply(days, function(day){
    which(mm * 33 + dd == day)
  })
  length(unlist(day_membership))

  ## Create full sine and cosine
  sine = sin(seq(from=0, to=2*pi, length=25))[-25]
  cosine = cos(seq(from=0, to=2*pi, length=25))[-25]

  ## Subset only the hours that exist.
  complete_sine = do.call(c, sapply(1:length(day_membership), function(ii){
    memb = day_membership[[ii]]
    hours = hh[memb]
    ## print(length(memb))
    ## print(length(hours))
    missing_hours_today = which(0:23 %ni% hours)
    if(length(hours) < 24){
      partial.sine = sine[-missing_hours_today]
      return(partial.sine)
    } else {
      return(sine)
    }
  }))

  ## Do the same for cosine.
  complete_cosine = do.call(c, sapply(1:length(day_membership), function(ii){
    memb = day_membership[[ii]]
    hours = hh[memb]
    missing_hours_today = which(0:23 %ni% hours)
    if(length(missing_hours_today)>0){
      partial.cosine = cosine[-missing_hours_today]
      return(partial.cosine)
    } else {
      return(cosine)
    }
  }))
  return(cbind(X, sine=complete_sine, cosine=complete_cosine))
}
