##' Add a growth variable to X.
##'
##' @param X Covariate matrix. Must contain \code{par} variable.
##'
##' @return New covariate matrix with one more column.
##'
add_growth <- function(X){

  ## Sunlight variable
  par = X[,"par"]
  par = par - min(par)

  ## Starting to make growth variable.
  growth0 = par
  rng = range(par)
  cap = rng[1] + 0.05 * (rng[2] - rng[1])
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
        growth[ii] = growth[ii-1] + growth0[ii]
      }
    }
  }
  return(cbind(X, growth=growth))
}


##' Add a division varable to X.
##'
##' @param X Covariate matrix. Must contain \code{par} variable.
##'
##' @return New covariate matrix with one more column.
##'
add_div <- function(X){

  ## Sunlight variable
  par = X[,"par"]
  par = par - min(par)

  ## cap = 100
  rng = range(par)
  cap = rng[1] + 0.05*(rng[2]-rng[1])

  ## Also make a nightime variable
  div0 = rep(0, length(par))
  nighttime = which(par < cap)
  div0[nighttime] = -cap

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
  sine = sin(seq(from = 0, to = 2*pi, length = 25))[-25]
  cosine = cos(seq(from = 0, to = 2*pi, length = 25))[-25]

  ## Subset only the hours that exist.
  complete_sine = do.call(c, sapply(1:length(day_membership), function(ii){
    memb = day_membership[[ii]]
    hours = hh[memb]
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
    if(length(missing_hours_today) > 0){
      partial.cosine = cosine[-missing_hours_today]
      return(partial.cosine)
    } else {
      return(cosine)
    }
  }))
  return(cbind(X, sine = complete_sine, cosine = complete_cosine))
}
