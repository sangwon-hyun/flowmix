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



########################################################################################
## All just convenience functions, temporarily. Nothing in this file will be retained ##
########################################################################################

## Let me reproduce the data here:
make_data75 <- function(){
  ## Load the hourly level X
  datadir = "~/repos/cruisedat/export"
  load(file.path(datadir, "MGL1704-hourly-only-binned.Rdata"))

  ## For back-compatibility; X now doesn't have any NA par values.
  if(any(is.na(X[,"par"]))){
    exclude_ind = which(is.na(X[,"par"]))
    X = X[-exclude_ind,]
  }

  ## Create the TF bases.
  TT = nrow(X)
  ind = rep(0, TT)
  north = which(X[,"lat"] > 37)
  regions = sapply(1:3, function(ii){
    endpt = c(0,range(north),TT)[ii:(ii+1)]
    (endpt[1]+1):endpt[2]
  })

  load(file.path(datadir, "MGL1704-hourly-only-binned-1d-diam.Rdata"))
  print("Hourly time resolution with bases and lagged PAR")
  print("Only sss, sst, par are being used.")
  X = X[,which(colnames(X) %in% c("sss", "sst", "par", "sss_cruise", "sst_cruise"))]
  print("Artificial covariates are being added.")

   ## Create the lagged sunlight variable
  lags = c(0,3,6,9,12)
  par = scale(X[,"par"])
  par = par - min(par)
  par = par / max(par)
  parlist = lapply(lags, function(lag)lagpad(par, lag))

  ## Make the additional columns
  dat = do.call(cbind, parlist)
  colnames(dat) = paste0("p", 0:4)

  ## Combine them with X
  X = X[,-which(colnames(X) == "par")]
  X = cbind(X, dat)
  X = scale(X)

  ## Add bases
  bases0 = lapply(regions, function(reg){
    vec = rep(0, TT)
    vec[reg] = 1
    vec
  })
  bases = do.call(cbind, bases0)
  ## bases = do.call(cbind, bases0[2:3]) ## Omitting one basis.
  ## colnames(bases) = c("b2", "b3")
  colnames(bases) = c("b1", "b2", "b3")
  X = cbind(X, bases)

  ## Rid of NA rows.
  na.rows = which(apply(dat, 1, function(myrow)any(is.na(myrow))))
  if(length(na.rows) > 0){
    X = X[-na.rows,]
    ylist = ylist[-na.rows]
    countslist = countslist[-na.rows]
    ylist = lapply(ylist, cbind)
  }

  maxdev = 0.1155245
  return(list(ylist=ylist,
              X=X,
              countslist=countslist))
}





##' Adding noise to covariates.
##'
##' @param Xorig Original covariate matrix
##' @param sigmalist Numeric vector of Gaussian noise standard deviations.
##' @param noise_ii Index of \code{sigmalist} to use.
##'
##' @return A noisy version of \code{Xorig}.
##'
##' @export
add_noise <- function(Xorig, sigmalist, noise_ii){

  X = Xorig

  if(noise_ii > 0){
    ## Add noise to the sunlight covariate
    X[,"par"] = Xorig[,"par"] + rnorm(TT, 0, sigmalist[noise_ii])

    ## Add noise to the spurious covariates by the same amount as well.
    for(icol in 3:ncol(X)){
      X[,icol] = Xorig[,icol] + rnorm(TT, 0, sigmalist[noise_ii])
    }
  }
  return(X)
}


##' Creates a directory \code{destin}, if it doesn't already exist.
##' @export
create_destin <- function(destin){
  if(!dir.exists(destin)){
    dir.create(destin, recursive = TRUE)
    cat("Creating destin: ", destin, fill=TRUE)
  } else {
    cat("All output goes out to destin: ", destin, fill = TRUE)
  }
}


##' Helper to divide up the jobs in \code{iimat} into a total of
##' \code{arraynum_max} jobs. The purpose is to divide the jobs, in order to run
##' this on a server.
##'
##' @param arraynum_max Maximum SLURM array number.
##' @param iimat matrix whose rows contain CV job indices.
##'
##' @export
make_iilist <- function(arraynum_max, iimat){
  iimax = nrow(iimat)
  if(arraynum_max > iimax){
    iilist = lapply(1:iimax, function(a)a)
  } else {
    ends = round(seq(from=0, to=iimax, length=arraynum_max+1))
    iilist = Map(function(a,b){ (a+1):b}, ends[-length(ends)], ends[-1])
    stopifnot(length(unlist(iilist)) == nrow(iimat))
  }
  stopifnot(length(unlist(iilist)) == nrow(iimat))
  stopifnot(all(sort(unique(unlist(iilist))) == sort(unlist(iilist))))
  return(iilist)
}




##' (From the "cruisedat" R package) Helper to get memberships. Assumes that the
##' times are in the format of "2017-06-13T17:46:23".
##'
##' @param times Vector of particularly formatted time strings.
##'
##' @return List of memberships (indices) of each of the entries in the |time|
##'   vector, divided into each unique hour.
##'
##' @export
get_membership <- function(times){
  ## All unique hours
  hours = sapply(times, function(mytime){
    substring(mytime, 1, 13)
  })
  all.hours = unique(hours)

  ## sapply(X$time[1:100], function(x){
  ##   lubridate::round_date(x, unit = "hour")
  ## })

  ## Cytogram membership in those unique hours
  membership = lapply(all.hours, function(one.hour){
    which(hours == one.hour)
  })

  ## Extract the hour strings
  strings = sapply(membership, function(a)names(a[1]))
  yyyy = as.numeric(sapply(strings, substr, 1, 4))
  mm = as.numeric(sapply(strings, substr, 6, 7))
  dd = as.numeric(sapply(strings, substr, 9, 10))
  hh = as.numeric(sapply(strings, substr, 12, 13))

  ## Make hourly time labels
  timemat = cbind(yyyy=yyyy, mm=mm, dd=dd, hh=hh)
  time = unlist(Map(function(y,m,d,h){paste0(y, "-0",
                                             m, "-",
                                             (if(d>=10) d else paste0(0, d)),
                                             "T",
                                             (if(h>=10) h else paste0(0, h)), ":00:00")}, yyyy, mm, dd,hh))

  return(list(membership = membership,
              time = time,
              timemat = timemat))
}




##' Create some artificial but realistic data.
##'
##' @param datadir Data directory
##' @param filename File name of the "cvres" object from the 2-76-5 experiment
##'
##' @return A list containing the generating coefficients, true means, and data
##'   (ylist, X, countslist=NULL for now).
generate_data_1d_pseudoreal_from_cv <- function(datadir, seed = NULL,
                                                nt = 100,
                                                ## Optionally binning the data
                                                bin = FALSE, dat.gridsize = 30){

  ## Load best 1d CV result
  ## cvres = blockcv_summary(2, 76, 5, 10, nrep = 5, datadir = datadir)##, subfolder="orig")
  ## saveRDS(cvres, file=file.path("~/repos/cruisedat/export", "1d-cvres.rds"))
  ## cvres = readRDS(file=file.path(datadir, filename)) #
  ## cvres = readRDS(file=file.path("~/repos/cruisedat/export", "1d-cvres.rds"))
  cvres = readRDS(file=file.path(datadir, "1d-cvres.rds"))
  ## Save this cvres and load it from datadir
  res = cvres$bestres
  X = res$X
  numclust = res$numclust
  TT = nrow(X)
  sigmas = (cvres$pretty.sigmas)

  ## Threshold beta at 0.01
  gen_beta = do.call(cbind, res$beta)
  small = which(abs(gen_beta[2:nrow(gen_beta),])<1E-2)
  gen_beta[2:nrow(gen_beta),][small] = 0

  ## Threshold alpha at 0.1
  gen_alpha = t(res$alpha)
  small = which(abs(gen_alpha[2:nrow(gen_alpha),])<1E-1)
  gen_alpha[2:nrow(gen_alpha),][small] = 0


  ## Beta coefficients
  mnmat = cbind(1, X) %*% gen_beta

  ## Alpha coefficients
  prob = exp(cbind(1,X) %*% gen_alpha)
  prob = prob/rowSums(prob)

  ## Particles per cytogram
  ntlist = rep(nt, TT)

  ## Samples |nt| memberships out of (1:numclust) according to the probs in prob.
  ## Data is a probabilistic mixture from these two means, over time.
  ylist = lapply(1:TT,
                 function(tt){
                   ## print(round(prob[[tt]],3))
                   draws = sample(1:numclust,
                                  size = ntlist[tt], replace = TRUE,
                                  ## prob = c(prob[[tt]], 1-prob[[tt]]))
                                  prob = c(prob[tt,]))
                   mns = mnmat[tt,]
                   means = mns[draws]
                   noises = sapply(draws, function(iclust){ stats::rnorm(1, 0, sigmas[iclust])})
                   datapoints = means + noises
                   cbind(datapoints)
                 })


  ## Make into countslist
  if(bin){
    dat.grid = flowmix::make_grid(ylist, gridsize = dat.gridsize) ## Having this to be common among all things is important
    obj = flowmix::bin_many_cytograms(ylist, dat.grid, mc.cores = 8, verbose = TRUE) ## This code needs to be made into 1d data
    ylist = lapply(obj$ybin_list, cbind)
    countslist = obj$counts_list
  } else {
    countslist = NULL
  }


  ## ## Other things about the true generating model (INCOMPLETE)
  ## sigma = array(NA, dim=c(2,1,1))
  ## sigma[1,1,1] = sigma[2,1,1] = 1
  ## mn = array(NA,dim=c(100,1,2))
  ## mn[,1,] = mnmat
  ## numclust=5


  return(list(ylist = ylist, X = X,
              countslist = countslist,
              mnmat = mnmat,
              prob = prob,
              alpha = gen_alpha,
              beta = gen_beta,
              numclust = numclust))
}




##' From simulation results, loads |bestreslist| containing the chosen models
##' from CV, then calculates the out-of-sample log likelihood of these models
##' with respect to a large new dataset (cytograms) generated from the same
##' mechanism.
##'
##' @param blocktype CV block type
##' @param datatypes 80:89
##' @param numclust 2
##' @param cv_gridsize 7
##' @param nrep 5
##' @param outputdir Output directory on client computer
##'
##' @return List of OOS scores.
outsample_loglik <- function(blocktype = 2,
                             datatypes = c(80:89),
                             numclust = 2,
                             cv_gridsize = 7,
                             nrep = 5,
                             datadir = "../paper-data"
                             ){

  ## Generate large data once
  datobj_large = generate_data_1d_pseudoreal(datadir = "~/repos/cruisedat/export",
                                             nt = 20000,
                                             beta_par = 0.3,
                                             p = 10)

  ## Calculate the out-of-sample scores
  ## datatypes = c(80:89)
  ## datatypes = c(80, 83, 86, 89)
  avg_oos = rep(NA, length(datatypes))
  quantile_oos = matrix(NA, nrow=length(datatypes), ncol = 4)
  nsim = 100
  all_oos = matrix(NA, nrow = length(datatypes), ncol = nsim)
  probs = c(0.05, 0.95, 0.25, 0.75)
  for(ii in 1:length(datatypes)){

    ## Load the best model from each noise level
    datatype = datatypes[ii]

    ## destin = file.path(outputdir, paste0("blockcv-", blocktype, "-", datatype, "-", numclust))
    ## load(file = file.path(destin, "summary", "bestreslist.Rdata"))

    ## destin = file.path(outputdir, paste0("blockcv-", blocktype, "-", datatype, "-", numclust))
    ## load(file = file.path(destin, "summary", "bestreslist.Rdata"))

    datadir = "../paper-data"
    load(file = file.path(datadir, "simulation",
                          paste0("bestreslist-noise-", datatype, ".Rdata")))

    nsim = length(bestreslist)
    for(isim in 1:nsim){
      ## lapply(1:9, function(isim){
      bestres = bestreslist[[isim]]
      if(is.null(bestres)) next
      class(bestres) = "flowmix"
      bestres$prob = bestres$pie
      oos_loglik = objective_newdat(ylist = datobj_large$ylist,
                                    countslist = datobj_large$countslist,
                                    res = bestres)
      all_oos[ii, isim] = oos_loglik
    }
  }

  ## Obtain the mean & quantiles
  avg_oos = apply(all_oos, 1, mean, na.rm = TRUE) ##mean(all_oos[ii,])
  quantile_oos = apply(all_oos, 1, function(myrow){
    quantile(myrow, probs = probs, na.rm=TRUE)
  })

  ## Return results
  return(list(all_oos = all_oos,
              avg_oos = avg_oos,
              quantile_oos = quantile_oos))
}





##' Plot 1 dimensional model results.
##'
##' @param ylist Data.
##' @param countslist Multiplicity.
##' @param res Optionally, |covarem| object of fitted model.
##' @param scale Defaults to TRUE. If TRUE, scale the data colors, relative to
##'   the largest count in |countslist|.
##' @param main Plot title.
##' @param time Vector of date strings.
##'
##'
##' @return NULL.
##' @export
plot_1d <- function(ylist,
                    countslist,
                    res = NULL,
                    scale = TRUE,
                    main = "",
                    time_axis = FALSE,
                    time = NULL,
                    mn.scale = 5,
                    cex_clust_label = 1.5,
                    omit_band = FALSE,
                    omit_label = FALSE,
                    reorder_clusters = TRUE,
                    cex_data = 1,##15
                    ylim = NULL,
                    cols = NULL){

  ## Plot ylist only
  stopifnot(ncol(ylist[[1]]) == 1)
  stopifnot(ncol(countslist[[1]]) == 1)
  plot_1d_ylist(ylist, countslist, scale = scale, main = main, cex = cex_data, ylim = ylim)
  if(!is.null(res)) assertthat::assert_that(methods::is(res, "flowmix"))

  ## Make the X ticks dates.
  if(time_axis){

    ## Find time.
    if(!is.null(time)){
      assertthat::assert_that(check_if_date(time))
    } else {
      if(!is.null(res)){
        check_if_date(rownames(res$X))
        assertthat::assert_that(check_if_date(rownames(res$X)))
        time = res$X %>% rownames() %>% lubridate::as_datetime()
      }
    }
    add_date_ticks_from_dates(time)
  } else {
    axis(1)
    axis(2)
  }

  ## Plot model means, if provided.
  if(!is.null(res)){
    stopifnot(res$dimdat==1)

    if(is.null(cols)){
    if(res$numclust <= 8){
      cols = RColorBrewer::brewer.pal(max(3,res$numclust), "Set2")
    } else {
      cols = RColorBrewer::brewer.pal(max(3,res$numclust), "Paired")
    }
    }

    ## Reorder clusters, according to total prob.
    if(reorder_clusters) res <- reorder_clust(res)

    ## Plot the means
    add_mn(res, cols, mn.scale)

    ## Add confidence band
    if(!omit_band){  add_band(res, cols)   }

    ## Make the remaining region outside of TT white, to prevent spillover of
    ## the cluster means' cex=17 points
    fill_both_sides(length(ylist))

    ## Add Cluster labeling
    if(!omit_label){
      add_cluster_label(res, cex_clust_label)
    }
  }
  return(list(res = res, cols = cols))
}




##' Make 1d data plot. No axes are added.
##'
##' @inheritParams plot_1d
##' @param cex point size.
##' @param ylim y limits for plot.
##'
##' @export
plot_1d_ylist <- function(ylist, countslist, res = NULL, scale = TRUE, main = "",
                       cex = 1, ylim = NULL){

  stopifnot(ncol(ylist[[1]]) == 1)
  TT = length(ylist)
  if(is.null(ylim)) ylim = range(unlist(ylist))
  graphics::matplot(NA,
          type = 'l',
          lty = 1,
          lwd = .1,
          ylim = ylim,
          xlim = c(1, TT),
          ylab = "",
          xlab = "",
          axes = FALSE)

  ## Add title
  title(main = main)

  ## Handle when countslist is NULL
  if(!is.null(countslist)){
    if(!scale) mx = 10 else mx = max(unlist(countslist))
    cols = lapply(countslist, function(counts){
      ct = counts / mx
      return(rgb(0, 0, 0, ct))
    })
    pch = 15
  } else {
    countslist = sapply(ylist, function(y) rep(1, nrow(y)))
    cols = lapply(1:TT, function(a) rgb(0,0,0,0.1))
    pch = 16
  }

  ## Visualize the data
  for(tt in 1:TT){
    y = ylist[[tt]]
    col = cols[[tt]]
    points(x = rep(tt, length(y)),
           y = y, col = col,
           pch = pch, cex = cex)
  }
}



##' Add mean (only for 1d data)
##' @noRd
add_mn <- function(res, cols, mn.scale = 5){
  for(iclust in 1:res$numclust){
    graphics::lines(res$mn[,1,iclust], type = 'o', lty = 1, lwd = .1,
                    cex = res$prob[,iclust] * mn.scale, pch = 15, col = cols[iclust])
  }
}

##' Add bands (only for 1d data).
##' @noRd
add_band <- function(res, cols){
  sigmas = sqrt(res$sigma[,1,])
  sigma_mat = matrix(sigmas,## res$sigma[,1,],
                     ncol=res$numclust,
                     nrow=res$TT, byrow=TRUE)
  for(iclust in 1:res$numclust){
    up = res$mn[,1,iclust] + 2*sigma_mat[,iclust]
    dn = res$mn[,1,iclust] - 2*sigma_mat[,iclust]
    graphics::polygon(c(1:length(up),
                        rev(1:length(up))),
                      c(up,rev(dn)),
                      col = grDevices::adjustcolor(cols[iclust], alpha.f = 0.2),
                      border = NA)
  }
}





##' Make the remaining region outside of TT white, to prevent spillover of the
##' cluster means' cex=17 points. (only for 1d data).
##' @param TT Integer.
##' @noRd
fill_both_sides <- function(TT){
  x = c(TT -1 + (1:100), TT -1 + (100:1))
  polygon(x = x, y = c(rep(-5, 100), rep(5, 100)),
          col = "white",
          border=FALSE)

  x = c(3-(1:100), 3-(100:1))
  polygon(x = x, y = c(rep(-5, 100), rep(5, 100)),
          col = "white",
          border=FALSE)
}




##' Check if time is TRUE.
##'
##' @param date String vector containing date strings
##'   e.g. "2017-06-13T11:00:00".
##'
##' @return TRUE only all are dates..
check_if_date <- function(date){
  date = date %>% lubridate::as_datetime()
  if(length(date) == 0) return(FALSE)
  return(all(sapply(date, lubridate::is.POSIXt)))
}





##' Add date ticks from string of dates.
##'
##' @param dates Vector of strings of the form "2017-05-29T00:00:00", or
##'   otherwise recognizeable using \code{lubridate::as_datetime()}.
##' @param empty_tick_marks Defaults to FALSE. If TRUE, use EMPTY tick marks.
##' @param ... Rest of arguments to both axes via \code{axis()}.
##'
##' @return No return.
##' @export
add_date_ticks_from_dates <- function(dates, empty_tick_marks=FALSE, ...){
  ## dates = sapply(as.Date(dates) %>% format("%B %d"), toString)

  ## Get dates and X coordinates
  dates = sapply(lubridate::as_datetime(dates) %>% format("%B %d"), toString)
  nums = as.numeric(as.factor(dates))

  ## Form the tick locations.
  left_ticks = sapply(sort(unique(nums)),function(ii){min(which(nums==ii))})
  left_ticks = c(left_ticks, length(dates))##res$TT)
  mid_ticks = sapply(sort(unique(nums)),function(ii){mean(which(nums==ii))})
  dates_mid_ticks = dates[round(mid_ticks)]
  if(empty_tick_marks) dates_mid_ticks = rep("", length(dates_mid_ticks))

  ## Place those ticks
  graphics::axis(1, at = left_ticks, labels = FALSE)
  graphics::axis(1, at = mid_ticks, labels = dates_mid_ticks, tick = FALSE, las=2, ...)
  graphics::axis(2, ...)
}





##' Add cluster label near beginning of plot (only for 1d data).
##' @param cex if NULL, don't make cluster label.
##' @param res flowmix object.
##' @export
add_cluster_label <- function(res, cex=1.5){
  if(!is.null(cex)){
    for(iclust in 1:res$numclust){
      y = res$mn[1, 1, iclust]
      graphics::text(x = res$TT/50, y = y,
                     label = paste0("Cluster ", iclust), cex = cex)
    }
  }
}


