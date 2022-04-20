##' Makes completely white (Gaussian) noise.
make_toy_dat1 <- function(type = c("raw", "bin")){

  ## Setup
  type = match.arg(type)

  ## Generate some data
  datlist = lapply(1:300, function(tt) rnorm(1000))##rnorm(10000))
  if(type == "raw"){ return(datlist) }
  if(type == "bin"){
    rng = range(unlist(datlist))
    breaks = seq(from = rng[1], to = rng[2], length = 30)
    allcounts = datlist %>%
      purrr::map(. %>% hist(, breaks = breaks, plot = FALSE) %>% purrr::pluck("counts") %>% t() %>% as_tibble()) %>%
      data.table::rbindlist() %>% as_tibble()
    colnames(allcounts) = breaks
    ## allcounts = allcounts %>% mutate_all(function(a){pmin(a, .01)})
    return(allcounts)
  }
}



##' From model and ylist & countslist, get deviance residuals.
get_devres <- function(res, ylist, countslist){
  mu = res$mn
  prob = res$prob
  sigma = res$sigma
  devres = objective(mu, prob, sigma, ylist, prob_lambda=0, mean_lambda=0, countslist = countslist, each = TRUE)
}



##' Given model and particle data, make residuals
get_resid <- function(res, dat){
  gateres = gate(res, dat$ylist, dat$countslist)
  TT = length(dat$ylist)
  resids = lapply(1:TT, function(tt){
    means_tt = res$mn[tt,1,]
    mem_tt = gateres[[tt]]
    particle_means = means_tt[mem_tt]
    particle = dat$ylist[[tt]]
    residuals = particle - particle_means
    return(tibble(x=tt, y=as.numeric(residuals)))
  }) %>% bind_rows()
}



##' Get scaled residuals
get_scaled_resid <- function(res, dat){
  gateres = gate(res, dat$ylist, dat$countslist)
  TT = length(dat$ylist)
  sigmas = sqrt(res$sigma[,1,1])
  resids = lapply(1:TT, function(tt){
    means_tt = res$mn[tt,1,]
    mem_tt = gateres[[tt]]
    particle_means = means_tt[mem_tt]
    particle_sigmas = sigmas[mem_tt]
    particle = dat$ylist[[tt]]

    ## Calculate residuals
    residuals = particle - particle_means
    scaled_residuals = residuals / particle_sigmas

    ## plot(residuals, col = mem_tt)
    ## plot(residuals/particle_sigmas, col = mem_tt)

    return(tibble(x=tt, y=as.numeric(scaled_residuals)))
  }) %>% bind_rows()
}



##' Plots |vec| vs |lagged vec|, for up to 25 lags.
plot_lagged <- function(vec, numlag = 18){

  lagged_vecs = lapply(1:numlag, function(kk){
    tibble(x = vec,
           y = lag(vec, kk),
           lag = paste0(kk))
  }) %>%
    bind_rows() %>%
    mutate(lag = factor(lag, levels = c(1:numlag)))


my.formula <- y ~ x
lagged_vecs %>% ggplot(aes(x=x, y=y)) +
  facet_wrap(~lag, labeller = label_both) +
  geom_smooth(method="lm") +
  geom_point() +
  stat_poly_eq(formula = my.formula,
                aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               parse = TRUE,
               col = 'blue')  +
  xlab("Lagged") +
  ylab("Original") +
  coord_fixed()
}



##' Helper function to get scaled residuals from real data's residuals.
get_scaled_resids_realdata_3d <- function(tt, residuals_obj, sigma){

  ## Setup
  numclust = length(residuals_obj$residuals_by_cluster[[1]])

  ## 1. Get the residuals and counts from each cluster
  resids = lapply(1:numclust, function(iclust){
    residuals_obj$residuals_by_cluster %>% .[[tt]] %>% .[[iclust]]
  })
  counts = lapply(1:numclust, function(iclust){
    residuals_obj$countslist_by_cluster %>% .[[tt]] %>% .[[iclust]]
  })

  ## 2. Scale residuals by the cluster's standard deviation
  sigma_sqrt_inv_list = lapply(1:numclust, function(iclust){
    onesigma = sigma %>% .[iclust,,]
    one_svd = onesigma %>% svd()
    onesigma_inv_sqrt = one_svd$u %*% diag(1/sqrt(one_svd$d)) %*% t(one_svd$v)
  })
  ## solve(onesigma_sqrt) %*% resids

  scaled_resids = Map(function(r, s){ r %*% s },
                      resids, sigma_sqrt_inv_list)

  ## 3. Combine and return them
  return(list(resids = scaled_resids, counts = counts))
}



##' Helper function to get scaled residuals from real data's residuals.
get_scaled_resids_realdata_1d <- function(tt, residuals_obj, sd_list){

  ## Setup
  numclust = length(residuals_obj$residuals_by_cluster[[1]])

  ## 1. Get the residuals and counts from each cluster
  resids = lapply(1:numclust, function(iclust){
    residuals_obj$residuals_by_cluster %>% .[[tt]] %>% .[[iclust]]
  })
  counts = lapply(1:numclust, function(iclust){
    residuals_obj$countslist_by_cluster %>% .[[tt]] %>% .[[iclust]]
  })

  ## 2. Scale residuals by the cluster's standard deviation
  scaled_resids = Map(function(r, s){ r / s },
                      resids, sd_list)

  ## 3. Combine and return them
  return(list(resids = scaled_resids, counts = counts))
}


get_resids_realdata_1d <- function(tt, residuals_obj){

  ## Setup
  numclust = length(residuals_obj$residuals_by_cluster[[1]])

  ## Get the residuals and counts from each cluster
  resids = lapply(1:numclust, function(iclust){
    residuals_obj$residuals_by_cluster %>% .[[tt]] %>% .[[iclust]]
  })
  counts = lapply(1:numclust, function(iclust){
    residuals_obj$countslist_by_cluster %>% .[[tt]] %>% .[[iclust]]
  })

  return(list(resids = resids, counts = counts))
}




##' Just bundling the analysis up.
##'
##' @param mat_oneclust This is one cluster's worth of residuals, in long format
##'   (columns are "time", "y", and "iclust").
##'
pipeline <- function(mat_oneclust, type = "median", truncate_distribution = FALSE){

  ## Basic checks
  stopifnot(setequal(colnames(mat_oneclust), c("time", "y", "iclust", "counts")))

  ## If we should truncate, let's do that..
  if(truncate_distribution){
    mat_oneclust = mat_oneclust %>% dplyr::filter(abs(y) < 2)
  }

  ## Apply two tools
  ## 1. Ljung-box on means of residuals.
  ## 2. (NOT DOING THIS ) Ljung-box on *deviance* residuals.
  ## 3. Binary Expansion Test (BET) on residuals.

  tool1(mat_oneclust, type)
  ## tool3(mat_oneclust)
}

## Ljung-box on medians of residuals
tool1 <- function(mat_oneclust, type = "median"){

  ## Calculate medians from a small sample.
  scaled_resids = mat_oneclust %>%
    .[sample(1:nrow(mat_oneclust), min(nrow(mat_oneclust), 1E6)),] %>%
    select(x = time, y, counts)
  if(type == "median") fn = weighted_median
  if(type == "mean") fn = weighted.mean
  medians = scaled_resids %>% group_by(x) %>% summarise(y = fn(y, counts))
  ## medians = resids %>% group_by(x) %>% summarise(y = weighted_median(y, counts))

  #########################
  #### Make five plots ####
  #########################

  p0 = scaled_resids %>%
    ggplot(aes(x = x, y = y, size= counts)) +
    geom_point(col = rgb(0, 0, 0, 0.1)) +
    scale_size_area() +
    scale_size_continuous(range = c(0,3))
    ylab("Scaled residuals")  +
    theme_bw()
  plot(p0)

  ## Plot the binned data
  p1 = scaled_resids %>% ggplot(aes(x=x, y=y, weight = counts)) +
    geom_hex(bins = 100) +
    scale_fill_continuous(type = "viridis") +
    ## scale_colour_manual(c("white", "black", "yellow")) +
    ## scale_fill_continuous(col=c("white", "black")) +
    ## scale_fill_continuous(type = "gradient") +
    ylab("Scaled residuals")  +
    theme_bw()
  plot(p1)

  ## Series of histograms at spaced out points
  p1_by_time =
    scaled_resids %>% dplyr::filter(x %% 10==0) %>%
    rename(time=x) %>%
    ggplot(aes(y, weight = counts)) +
    facet_wrap(~time, scale = "free_y", labeller = label_both) +
    geom_histogram(aes(y=..density..)) + ##, bins=20)
    geom_density(col='red', bw=.5) +
    theme_bw() +
    ggtitle("Histograms at every 10'th time point")
  print(p1_by_time)

  ## Plot the medians
  if(type == "mean") ylab = "Means of residuals"
  if(type == "median") ylab = "Medians of residuals"
  p2 = medians %>% ggplot(aes(x=x, y=y)) +
    geom_line() +
    geom_point() +
    ylab(ylab) + xlab("Time")
  plot(p2)

  ## Ljung-box on Medians
  TT = nrow(medians)
  mq_table(medians, round(sqrt(TT)))

  ## Autocorrelation plots
  acf(medians$y, main = paste0("ACF of ", ylab))
  plot_lagged(medians$y) %>% plot()

  ## ## BET plots
  ## library(BET)
  ## set.seed(82443684)
  ## resids = scaled_resids %>% .[sample(1:nrow(scaled_resids), 1E3),] %>% select(x = x, y)
  ## betres = BEAST(as.matrix(resids), 3, test.independence = TRUE)
  ## ## print(kable(betres))
  ## print("BET test p-value is")
  ## print(betres[["p-value"]])
  ## bet.plot(X=as.matrix(resids), d=3)
}


## ##' Apply BET test on residuals form one cluster.
## tools3 <- function(mat_oneclust){
##   library(BET)
##   set.seed(82443684)
##   resids = mat_oneclust %>% .[sample(1:nrow(mat_oneclust), 1E3),] %>% select(x = time, y)
##   betres = BEAST(as.matrix(resids), 3, test.independence = TRUE)
##   print(betres)
##   bet.plot(X=as.matrix(resids), d=3)
## }


mq_table <- function(x, lag, caption = NULL){
  capture.output(mq_invisible(x, lag = lag))  %>%
    strsplit("\\s+") %>% purrr::map(. %>% strsplit("\\s+")) %>%
    purrr::map(.%>%tail(4) %>% unlist()) %>% do.call(rbind,.) %>%
    knitr::kable(caption = caption) %>% print()
  cat(fill=TRUE)
}



## Weighted median function from here: https://stackoverflow.com/questions/2748725/is-there-a-weighted-median-function
weighted_median <- function(x, w) {
  w <- w[order(x)]
  x <- x[order(x)]

  prob <- cumsum(w)/sum(w)
  ps <- which(abs(prob - .5) == min(abs(prob - .5)))
  return(x[ps])
}




##' Use when nested list is in time->cluster order, and you'd like it to be
##' cluster->time order.
cluster_then_time <- function(residuals){
  resids_by_clust = list()
  TT = length(residuals)
  numclust = length(residuals[[1]])
  for(iclust in 1:numclust){
    resids_thistime = list()
    for(tt in 1:TT){
      resids_thistime[[tt]] = residuals[[tt]][[iclust]]
    }
    resids_by_clust[[iclust]]  = resids_thistime
  }
  return(resids_by_clust)
}


mq_invisible <- function(...){
    ff <- tempfile()
    png(filename=ff)
    res <- mq(...)
    dev.off()
    unlink(ff)
    res
}




##' Performs the permuation test.
##' @param means A 3d matrix of means.
permtest <- function(means, nperm = 100, plot = FALSE, title=NULL){

  test_stats = sapply(1:nperm, function(iperm){
    rownums = sample(1:nrow(means), replace = FALSE)
    permuted_means = means[rownums,,drop=FALSE]
    test_stat(permuted_means)
  })
  orig_stat = test_stat(means)
  if(plot){
    hist(test_stats, xlim = c(-1, 1),
         xlab = "Test statistic under\nPermutation null", main = title, border = "grey")
    abline(v = orig_stat, col = 'red', lwd=3)
    ## return(p)
  }

  ## How often is are the null samples larger than what's observed?
  numer = (sum(test_stats > abs(orig_stat)) + sum(test_stats < -abs(orig_stat)))
  pv = numer / nperm
  return(pv)
}

##' Calculate the test statistic (corr(y, lag(y,1))), given a matrix whose rows
##' are histogram density values.
##'
##' @param means A 3d matrix of means.
##'
test_stat <- function(means){
  threestats = sapply(1:ncol(means), function(idim){
    mns = means[,idim, drop=TRUE]
    mnmat = cbind(orig = mns, lagged = lag(mns,1)) %>% as_tibble() %>%  na.omit()
    onestat = cor(mnmat %>% pull(orig), mnmat %>% pull(lagged))
    return(onestat)
  })
  max(threestats)
}







############ TEMPORARY ###############################################
############ TEMPORARY ###############################################
############ TEMPORARY ###############################################
############ TEMPORARY ###############################################
############ TEMPORARY ###############################################
############ TEMPORARY ###############################################

##' @param bin1 bin centers of 1 deconvolved distribution.
##' @param count1 counts in each bins.
##' @param bin2 counterpart to bin1
##' @param count2 counterpart to count1
emd <- function(bin1, count1, bin2, count2){


  ## Make cost matrix
  costm = make_cost_matrix(bin1, bin2, verbose = FALSE)

  ## Do the transport calculation
  count1 = count1/sum(count1)
  count2 = count2/sum(count2)
  res = transport::transport(count1, count2, costm = costm, p = 2)$dist
  pp = 2
  dist = transport::wasserstein(count1,
                                count2,
                                costm = costm^pp,
                                tplan = res)

  return(dist^(1/pp))
}


##' @param bin1 3d matrix
##' @param bin2 3d matrix
##' @param verbose TRUE if you want to see stuff
make_cost_matrix <- function(bin1, bin2, verbose = FALSE){
  n1 = nrow(bin1)
  n2 = nrow(bin2)
  dmat = matrix(NA, n1, n2)
  for(i1 in 1:n1){
    if(verbose) printprogress(i1, n1)
    for(i2 in 1:n2){
      dmat[i1, i2] = sqrt(sum((bin1[i1,] - bin2[i2,])^2))
    }
  }
  return(dmat)
}




##' Taken from 01-helpers in ~/repos/omd/main/.
make_mds_plot_temp <- function(distmat,
                          shortnames = c(paste0("real-", 1:12), paste0("darwin-", 1:12)),
                          longnames  = c(rep("Remote sensing", 12), rep("Darwin model", 12)),
                          angle = -pi * 0.2,
                          path = FALSE
                          ){

  ## Basic checks

  diag(distmat) = NA
  colnames(distmat) = rownames(distmat) = shortnames

  ## Turn into "dist" class object
  distmat <- as.dist(distmat, diag = TRUE)

  ## Calculating an MDS plot
  fit <- cmdscale(distmat, eig=TRUE, k=2) # k is the number of dim

  # Plot solution
  x <- fit$points[,1]
  y <- fit$points[,2]

  ## Define rotation
  get_rotation_matrix = function(aa){
    cbind(c(cos(aa), sin(aa)),
        c(-sin(aa), cos(aa)))
  }
  rot = get_rotation_matrix(angle)

  ## Plot the points on a 2d plot
  dt = data.frame(x = x, y = y)
  dt = as.matrix(dt) %*% rot
  dt = data.frame(x= dt[,1], y = dt[,2],
                  labels = 1:nrow(dt),
                  dat_type = longnames) %>%
    arrange(labels)

  ## some_shortnames = shortnames[c(12, 1, 24, 13)]
  ## ## some_shortnames = c("real-12", "real-1", "darwin-12", "darwin-1"),])
  ## dt_extra = rbind(dt[some_shortnames,])

  p =
    dt %>% ggplot(aes(x = x, y = y, label = labels, col = dat_type)) +
    geom_point(cex = 2) +
    ## geom_path(aes(x=x, y=y),  data =dt) +
    ## geom_path(aes(x = x, y = y, label = labels, col = dat_type), alpha=.5, lty=2, data = dt_extra) +
    ## geom_text_repel(cex = rel(5), ##fontface = "bold",
    ##                  alpha = .8, show.legend = F) +
                     ## label.padding=.05) +
    labs(col = "Data source") +
    theme(legend.title= element_blank()) +
    theme_minimal() +
    xlab("Coordinate 1") +
    ylab("Coordinate 2")
  if(path) p = p + geom_path(aes(x=x, y=y),  data =dt)
  ##   ## ggtitle("Multidimensional scaling of climatology data (1998-2018)") +
  ##   scale_color_manual(values = c("black","blue")) +
  ##   ## scale_color_manual(values=c("black","blue")) +
  ##   coord_fixed()  +
  ##   theme(legend.text = element_text(size = rel(1), face = 1))
  ## p = p + theme(legend.justification = c(1, 0), legend.position = c(1, 0))
  return(p)
}
