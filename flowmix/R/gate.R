##' Gate (classify) the particles in \ref{ylist_particle} according to random
##' draws from the particle responsibilities -- from the model in \code{res}.
##'
##' @param res flowmix class object.
##' @param ylist_particle Particles
##' @param countslist_particle Counts or biomass of the particles.
##' @param seed Optionally set a seed for drawing particle membership.
##'
##' @return Memberships.
##'
##' @export
gate <- function(res, ylist_particle, countslist_particle, seed = NULL){

  ## Setup
  TT = length(ylist_particle)
  times = names(ylist_particle)

  ## Conduct the E-step once to calculate responsibilities
  resp <- Estep(res$mn, res$sigma, res$prob, ylist = ylist_particle,
                numclust = res$numclust, first_iter = TRUE)

  ## Draw cluster membership.
  if(!is.null(seed)) set.seed(seed)
  drawslist = draw_membership(resp)

  ## Calculate memberships (a vector consisting of (1,.., res$numclust))
  start.time = Sys.time()
  memberships = list()
  for(tt in 1:TT){
    memberships[[tt]] = drawslist[[tt]] %>% apply(1, function(rr) which(rr == 1))
  }
  cat(fill = TRUE)
  return(memberships)
}


##' Using gated particles (from \code{gate()}) , create one plot of the top
##' levels.
##'
##' @param tt Time point (one of \code{1:length(y_gated)}).
##' @param dims Which dimensions to plot.
##' @param y Particle data.
##' @param counts Particle multiplicities.
##' @param memberships Particle memberships.
##' @param datetime An object (string or datetime) containing the time and date,
##'   to go in the title.
##' @param top_clusters Memberships (cluster numbers) you'd like to highlight.
##' @param top_cluster_names (optional) Labels for the clusters in
##'   \code{top_clusters}.
##'
##' @return One ggplot object
gate_plot <- function(tt, dims, y, counts, memberships, datetime = "",
                               top_clusters, top_cluster_names = NULL){

  ## Setup
  dimdat = ncol(y)
  varnames = colnames(y)[dims]

  ## Set up cluster labels
  top_clusters = top_clusters %>% as.character()
  my_labeller = NULL
  if(!is.null(top_cluster_names)){
    names(top_cluster_names) = top_clusters
    my_labeller = labeller(mem = top_cluster_names)
  }

  ## Some aesthetics
  mygrey = rgb(0,0,0, 0.05) ## Plotting the background points ##"grey80"##
  myred = rgb(1,0,0, 0.2)

  ## Basic checks
  stopifnot(nrow(y) == length(counts))
  stopifnot(length(counts) == length(memberships))

  ## Gather the data
  y_aug = y %>% as_tibble() %>%
    dplyr::select(!!!varnames) %>%
    add_column(counts = counts) %>%
    add_column(mem = as.factor(memberships))


  ## Create ggplot
  y_aug_subset = y_aug %>% subset(mem %in% top_clusters)
  g = y_aug_subset  %>%
    group_by(mem) %>%
    ggplot() +
    geom_point(aes_(x = as.name(varnames[1]), y = as.name(varnames[2]), cex = ~sqrt(counts)), colour = mygrey, data = y_aug %>% select(-mem)) +
    geom_point(aes_(x = as.name(varnames[1]), y = as.name(varnames[2]), cex = ~sqrt(counts)), colour = myred) +
    facet_wrap(~mem, labeller = my_labeller) +

    ## Plot formatting
    xlim(c(0,8)) + ylim(c(0,8)) +
    ggtitle(datetime) +
    theme_bw() +
    theme(legend.position = "none") +
    theme(strip.text = element_text(colour = "black", face = "bold",  size = rel(1.5)),
          plot.title = element_text(size = rel(1.5), face = "bold")) +
    theme(aspect.ratio = 1) +
    theme(axis.title.y = element_text(size = rel(2)),
          axis.title.x = element_text(size = rel(2)),
          axis.text.x = element_text(face="bold", size=rel(1.5)),
          axis.text.y = element_text(face="bold", size=rel(1.5)))

  return(g)
}
