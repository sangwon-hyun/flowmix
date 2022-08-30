##' Make a 1d binned data plot, using ggplot.
##'
##' @param ylist Must be numeric vectors containing bin centers.
##' @param qclist Corresponding QC in those bins.
##' @param col Vector of colors.
##'
##' @return ggplot object.
##'
##' @export
##' @import ggplot2
##' @importFrom lubridate as_datetime
bin_plot_1d <- function(ylist, qclist, col=NULL){

  time = NULL ## fixing check()

  if(is.null(col)) col = c("white", "black", "yellow", "red")

  stopifnot(all(sapply(ylist, length) == sapply(qclist, length)))

  combined_dat = lapply(names(ylist), function(onedatetime){
    tibble(y = ylist[[onedatetime]],
           qc = qclist[[onedatetime]],
           time = onedatetime)
  }) %>% bind_rows()
  ## combined_dat = combined_dat %>% mutate(time = lubridate::as_datetime(time))
  combined_dat = combined_dat %>% mutate(time = as_datetime(time))

  p = combined_dat %>% ggplot() + geom_raster(aes_string(x = 'time', y = 'y', fill = 'qc')) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_fill_gradientn(colours = col)
  return(p)
}
