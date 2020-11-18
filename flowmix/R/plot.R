#######################
## Function generics ##
#######################

plot2d <- function (x, ...) {
   UseMethod("plot2d", x)
}
plot2d_generic <- function (x, ...) {
   UseMethod("plot2d_generic", x)
}
