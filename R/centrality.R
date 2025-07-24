#' @title Spatial centrality
#' @author Gabriel Gaona
#' @description Functions to find spatial measures of gravity centers.
#' @note Inspired on `aspace::*()` from Ron Buliung & Randy Bui (2012)
#' @param .x,.y  \code{sf} points 2D or 3D
#' @param weights Numeric. Used in for weigthed Mean Center. Has to be same length
#'                as number of points.
#' @param method Character. Type of center point to calculate
#' @param dist Atomic numeric, Default 100. Starting distance value for center
#'             moving during iterations.
#' @param ... arguments to be passed to or from other methods
#' @details Spatial centers are spatial measures of the gravity center.
#'
#' `methods` options are:
#'  `"mean"` is the mean center (equivalent to centroid of the points) calculated by
#'  the arithmetic mean of each axis;
#'  `"geometric"`, is the corresponding geometric mean of each axis;
#'  `"median"`, is the median center, a pair of c(median(x), median(y)) coordinates;
#'  `"feature"`, is a minimization of the sum of distances from _ith_ point to every point;
#'  `"min.dist"`, is iterative looking for the closest point in bbox of .x
#'  that minimizes the sum of distances from _ith_ point to every point in the dataset
#'
#' @return \code{"Simple Features"} of lenght 1.
#' @importFrom stats median sd
#' @examples
#' requireNamespace("ggplot2", quietly = TRUE)
#' library(sf, quietly = TRUE)
#' library(ggplot2)
#' bbx <- matrix(c(697047,9553483,
#'                 696158,9560476,
#'                 700964,9561425,
#'                 701745,9555358),
#'               byrow = TRUE,
#'               ncol = 2)
#' bbx <- st_multipoint(bbx)
#' bbx <- st_cast(bbx,"POLYGON")
#' bbx <- st_sfc(bbx, crs = 31992)
#' set.seed(1234)
#' points <- st_sf(geometry = st_sample(bbx, 100))
#' mean_center <- st_central_point(points, method = "mean")
#' median_center <- st_central_point(points, method = "median")
#' geom_center <- st_central_point(points, method = "geometric")
#' central_feature <- st_central_point(points, method = "feature")
#' min_dist_center <- st_central_point(points, method = "min.dist")
#' ggplot() +
#'   geom_sf(data = points, color = "steelblue", size = 0.5) +
#'   geom_sf(data = mean_center, color = "blue", size = 3) +
#'   geom_sf(data = median_center, color = "red") +
#'   geom_sf(data = geom_center, color = "grey80") +
#'   geom_sf(data = central_feature, color = "orange") +
#'   geom_sf(data = min_dist_center, color = "green")

#' @export
#' @rdname centrality
st_central_point = function(.x, .y, ...) UseMethod("st_central_point")

#' @export
#' @rdname centrality
st_central_point.sfg <- function(.x, .y = NULL,
                                 weights = NULL,
                                 method = c("mean",
                                            "median",
                                            "geometric",
                                            "feature",
                                            "min.dist"), ...){
  .x <- st_geometry(.x)
  if(!is.null(.y)) .y <- st_geometry(.y)
  st_central_point.sfc(.x, .y, weights, method, ...)
}

#' @export
#' @rdname centrality
st_central_point.sf <- function(.x, .y = NULL,
                                weights = NULL,
                                method = c("mean",
                                           "median",
                                           "geometric",
                                           "feature",
                                           "min.dist"), ...){
  .x <- st_geometry(.x)
  if(!is.null(.y)) .y <- st_geometry(.y)
  st_central_point.sfc(.x, .y, weights, method, ...)
}

#' @export
#' @rdname centrality
st_central_point.sfc <- function(.x, .y = NULL,
                                     weights = NULL,
                                     method = c("mean",
                                                "median",
                                                "geometric",
                                                "feature",
                                                "min.dist"), ...){

  wstring <- "weighted "
  if(is.null(weights)) {
    wstring <- ""
    weights <- rep(1, length(.x))
  }

  method <- match.arg(method)

  if(method == "median")
    cp <- .median_centre(.x)

  if(method == "mean")
    cp <- .mean_centre(.x, weights = weights, ...)

  if(method == "geometric")
    cp <- .geom_mean_centre(.x, weights = weights, ...)

  if(method == "feature")
    cp <- .central_feature(.x)

  if(method == "min.dist")
    cp <- .centre_min_dis(.x, ...)

  # if(method == "feature.2sets")
  #   if(is.null(.y)) stop("Required .y dataset for method feature.2sets")
  #   cp <- .central_feature2(.x, .y, ...)

  return(
    st_sf(feature = paste0(wstring, method),
          geometry = cp)
  )
}

#' @param .x [sf] `geometry`
#' @param ... additional parameters for [`stats::median()`]
#' @noRd
.median_centre <- function(.x, ...) {
  crd <- apply(st_coordinates(.x), 2, median, ...)
  st_sfc(st_point(crd), crs = st_crs(.x))
}

#' @param .x [sf] `geometry`
#' @param ... additional parameters for [`stats::weighted.mean()`]
#' @noRd
.mean_centre <- function(.x, weights = NULL, ...) {
  if(is.null(weights)) weights <- rep(1, length(.x))
  crd <- apply(st_coordinates(.x), 2, weighted.mean, w = weights, ...)
  st_sfc(st_point(crd), crs = st_crs(.x))
}

#' @param .x [sf] `geometry`
#' @noRd
.geom_mean_centre <- function(.x, weights = NULL, ...) {
  if(is.null(weights)) weights <- rep(1, length(.x))

  .f <- function(x, w){
    from <- range(x)
    to <- 0:1
    crd <- prod((w * scales::rescale(x, to))^(1/length(x)))
    scales::rescale(crd, from)
  }
  crd <- apply(st_coordinates(.x), 2, .f, w = weights, ...)
  st_sfc(st_point(crd), crs = st_crs(.x))
  }

#' @param .x [sf] `geometry`
#' @param ... additional parameters for [`sf::st_distance()`]
#' @noRd
.central_feature <- function(.x, ...) {
  d <- st_distance(.x, ...)
  ds <- rowSums(d)
  .x[which.min(ds)]
}


#' @param .x [sf] `geometry`
#' @param ... ignored
#' @noRd
.centre_min_dis <- function(.x, dist = 100, tol = 0.1, ...) {
  crs <- st_crs(.x)
  units(dist) <- st_crs(.x)$units
  units(tol) <- st_crs(.x)$units
  d <- dist * 2
  xy <- st_sfc(.mean_centre(.x), crs = crs)

  while (d >= dist & dist > tol) {
    grid <- .point_move(xy, dist = dist, crs = crs)
    sumdist.cmd <- colSums(st_distance(.x, grid))
    CMD <- grid[which.min(sumdist.cmd)]
    d <- st_distance(CMD, xy, by_element = TRUE)
    if(is.null(st_crs(.x)$units)){
      units(dist) <- units(d)
      units(tol) <- units(d)
    }

    xy <- CMD
    if(dist > tol){
      dist <- dist/2
    }
  }
  xy
}


#' @param .x [sf] `geometry`
#' #' @param ... additional parameters for [`sf::st_distance()`]
#' @noRd
.central_feature2 <- function(.x, .y, ...) {
  d <- colSums(st_distance(.x, .y))
  .x[which.min(d)]
}

#' @param .x [sf] `geometry`
#' @noRd
.point_move <- function(.x, dist = dist, crs = crs) {
  st_cast(st_sfc(st_buffer(.x, dist = dist), crs = crs), "POINT")
}

