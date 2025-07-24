#' @title Standard deviation distance calculator
#' @description Calculate the spatial deviaction distance from a points sf dataset.
#' @author Gabriel Gaona
#' @param .x  \code{sf} points 2D or 3D
#' @param centre  One central point of class _sf_, _sfc_, _numeric_
#'   (length 2), _matrix_ (2 col, 1 row), _data.frame_  (2 col, 1 row),
#'   or _list_ (length 2). Default `NULL`, means a calculation of the [`st_central_point()`]
#'   from `.x` localities.
#' @param weights Numeric. Same length as number of points in `.x`.
#' @param ... other parameters for [`sf::st_distance()`]
#' @return A sf `"POLYGON"` with atributes:
#'
#' - `radius` (standard deviation distance)
#' - `area` surrounding,
#' - `perimeter`,
#' - `center` coordinates,
#' - `weigted` indicator if weights were used or not in the calculaton.
#' @examples
#'   requireNamespace("ggplot2", quietly = TRUE)
#'   library(sf, quietly = TRUE)
#'   library(ggplot2)
#'   bbx <- matrix(c(697047,9553483,
#'                   696158,9560476,
#'                   700964,9561425,
#'                   701745,9555358),
#'                 byrow = TRUE,
#'                 ncol = 2)
#'   bbx <- st_multipoint(bbx)
#'   bbx <- st_cast(bbx,"POLYGON")
#'   bbx <- st_sfc(bbx, crs = 31992)
#'   set.seed(1234)
#'   points <- st_sf(geometry = st_sample(bbx, 100))
#'   SDD <- st_sd_distance(points)
#'   ggplot() +
#'     geom_sf(data = SDD, fill = NA, color = "darkolivegreen") +
#'     geom_sf(data = points, color = "steelblue", size = 0.5)
#' @export
#' @rdname sd_distance
st_sd_distance = function(.x, centre = NULL, weights = NULL, ...) UseMethod("st_sd_distance")

#' @export
#' @rdname sd_distance
st_sd_distance.sfg <- function(.x, centre = NULL, weights = NULL, ...){
  .x <- st_geometry(.x)
  st_sd_distance.sfc(.x, centre = centre, weights = weights, ...)
}

#' @export
#' @rdname sd_distance
st_sd_distance.sf <- function(.x, centre = NULL, weights = NULL, ...){
  .x <- st_geometry(.x)
  st_sd_distance.sfc(.x, centre = centre, weights = weights, ...)
}

#' @export
#' @rdname sd_distance
st_sd_distance.sfc <- function (.x,
                            centre = NULL,
                            weights = NULL,
                            ...) {
  if(is.na(sf::st_crs(.x))){
    warning("st_crs(.x) returned NA, asuming EPSG:4326")
    sf::st_crs(.x) <- 4326
  }

  if (is.null(weights)) {
    weights <- rep(1, nrow(st_coordinates(.x)))}

  if(is.null(centre)) {
    centre <- .mean_centre(.x = .x, weights = weights)
  } else {
    centre_class <- class(centre)[1]
    centre <- switch(centre_class,
                     "matrix" = ,
                     "numeric" = st_sfc(st_point(centre), crs = st_crs(.x)),
                     "data.frame" = st_sfc(st_point(as.matrix(centre[1,1:2])),
                                           crs = st_crs(.x)),
                     "list" = st_sfc(st_point(as.matrix(as.data.frame(centre))),
                                     crs = st_crs(.x)),
                     "sfc_POINT" = centre[1],
                     "sf" = st_geometry(centre),
                     centre
                     )
  }

  dist <- st_distance(.x, centre, ...)[,1]
  SDD <- sqrt(sum((weights * dist ^ 2) / (sum(weights) - 2)))
  SDD.sfc <-  st_buffer(centre, dist = SDD)

  perim <- if(sf::st_is_longlat(SDD.sfc)){
    sf::st_length(SDD.sfc)
  } else {
    lwgeom::st_perimeter_lwgeom(SDD.sfc)
  }
  st_sf(feature = "Standard deviation distance",
        radius = SDD,
        area = st_area(SDD.sfc),
        perimeter = perim,
        centre = st_coordinates(centre),
        weigted = ifelse(all(weights == 1 ), FALSE, TRUE),
        geometry = st_buffer(centre, dist = SDD)
  )
}

