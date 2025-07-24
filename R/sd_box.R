#' @title Standard deviation box calculator in 2D or 3D
#' @description Calculate the spatial deviaction box from a points sf dataset.
#' #' @author Gabriel Gaona
#' @param .x  \code{sf} points 2D or 3D
#' @param centre  Numeric. Coordinates 2D or 3D of central point. Default NULL,
#'   performs a calculation of mean_centre() from point localities
#' @param weights Numeric. Same length of number of .x.
#' @param ... ignored
#' @return Depends on input, "coords" returns a data.frame of 2 or 3 columns and
#'   4 or 8 point coordinates. "param" returns a data.frame with centre
#'   coordinates, standard deviation in each axis, space(area for 2D, volume for
#'   3D) and number of dimensions in coordinates.
#' @importFrom Hmisc wtd.var
#' @examples
#'  requireNamespace("ggplot2", quietly = TRUE)
#'  library(sf, quietly = TRUE)
#'  library(ggplot2)
#'  bbx <- matrix(c(697047,9553483,
#'                  696158,9560476,
#'                  700964,9561425,
#'                  701745,9555358),
#'                byrow = TRUE,
#'                ncol = 2)
#'  bbx <- st_multipoint(bbx)
#'  bbx <- st_cast(bbx,"POLYGON")
#'  bbx <- st_sfc(bbx, crs = 31992)
#'  set.seed(1234)
#'  points <- st_sf(geometry = st_sample(bbx, 100))
#'  SD_BOX <- st_sd_box(points)
#'  ggplot() +
#'    geom_sf(data = SD_BOX, fill = NA, color = "darkolivegreen") +
#'    geom_sf(data = points, color = "steelblue", size = 0.5)
#' @export
#' @rdname sd_box
st_sd_box <- function(.x, centre = NULL, weights = NULL, ...) UseMethod("st_sd_box")

#' @export
#' @rdname sd_box
st_sd_box.sfg <- function(.x, centre = NULL, weights = NULL, ...){
  .x <- st_geometry(.x)
  st_sd_box.sfc(.x, centre = centre, weights = weights, ...)
}

#' @export
#' @rdname sd_box
st_sd_box.sf <- function(.x, centre = NULL, weights = NULL, ...){
  .x <- st_geometry(.x)
  st_sd_box.sfc(.x, centre = centre, weights = weights, ...)
}

#' @export
#' @rdname sd_box
st_sd_box.sfc <- function(.x, centre = NULL, weights = NULL, ...) {

  if(is.na(sf::st_crs(.x))){
    warning("st_crs(.x) returned NA, asuming EPSG:4326")
    sf::st_crs(.x) <- 4326
  }

  if (is.null(weights)) {
    weigthed <- FALSE
    weights <- rep(1, nrow(st_coordinates(.x)))
  }

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

  SD <- apply(st_coordinates(.x), 2, Hmisc::wtd.var, weights = weights) ^ (1/2)

  dirs <- t(expand.grid(rep(list(c(1, -1)), each = length(SD))))
  dirs <- dirs[,rank(apply(dirs, 2, function(x){atan2(x[2], x[1])}))]
  sd_bbox <- t(as.numeric(st_coordinates(centre)) + dirs * SD)
  sd_bbox <- st_multipoint(sd_bbox)
  sd_bbox <- st_cast(sd_bbox, "POLYGON")
  sd_bbox <- st_sf(feature = "Standard distance box",
                   geometry = st_sfc(sd_bbox, crs = st_crs(.x)))

  param <- data.frame(
    centre = st_coordinates(centre),
    sd = t(SD),
    area = st_area(sd_bbox)
  )

  return(cbind(sd_bbox, param))
}


