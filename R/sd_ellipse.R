#' @title Standard deviation ellipse calculator
#' @description Calculate the spatial deviaction ellipse from a points sf dataset.
#' @author Gabriel Gaona
#' @param .x  \code{sf} points 2D or 3D
#' @param centre  Numeric. Coordinates 2D of central point. Default NULL,
#'   performs a calculation of `mean_centre()` from point localities
#' @param weights Numeric. Same length of number of points.
#' @param ... ignored
#' @return simple features as "POLYGON"  with atributes:
#' centre coordinates, values for mayor and minor axis radius (sigma.x and sigma.y),
#' rotation (theta and theta_corrected) and geometry properties (eccentricity, area
#' and perimeter)
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
#'   SDE <- st_sd_ellipse(points)
#'   ggplot() +
#'     geom_sf(data = SDE, fill = NA, color = "darkolivegreen") +
#'     geom_sf(data = points, color = "steelblue", size = 0.5)
#' @export
#' @rdname sd_ellipse
st_sd_ellipse = function(.x, centre = NULL, weights = NULL, ...) UseMethod("st_sd_ellipse")

#' @export
#' @rdname sd_ellipse
st_sd_ellipse.sfg <- function(.x, centre = NULL, weights = NULL, ...){
  .x <- st_geometry(.x)
  st_sd_ellipse.sfc(.x, centre = centre, weights = weights, ...)
}

#' @export
#' @rdname sd_ellipse
st_sd_ellipse.sf <- function(.x, centre = NULL, weights = NULL, ...){
  .x <- st_geometry(.x)
  st_sd_ellipse.sfc(.x, centre = centre, weights = weights, ...)
}

#' @export
#' @rdname sd_ellipse
st_sd_ellipse.sfc <- function (.x,
                       centre = NULL,
                       weights = NULL,
                       ...) {
  if(is.na(sf::st_crs(.x))){
    warning("st_crs(.x) returned NA, asuming EPSG:4326")
    sf::st_crs(.x) <- 4326
  }

  # Control for weights
  if (is.null(weights)) {
    weigthed <- FALSE
    weights <- rep(1, nrow(st_coordinates(.x)))
  }
  # Calculate centre if is missing
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

  #prepare coordinates for calcules

  crds  <- st_coordinates(.x)
  n <- nrow(crds)
  cen <- as.numeric(st_coordinates(centre))
  crds2 <- crds ^ 2
  crdsc <- t(t(crds) - cen)
  crdsc2 <- crdsc ^ 2
  crdscp <- as.numeric(apply(crdsc, 1, prod))

  top1 <- sum(rowSums(t(t(weights * crdsc2) * c(1,-1))))
  top2 <- sqrt(sum(top1) ^ 2 + 4 * sum(weights * crdscp) ^ 2)
  bottom <- (2 * sum(weights * crdscp))
  tantheta <- (top1 + top2) / bottom

  if (tantheta < 0) {
    theta <- 180 + (tantheta * 180 / pi)
  } else {
    theta <- (tantheta * 180 / pi)
  }
  sintheta <- sin(theta * pi / 180)
  costheta <- cos(theta * pi / 180)
  sin2theta <- sintheta ^ 2
  cos2theta <- costheta ^ 2
  sinthetacostheta <- sintheta * costheta
  if (weigthed) {
    sigmax <-
      sqrt(2) * sqrt((sum(colSums(weights * crdsc2) * c(cos2theta, sin2theta)) -
        2 * sum(weights * crdscp) * sinthetacostheta) / sum(weights))

    sigmay <-
      sqrt(2) * sqrt((sum(colSums(weights * crdsc2) * c(sin2theta, cos2theta)) +
        2 * (sum(weights * crdscp)) *  sinthetacostheta) / sum(weights))
  } else {
    sigmax <-
      sqrt(2) * sqrt((sum(colSums(crdsc2) * c(cos2theta, sin2theta)) -
        2 * sum(crdscp) * sinthetacostheta) / (n - 2))
    sigmay <-
      sqrt(2) * sqrt((sum(colSums(crdsc2) * c(sin2theta, cos2theta)) +
        2 * (sum(crdscp)) *  sinthetacostheta) /
        (n - 2))
  }

  sigmas <- sort(c(sigmax, sigmay))
  lengthsigmax <- 2 * sigmax
  lengthsigmay <- 2 * sigmay
  areaSDE <- pi * sigmax * sigmay
  eccentricity <- sqrt(1 - ((min(sigmax, sigmay) ^ 2) / (max(sigmax, sigmay) ^ 2)))
  B <- min(sigmax, sigmay)
  A <- max(sigmax, sigmay)
  d2 <- (A - B) * (A + B)
  phi <- 2 * pi * seq(0, 1, len = 360)
  sp <- sin(phi)
  cp <- cos(phi)
  r <- sigmax * sigmay / sqrt(B ^ 2 + d2 * sp ^ 2)
  xy <- r * cbind(cp, sp)
  al <- (90 - theta) * pi / 180
  ca <- cos(al)
  sa <- sin(al)
  coordsSDE <- xy %*% rbind(c(ca, sa), c(-sa, ca)) +
    matrix(cen, nrow = 360, ncol = 2, byrow = TRUE)



  if (sigmax < sigmay) {
    Theta.Corr <- theta
  }  else {
    Theta.Corr <- theta + 90
  }

  SDE.sfc <- st_sfc(
    st_cast(st_multipoint(coordsSDE), "POLYGON"),
    crs = st_crs(.x)
    )

  perim <- if(sf::st_is_longlat(SDE.sfc)){
    sf::st_length(SDE.sfc)
  } else {
    lwgeom::st_perimeter_lwgeom(SDE.sfc)
  }

  st_sf(feature = "Standard deviation ellipse",
        centre = st_coordinates(centre),
        sigma.mayor = sigmas[2],
        sigma.minor = sigmas[1],
        theta = theta,
        theta_corrected = Theta.Corr,
        eccentricity = eccentricity,
        area = st_area(SDE.sfc),
        perimeter = perim,
        weighted = ifelse(all(weights == 1), FALSE, TRUE),
        geometry = st_sfc(st_cast(st_multipoint(coordsSDE), "POLYGON"),
                          crs = st_crs(.x)))

}


