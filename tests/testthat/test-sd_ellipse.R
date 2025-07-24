test_that("Spatial Standard Ellipse works", {
  x <- y <- sf::st_sf(geometry = sf::st_sfc(list(sf::st_point(c(1, 1)),
                                                 sf::st_point(c(2, 2)),
                                                 sf::st_point(c(3, 3)))))
  sf::st_crs(x) <- 4326
  sde <- st_sd_ellipse(x)
  expect_match("Standard deviation ellipse", sde$feature)
  expect_s3_class(sde, class(x))
  expect_no_message(st_sd_ellipse(x))
  expect_no_warning(st_sd_ellipse(x))
  expect_warning(st_sd_ellipse(y))
})
