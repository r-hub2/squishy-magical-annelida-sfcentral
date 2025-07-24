test_that("Spatial Standard Distance works", {
  x <- y <- sf::st_sf(geometry = sf::st_sfc(list(sf::st_point(c(1, 1)),
                                            sf::st_point(c(2, 2)),
                                            sf::st_point(c(3, 3)))))
  sf::st_crs(x) <- 4326
  sdd <- st_sd_distance(x)
  expect_match("Standard deviation distance", sdd$feature)
  expect_s3_class(sdd, class(x))
  expect_no_message(st_sd_distance(x))
  expect_no_warning(st_sd_distance(x))
  expect_warning(st_sd_distance(y))
})
