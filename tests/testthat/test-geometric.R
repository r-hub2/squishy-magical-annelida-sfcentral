test_that("geometric works", {
  x <- sf::st_sf(geometry = sf::st_sfc(list(sf::st_point(c(1, 1)),
                                            sf::st_point(c(2, 2)),
                                            sf::st_point(c(3, 3)))))
  y <- sf::st_coordinates(x)
  expect_equal(sf::st_coordinates(st_central_point(x, method = "geometric")),
               matrix(apply(y, 2,
                      FUN = function(x){
                        from <- range(x)
                        to <- 0:1
                        crd <- prod((scales::rescale(x, to))^(1/length(x)))
                        scales::rescale(crd, from)
                      }),
                      ncol = 2,
                      dimnames = list(NULL, c("X", "Y"))))
  expect_equal(nrow(st_central_point(x, method = "geometric")), 1)
  expect_error(st_central_point(y, method = "median"))
  expect_s3_class(st_central_point(x, method = "median"), c("sf", "data.frame"))
})
