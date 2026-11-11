library("tidyverse")
library("here")
library("readxl")
library("ape")

site_ref <- read_excel(here("data_raw", "Plot-IDs_Abgleich_alleUntersuchungen.xlsx")) |>
  select(ID, Site_ID = Sample_ID)

coords_df <- read_excel(here("data_raw", "mastertable.xlsx")) |>
  select(ID, x_coord = Lon, y_coord = Lat) |>
  left_join(site_ref, by = "ID") |>
  select(Site_ID, x_coord, y_coord) |>
  filter(!is.na(x_coord), !is.na(y_coord))

stopifnot("x_coord not in UTM 32N easting range" = all(between(coords_df$x_coord, 3e5, 7e5)),
          "y_coord not in UTM 32N northing range" = all(between(coords_df$y_coord, 5.5e6, 5.8e6)))

moran_resid <- function(mod, site_ids) {
  xy <- coords_df[match(site_ids, coords_df$Site_ID), c("x_coord", "y_coord")]
  w <- 1 / as.matrix(dist(xy))
  diag(w) <- 0
  Moran.I(residuals(mod, type = "deviance"), w)
}