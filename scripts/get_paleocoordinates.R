# Generate paleocoordinates for Sargentodoxa fossil records
#
# Input:
#   data/fossil_record/Fossil_Record_Formal_v2.1_for_coding.csv
#
# Output:
#   data/fossil_record/FossilRecord-paleocor-v2.1.csv
#
# Required input columns (v2.1, trimmed):
#   No., Species, Organ, Country/Region, Latitude, Longitude, Age, Age (Ma)
#
# Output columns (only those used by the visualization scripts):
#   number, fossil_sp, organ, country_region, Age_stage, Age, paleoLong, paleoLat, age, simulation_name

library(dplyr)
library(readr)
library(httr2)
library(jsonlite)

project_root <- normalizePath(getwd(), mustWork = TRUE)
if (basename(project_root) == "scripts" && dir.exists(file.path(dirname(project_root), "data"))) {
  project_root <- dirname(project_root)
}

input_file <- file.path(project_root, "data", "fossil_record", "Fossil_Record_Formal_v2.1_for_coding.csv")
simulation_file <- file.path(project_root, "documents", "teyf_simulations.csv")
output_file <- file.path(project_root, "data", "fossil_record", "FossilRecord-paleocor-v2.1.csv")

if (!file.exists(input_file)) {
  stop("Input fossil record not found: ", input_file)
}
if (!file.exists(simulation_file)) {
  stop("TEYF simulation metadata not found: ", simulation_file)
}

cat("Reading fossil records:", input_file, "\n")
fossil_raw <- read_csv(input_file, show_col_types = FALSE)

required_cols <- c("No.", "Species", "Organ", "Country/Region",
                   "Latitude", "Longitude", "Age", "Age (Ma)")
missing_cols <- setdiff(required_cols, names(fossil_raw))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

fossil_for_reconstruction <- fossil_raw %>%
  transmute(
    number = `No.`,
    fossil_sp = Species,
    organ = Organ,
    country_region = `Country/Region`,
    decimalLongitude = as.numeric(Longitude),
    decimalLatitude = as.numeric(Latitude),
    Age_stage = Age,
    Age = as.numeric(`Age (Ma)`)
  )

if (any(!is.finite(fossil_for_reconstruction$decimalLongitude)) ||
    any(!is.finite(fossil_for_reconstruction$decimalLatitude)) ||
    any(!is.finite(fossil_for_reconstruction$Age))) {
  bad_rows <- which(
    !is.finite(fossil_for_reconstruction$decimalLongitude) |
      !is.finite(fossil_for_reconstruction$decimalLatitude) |
      !is.finite(fossil_for_reconstruction$Age)
  )
  stop("Invalid longitude, latitude, or age in rows: ", paste(bad_rows, collapse = ", "))
}

reconstruction_age <- round(fossil_for_reconstruction$Age)
coords <- fossil_for_reconstruction[, c("decimalLongitude", "decimalLatitude")]

cat("Reconstructing", nrow(coords), "records with PALEOMAP model...\n")
paleo_coords <- matrix(NA_real_, nrow = nrow(coords), ncol = 2)

reconstruct_point_gws <- function(lon, lat, age, model = "PALEOMAP") {
  req <- request("https://gws.gplates.org/reconstruct/reconstruct_points/") %>%
    req_body_form(
      lats = as.character(lat),
      lons = as.character(lon),
      time = as.character(age),
      model = model,
      reverse = "false",
      anchor_plate_id = "0",
      ignore_valid_time = "false"
    ) %>%
    # Some Windows/R setups fail the default rgplates SSL handshake even when
    # the GPlates endpoint is reachable. These options keep the script runnable.
    req_options(ssl_verifypeer = 0L, ssl_verifyhost = 0L)

  resp <- req_perform(req)
  body <- fromJSON(resp_body_string(resp))
  as.numeric(body$coordinates[1, ])
}

for (i in seq_len(nrow(coords))) {
  paleo_coords[i, ] <- reconstruct_point_gws(
    lon = coords$decimalLongitude[i],
    lat = coords$decimalLatitude[i],
    age = reconstruction_age[i],
    model = "PALEOMAP"
  )
}

paleo_coords <- as.data.frame(paleo_coords)
names(paleo_coords) <- c("paleoLong", "paleoLat")

teyf_simulations <- read_csv(simulation_file, show_col_types = FALSE)
if (!all(c("Simulation Name", "Ma") %in% names(teyf_simulations))) {
  stop("TEYF simulation metadata must include columns: Simulation Name, Ma")
}

nearest_idx <- vapply(
  reconstruction_age,
  function(x) which.min(abs(teyf_simulations$Ma - x)),
  FUN.VALUE = integer(1)
)

fossil_paleocor <- bind_cols(fossil_for_reconstruction, paleo_coords) %>%
  mutate(
    age = teyf_simulations$Ma[nearest_idx],
    simulation_name = teyf_simulations$`Simulation Name`[nearest_idx]
  ) %>%
  select(number, fossil_sp, organ, country_region, Age_stage, Age,
         paleoLong, paleoLat, age, simulation_name)

cat("Writing paleocoordinate CSV:", output_file, "\n")
write_csv(fossil_paleocor, output_file)

cat("Done. Rows:", nrow(fossil_paleocor), "\n")
cat("Mapped TEYF ages:", paste(sort(unique(fossil_paleocor$age)), collapse = ", "), "\n")
