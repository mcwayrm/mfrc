argument <- grep("^--file=", commandArgs(), value = TRUE)
if (length(argument) != 1L) stop("Run with Rscript.")
root <- dirname(dirname(normalizePath(sub("^--file=", "", argument))))
source(file.path(root, "scripts", "task-value-EVSD-implied.R"))

record <- function(id, service, label, biome, value, change = "Total", include = NA_character_) {
  data.frame(
    ValueId = id, StudyId = paste0("study-", id), ESVD2.0_Biome_Codes = biome,
    ESVD2.0_Biome = biome, TEEB_ES = service, ES_1 = label,
    Value_Change = change, Incl_Excl = include,
    `Int$ Per Hectare Per Year` = value, Reference = "test",
    check.names = FALSE, stringsAsFactors = FALSE
  )
}
forest <- do.call(rbind, list(
  record("a", "8", "Climate", "6", "10"),
  record("b", "8", "Climate", "6; 7", "30"),
  record("c", "8", "Climate", "6; 4", "999"),
  record("d", "8; 9", "Climate", "6", "999"),
  record("e", "8", "Climate", "6", "999", change = "Positive"),
  record("f", "8", "Climate", "6", "999", include = "0"),
  record("g", "15", "Biological control", "6", "7")
))
all_biomes <- do.call(rbind, list(
  record("h", "4", "Genetic resources", "2", "4"),
  record("i", "15", "Biological control", "2", "999"),
  record("j", "16", "Life cycles", "3", "16"),
  record("k", "20", "Inspiration", "1", "20"),
  record("l", "21", "Spiritual", "8", "21"),
  record("m", "16", "Life cycles", "3; 4", "999")
))
selected <- select_esvd_values(forest, all_biomes)
stopifnot(setequal(selected$value_id, c("a", "b", "g", "h", "j", "k", "l")),
          selected$source[selected$value_id == "g"] == "forest biomes",
          selected$source[selected$value_id == "h"] == "all biomes (fallback)")
prices <- summarize_esvd(selected)
climate <- prices[prices$service_code == "8", ]
stopifnot(climate$observations == 2, climate$studies == 2,
          climate$mean_int2020_per_ha_per_year == 20,
          climate$minimum == 10, climate$maximum == 30)

invalid <- forest
invalid[["Int$ Per Hectare Per Year"]][1] <- "not-a-number"
stopifnot(inherits(try(select_esvd_values(invalid, all_biomes), silent = TRUE), "try-error"))
missing_fallback <- all_biomes[all_biomes$TEEB_ES != "21", ]
stopifnot(inherits(try(select_esvd_values(forest, missing_fallback), silent = TRUE), "try-error"))
cat("ESVD filtering, fallback precedence, arithmetic means and invalid-input checks passed.\n")
