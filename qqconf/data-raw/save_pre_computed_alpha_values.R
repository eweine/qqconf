alpha_01_df <- readr::read_csv("~/downloads/alpha_01_500k.csv")
alpha_05_df <- readr::read_csv("~/Documents/local_level_results.csv")
usethis::use_data(alpha_05_df, alpha_01_df, internal = TRUE, overwrite = TRUE)
