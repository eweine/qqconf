library(dplyr)
alpha_01_df <- readr::read_csv("~/downloads/alpha_01_500k.csv")
alpha_05_df <- readr::read_csv("~/Documents/local_level_results.csv")
df_100k <- data.frame(n=c(100000), local_level=c(.0004781731))
alpha_05_df <- rbind(alpha_05_df, df_100k)

# edits to .01 df based on more accurate calculations
exclude_n <- c(15000, seq(from = 50000, to = 500000, by = 50000))
alpha_01_df <- alpha_01_df %>%
  filter(!(n %in% exclude_n))

new_loc_lev <- c(
  9.822814e-05, 8.200296e-05, 7.476107e-05, 7.105596e-05, 6.862968e-05,
  6.685227e-05, 6.546293e-05, 6.432979e-05, 6.337766e-05, 6.255951e-05,
  6.184428e-05
)

new_alpha_01_df <- data.frame(n = exclude_n, local_level = new_loc_lev)
alpha_01_df <- rbind(alpha_01_df, new_alpha_01_df) %>% arrange(n)

usethis::use_data(alpha_05_df, alpha_01_df, internal = TRUE, overwrite = TRUE)
