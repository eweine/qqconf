library(readr)
library(dplyr)

approx_table_df <- readr::read_csv(
  "~/downloads/alpha_01_500k.csv"
)

# now, want to build a function that can do the approximation with the correction
# this should be vectorizable

# have to figure out what to feed in here as input
get_approx <- function(alpha, n, c_alpha) {

  eta_approx = -log(1-alpha)/(2*log(log(n))*log(n))*(1-c_alpha*log(log(log(n)))/log(log(n)))

}

approx_alphas <- get_approx(
  alpha = .01,
  n = approx_table_df$n,
  c_alpha = 1.61
)

approx_table_df <- approx_table_df %>%
  mutate(local_level_approx = approx_alphas)

approx_table_df <- approx_table_df %>%
  mutate(error = local_level_approx - local_level)

approx_table_df <- approx_table_df %>%
  mutate(relative_error = error / local_level)

approx_table_df <- approx_table_df %>%
  filter(n >= 1000)

plot(approx_table_df$n, approx_table_df$relative_error, xlab = "n", ylab = "Relative Error", main = "C = 1.61")
abline(h = 0)

readr::write_csv(approx_table_df, "~/Documents/local_levels_relative_error.csv")
