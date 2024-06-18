read_clariostar <- function(file, readings = 1L) {
  
  res <- read_delim(
    file = file,
    skip = 8L,
    n_max = 16L,
    delim = ",",
    col_names = c("Row", 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12),
    col_types = "cdddddddddddd",
    show_col_types = FALSE
  ) |>
  pivot_longer(
    data = _,
    cols = !Row,
    names_to = "Column",
    values_to = "A562"
  ) |>
  drop_na(data = _) |>
  mutate(.data = _, Well = paste0(Row, Column), .before = 1L) |>
  select(.data = _, -Row, -Column)
  
  return(res)
  
}


