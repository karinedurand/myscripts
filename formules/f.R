f <- function(x){
  m <- mean(x, na.rm = TRUE)
  x[is.na(x)] <- m
  x
}