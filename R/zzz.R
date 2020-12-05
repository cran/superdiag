# HELPER: burnin
burn <- function(input.matrix, burnin) {
  out <- tail(x = input.matrix, n = (nrow(input.matrix)-burnin))
  return(out)
}


