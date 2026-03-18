AUC <- function(score, truth) {
  n1 <- sum(truth == 1)
  n0 <- sum(truth == 0)
  if (n1 == 0 || n0 == 0) return(NA)
  r <- rank(score, ties.method = "average")
  (sum(r[truth == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
