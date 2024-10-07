stratify <- function(...) {
  return(as.numeric(interaction(...)))
}

usual.test <- function(Y, A, S, X, pi, q) {
  A = cat.to.int(A)
  S = cat.to.int(S) + 1
  X = cat.to.int(X)

  if (max(X) == 1) {
    return(old.test(Y, A, S, X, pi, q))
  } else {
      return(old.test.multi(Y, A, S, X, pi, q))
    }
}

modified.test <- function(Y, A, S, X, pi, q) {
  A = cat.to.int(A)
  S = cat.to.int(S) + 1
  X = cat.to.int(X)

  if (max(X) == 1) {
    return(new.test(Y, A, S, X, pi, q))
  } else {
    return(new.test.multi(Y, A, S, X, pi, q))
  }
}

stratified.adjusted.test <- function(Y, A, S, X, pi, q) {
  A = cat.to.int(A)
  S = cat.to.int(S) + 1
  X = cat.to.int(X)

  if (max(X) == 1) {
    return(strata.test(Y, A, S, X, pi, q))
  } else {
    return(strata.test.multi(Y, A, S, X, pi, q))
  }
}


