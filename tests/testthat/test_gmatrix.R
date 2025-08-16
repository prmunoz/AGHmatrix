test_that("VanRaden method computes correctly", {
  M <- matrix(c(0, 1, 2,
                1, 0, 2), nrow = 2, byrow = TRUE)
  colnames(M) <- c("M1", "M2", "M3")
  G <- Gmatrix(M, method = "VanRaden", ploidy = 2)
  expected <- matrix(c(0.6666667, -0.6666667,
                       -0.6666667, 0.6666667), 2)
  expect_equal(G[1:2,1:2], expected, tolerance = 1e-6, ignore_attr = TRUE)
})

test_that("Yang method computes expected G matrix values", {
  M <- matrix(c(0, 1, 2,
                1, 0, 2), nrow = 2, byrow = TRUE)
  colnames(M) <- c("M1", "M2", "M3")
  
  G <- Gmatrix(M, method = "Yang", ploidy = 2)
  
  # M1 and M2 have allele frequency p = 0.25 -> 2p(1-p) = 0.375
  # M3 is monomorphic so ignored in variance
  # Centered Z matrix for polymorphic markers:
  # M1: [-0.5, 0.5], M2: [0.5, -0.5]
  # Product terms:
  # G[1,1] = (0.25/0.375 + 0.25/0.375)/3 = 0.4444 + adjustment = 0.7778
  # G[1,2] = (-0.25/0.375 + -0.25/0.375)/3 = -0.4444
  expected <- matrix(c(0.7777778, -0.4444444,
                       -0.4444444, 0.7777778), nrow = 2, byrow = TRUE)
  
  expect_true(isSymmetric(G))
  expect_equal(round(G[1:2,1:2], 7), round(expected, 7), ignore_attr = TRUE)
})


test_that("Su method computes expected matrix", {
  M <- matrix(c(0, 1, 2,
                2, 1, 0), nrow = 2, byrow = TRUE)
  colnames(M) <- c("M1", "M2", "M3")
  
  G <- Gmatrix(M, method = "Su", ploidy = 2)
  expected <- matrix(c(1, 1,
                       1, 1), 2, 2)
  
  expect_equal(G[1:2,1:2], expected, ignore_attr = TRUE)
})

test_that("Vitezica method computes expected matrix", {
  M <- matrix(c(0, 1, 2,
                2, 1, 0), nrow = 2, byrow = TRUE)
  colnames(M) <- c("M1", "M2", "M3")
  
  G <- Gmatrix(M, method = "Vitezica", ploidy = 2)
  expected <- matrix(c(1, 1,
                       1, 1), 2, 2)
  
  expect_equal(G[1:2,1:2], expected, ignore_attr = TRUE)
})

test_that("Slater method handles zero-variance case gracefully", {
  M <- matrix(c(0, 2, 4,
                4, 2, 0), nrow = 2, byrow = TRUE)
  colnames(M) <- c("M1", "M2", "M3")
  
  expect_warning({
    G <- Gmatrix(M, method = "Slater", ploidy = 4)
    expect_false(any(is.nan(G[1:2,1:2])))
    expect_true(isSymmetric(G[1:2,1:2]))
  }, regexp = "zero variance")
})

test_that("Slater method works with non-zero variance markers", {
  # Input
  M <- matrix(c(0, 2, 4,
                4, 2, 0), nrow = 2, byrow = TRUE)
  colnames(M) <- c("M1","M2","M3")
  ploidy <- 4
  
  ## Legacy Slater steps
  
  # Flip dosages where mean > ploidy/2 
  P <- colMeans(M, na.rm = TRUE)
  M_flip <- M
  if (any(P > ploidy/2)) {
    idx <- which(P > ploidy/2)
    M_flip[, idx] <- ploidy - M_flip[, idx]
  }
  
  # presence/absence encode per Slater
  M_sl <- slater_par(M_flip, ploidy = ploidy)
  
  # denominator used by legacy code
  NumberMarkers_legacy <- ncol(M_sl)
  
  # Presence/absence frequency
  Frequency <- colMeans(M_sl, na.rm = TRUE)
  
  # drop columns with zero frequency
  drop <- which(Frequency == 0)
  if (length(drop) > 0L) {
    M_use     <- M_sl[, -drop, drop = FALSE]
    Frequency <- Frequency[-drop]
  } else {
    M_use <- M_sl
  }
  
  # FreqP and FreqPQ
  n <- nrow(M_use)
  FreqP  <- matrix(rep(Frequency, each = n), nrow = n)
  FreqPQ <- matrix(rep(Frequency * (1 - Frequency), each = n), nrow = n)
  
  # Compute G
  Z <- (M_use - FreqP) / sqrt(FreqPQ)
  Z[!is.finite(Z)] <- 0  # legacy sets NaN from 0/0 to 0
  
  G_expected <- (Z %*% t(Z)) / NumberMarkers_legacy
  
  # Diagonal
  diag_terms <- (M_use^2 - 2 * FreqP * M_use + FreqP^2) / FreqPQ
  Gii <- 1 + rowSums(diag_terms, na.rm = TRUE) / NumberMarkers_legacy
  diag(G_expected) <- Gii
  
  #      [,1]      [,2]
  # [1,] 1.266667 -0.266667
  # [2,] -0.266667 1.266667
  
  # Actual G from C++
  G <- Gmatrix(M, method = "Slater", ploidy = 4)
  
  expect_equal(round(G[1:2, 1:2], 6), round(G_expected, 6))
  expect_true(isSymmetric(G[1:2, 1:2]))
})