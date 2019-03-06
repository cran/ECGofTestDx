# Some definitions
.e <- function(m, k) { 
  if (k == 0) return(1)
  if (k == 1) return(m)
  if (k >= 2) return(choose(m + k - 1, m - 1) - choose(m + k - 3, m - 1))
}
.GegenbauerC <- function(n, m, x) {
  poly <- orthopolynom::gegenbauer.polynomials(n, m, normalized = FALSE)[n + 1]
  # This function returns a list with n + 1 elements containing the order k Gegenbauer polynomials,
  # C_k^{(\alpha)}(x), for orders k = 0, 1, . . . , n
  return(orthopolynom::polynomial.values(poly, x)[[1]])
}
.ChebyshevT <- function(n, x) {
  poly <- orthopolynom::chebyshev.t.polynomials(n, normalized = FALSE)[n + 1]
  # This function returns a list with n + 1 elements containing the order k Chebyshev polynomials of
  # the first kind, T_k(x), for orders k = 0, 1, . . . , n.
  return(orthopolynom::polynomial.values(poly, x)[[1]])
}
.ChebyshevT.Cpp <- function(n, x) {
  poly <- .C("chebyshev_t_polynomials", res = as.double(x), as.integer(length(x)), as.integer(n), PACKAGE = "ECGofTestDx")
  # This function returns a list with n + 1 elements containing the order k Chebyshev polynomials of
  # the first kind, T_k(x), for orders k = 0, 1, . . . , n.
  return(poly$res)
}
.ChebyshevU <- function(n, x) {
  poly <- orthopolynom::chebyshev.u.polynomials(n, normalized = FALSE)[n + 1]
  # This function returns a list with n + 1 elements containing the order k Chebyshev polynomials of
  # the second kind, U_k(x), for orders k = 0, 1, . . . , n.
  return(orthopolynom::polynomial.values(poly, x)[[1]])
} 
.ChebyshevU.Cpp <- function(n, x) {
  poly <- .C("chebyshev_u_polynomials", res = as.double(x), as.integer(length(x)), as.integer(n), PACKAGE = "ECGofTestDx")
  # This function returns a list with n + 1 elements containing the order k Chebyshev polynomials of
  # the second kind, U_k(x), for orders k = 0, 1, . . . , n.
  return(poly$res)
}
.LaguerreL <- function(n, a, x) {
  poly <- orthopolynom::glaguerre.polynomials(n, a, normalized = FALSE)[n + 1]
  # This function returns a list with n + 1 elements containing the order k Chebyshev polynomials of
  # the second kind, U_k(x), for orders k = 0, 1, . . . , n.
  return(orthopolynom::polynomial.values(poly, x)[[1]])
}
.s <- function(k, j, m, z) { 
  res <- .GegenbauerC(k, (m + 2 * j - 2) / 2, z) *
    sqrt(gamma((m - 1) / 2) *
           (if((m / 2 + j - 1) == 0) 1 else gamma(m / 2 + j - 1)) ^ 2 *
           factorial(k) * (k + m/2 + j - 1) / (gamma(m / 2) * sqrt(pi) * 2 ^ (3 - m - 2 * j) *
                                                 gamma(k + m + 2 * j - 2)))
}
.norm <- function(X) return(sqrt(sum(X ^ 2)))
.PseudoInverse <- function(M, eps = 1e-13) {
  if (prod(dim(M)) == 1) return(1 / M)
  s <- svd(M)
  e <- s$d
  e[e>eps] <- 1/e[e>eps]
  return(s$v %*% diag(e) %*% t(s$u))
}

# recursive definition of Spherical harmonics
.psi <- function(m, k, pp, xlist, Cpp = TRUE) {
#  if (k > 100) browser()
  if (m == 1 && k == 0) return(1)
  if (m == 1 && k == 1) return(xlist[1])
  if (m == 2 && k == 0 && pp == 1) return(1)
  if (m == 2 && pp == 1) if (Cpp) return(sqrt(2) * .ChebyshevT.Cpp(k, xlist[1])) else return(sqrt(2) * .ChebyshevT(k, xlist[1]))
  if (m == 2 && pp == 2) if (Cpp) return(sqrt(2) * xlist[2] * .ChebyshevU.Cpp(k - 1, xlist[1])) else return(sqrt(2) * xlist[2] * .ChebyshevU(k - 1, xlist[1]))
  list1 <- as.list(NULL)
  for (jj in 0:k) {
    deb <- if (length(list1) == 0) 0 else rev(unlist(list1))[[1]]
    list1[[jj + 1]] <- (deb + 1):(deb + .e(m - 1, jj))
  }
  for (jj in 1:length(list1)) {
    if (pp %in% list1[[jj]]) break()
  }
  l <- pp - list1[[jj]][1] + 1
  jj <- jj - 1
  denom <- .norm(xlist[-1]) ^ 2
  res <- denom ^ (jj / 2) * .s(k - jj, jj, m, xlist[1]) * .psi(length(xlist) - 1, jj, l, xlist[-1] / sqrt(denom), Cpp)
  return(res)
}

.Phifunc <- function(k, j, l, m, r, xlist, Cpp) {
  res <- r ^ (k - 2 * j) * 
    ((-1) ^ j * sqrt((factorial(j) * gamma(m / 2)) / (2 ^ (k - 2 * j) * gamma(m / 2 + k - j))) * .LaguerreL(j, m / 2 + k - 2 * j - 1, r ^ 2 / 2)) * 
    .psi(length(xlist), k - 2 * j, l, xlist, Cpp)
  return(res)
}

SmoothECTest <- function(data, K = 7, family = "MVN", Est.Choice = "", Cpp = TRUE) {

  dataX <- as.matrix(na.omit(data))
  Kmax <- K
  
  # Initialisation 
  n <- nrow(dataX)
  m <- ncol(dataX)
  x <- paste("x", 1:m, sep = "")
  muhat <- colMeans(dataX)
  tmp <- svd((n - 1) * cov(dataX) / n)
  sqrtSigmahat <- tmp$u %*% diag(1 / sqrt(tmp$d)) %*% t(tmp$v)
  QUIRstructure <- list()
  ind <- 1
  for (k in 3:Kmax) {
      QUIRstructure[[ind]] <- rep(0, floor(k / 2) + 1)
      ind <- ind + 1
  }
  QUIRstructurescaled <- QUIRstructure

  # Tidy the data
  dataYhat <- t(sqrtSigmahat %*% (t(dataX) - muhat))
  dataRhat <- rep(NA, n)
  for (i in 1:n) dataRhat[i] <- .norm(dataYhat[i, ])
  dataUhat <- matrix(NA, nrow = n, ncol = m)
  for (i in 1:n) dataUhat[i, ] <- dataYhat[i, ] / dataRhat[i]
  
  # Computation of test statistics
  for (k in 3:Kmax) {
    for (j in 1:(floor(k / 2) + 1)) {
      ind <- .e(m, k - 2 * (j - 1))
      datatrans <- matrix(NA, nrow = n, ncol = ind)
      for (i in 1:n) {
        for (l in 1:ind) {
          datatrans[i, l] <- .Phifunc(k, j - 1, l, m, dataRhat[i], dataUhat[i, ], Cpp)
        }
      }
      QUIRstructure[[k - 2]][j] <- n * .norm(colMeans(datatrans)) ^ 2
      QUIRstructurescaled[[k - 2]][j] <- if(is.null(datatrans[[1]])) 0 else n * t(as.matrix(colMeans(datatrans))) %*% 
        .PseudoInverse(cov(datatrans)) %*% as.matrix(colMeans(datatrans))
    }
  }
  Q <- 0
  for (k in 3:Kmax) {
    for (j in 1:(floor(k / 2) + 1)) {
      Q <- Q + QUIRstructure[[k - 2]][j]
    }
  }
  R <- 0
  if (Kmax > 3) {
   for (kp in 2:(Kmax / 2)) {
     k <- 2 * kp
     j <- floor(k / 2)
     R <- R + QUIRstructure[[k - 2]][j + 1]
   }
  }
  U <- 0
  for (k in 3:Kmax) {
    j <- 0
    U <- U + QUIRstructure[[k - 2]][j + 1]
  }
  I <- 0
  for (k in 3:Kmax) {
    for (j in 1:(floor((k - 1) / 2))) {
      I <- I + QUIRstructure[[k - 2]][j + 1]
    }
  }
  Qscaled <- 0
  for (k in 3:Kmax) {
    for (j in 1:(floor(k / 2) + 1)) {
      Qscaled <- Qscaled + QUIRstructurescaled[[k - 2]][j]
    }
  }
  Rscaled <- 0
  if (Kmax > 3) {
   for (kp in 2:(Kmax / 2)) {
     k <- 2 * kp
     j <- floor(k / 2)
     Rscaled <- Rscaled + QUIRstructurescaled[[k - 2]][j + 1]
   }
  }
  Uscaled <- 0
  for (k in 3:Kmax) {
    j <- 0
    Uscaled <- Uscaled + QUIRstructurescaled[[k - 2]][j + 1]
  }
  Iscaled <- 0
  for (k in 3:Kmax) {
    for (j in 1:(floor((k - 1) / 2))) {
      Iscaled <- Iscaled + QUIRstructurescaled[[k - 2]][j + 1]
    }
  }
  
  dfQ <- 0
  for (k in 3:Kmax) dfQ <- dfQ + choose(m + k - 1, k)
  dfR <- 0
  if (Kmax > 3) {
   for (kp in 2:(Kmax / 2)) {
     k <- 2 * kp
     dfR < dfR + .e(m ,0)
   }
  }
  dfU <- 0
  for (k in 3:Kmax) dfU <- dfU + .e(m, k)
  dfI <- dfQ -dfR -dfU
  
  # Output is under the form : 
  # (Q, Qscaled, dfQ) : dimensions Kmax - 3 + 1 
  # (U, Uscaled, dfU) : dimensions Kmax - 3 + 1 
  # (I, Iscaled, dfI) : dimensions Kmax - 3 + 1 
  # (R, Rscaled, dfR) : dimensions floor((Kmax - 4) / 2) + 1 

  resdfI <- resIscaled <- resI <- resdfU <- resUscaled <- resU <- resdfQ <- resQscaled <- resQ <- rep(NA, Kmax - 3 + 1)
  id <- 1
  sommedfI <- sommeIscaled <- sommeI <- sommedfU <- sommeUscaled <- sommeU <- sommedfQ <- sommeQscaled <- sommeQ <- 0
  for (k in 3:Kmax) {
    for (jj in 1:(floor(k / 2) + 1)) {
      sommeQ <- sommeQ + QUIRstructure[[k - 2]][jj]
      sommeQscaled <- sommeQscaled + QUIRstructurescaled[[k - 2]][jj]
    }
    resQ[id] <- sommeQ
    resQscaled[id] <- sommeQscaled
    sommedfQ <- sommedfQ + choose(m + k - 1, k)
    resdfQ[id] <- sommedfQ
    
    sommeU <- sommeU + QUIRstructure[[k - 2]][1]
    resU[id] <- sommeU
    sommeUscaled <- sommeUscaled + QUIRstructurescaled[[k - 2]][1]
    resUscaled[id] <- sommeUscaled
    sommedfU <- sommedfU + .e(m, k)
    resdfU[id] <- sommedfU
    
    for (jj in 2:(floor((k -1) / 2) + 1)) {
      sommeI <- sommeI + QUIRstructure[[k - 2]][jj]
      sommeIscaled <- sommeIscaled + QUIRstructurescaled[[k - 2]][jj]
    }
    resI[id] <- sommeI
    resIscaled[id] <- sommeIscaled
    sommedfI <- sommedfI + choose(m + k - 1, k)- (if(k %% 2 == 0) .e(m, 0) else 0) - .e(m, k)
    resdfI[id] <- sommedfI

    id <- id + 1
  }
  
  resdfR <- resRscaled <- resR <- rep(NA, (Kmax/ 2) - 2 + 1)
  id <- 1
  sommedfR <- sommeRscaled <- sommeR <- 0
  if (Kmax > 3) {
    for (kp in 2:(Kmax / 2)) {
    k <- 2 * kp
    jj <- floor(k / 2)
    sommeR <- sommeR + QUIRstructure[[k - 2]][jj + 1]
    resR[id] <- sommeR
    sommeRscaled <- sommeRscaled + QUIRstructurescaled[[k - 2]][jj + 1]
    resRscaled[id] <- sommeRscaled
    sommedfR <- sommedfR + .e(m, 0)
    resdfR[id] <- sommedfR
    
    id <- id + 1
   }
  }
  
  Q <- resQ[length(resQ)]
  dfQ <- resdfQ[length(resdfQ)]
  Uscaled <- resUscaled[length(resUscaled)]
  dfU <- resdfU[length(resdfU)]
  Iscaled <- resIscaled[length(resIscaled)]
  dfI <- resdfI[length(resdfI)]
  Rscaled <- resRscaled[length(resRscaled)]
  dfR <- resdfR[length(resdfR)]

  if (length(Rscaled) != 0) {  
   res <- matrix(c(Q = round(Q, 3), dfQ = dfQ, pval.asymp.Q = round(1- pchisq(Q, dfQ), 4), "Global",
              Uscaled = round(Uscaled, 3), dfU = dfU, pval.asymp.U = round(1- pchisq(Uscaled, dfU), 4), "U dist.",  
              Iscaled = round(Iscaled, 3), dfI = dfI, pval.asymp.I = round(1- pchisq(Iscaled, dfI), 4), "cor(R,U)", 
              Rscaled = round(Rscaled, 3), dfR = dfR, pval.asymp.R = round(1- pchisq(Rscaled, dfR), 4), "R dist."), nrow = 4, ncol = 4, byrow = TRUE,
              dimnames = list(c("Q", "U-scaled", "I-scaled", "R-scaled"), c("Stat", "df", "p.val.Chi2", "Interpretation")))
  } else {
   res <- matrix(c(Q = round(Q, 3), dfQ = dfQ, pval.asymp.Q = round(1- pchisq(Q, dfQ), 4), "Global",
              Uscaled = round(Uscaled, 3), dfU = dfU, pval.asymp.U = round(1- pchisq(Uscaled, dfU), 4), "U dist.",  
              Iscaled = round(Iscaled, 3), dfI = dfI, pval.asymp.I = round(1- pchisq(Iscaled, dfI), 4), "cor(R,U)", 
              Rscaled = 0, dfR = 0, pval.asymp.R = NA, "R dist."), nrow = 4, ncol = 4, byrow = TRUE,
              dimnames = list(c("Q", "U-scaled", "I-scaled", "R-scaled"), c("Stat", "df", "p.val.Chi2", "Interpretation")))
    }

  cat(paste("\n\n Smooth Test of H0 : X ~ MVN   (K = ", Kmax, ", n * p = ", n, " * ", m, ")\n\n", sep = ""))
  
return(as.data.frame(res))
  
#  return(list(Q = Q, dfQ = dfQ, pval.asymp.Q = 1- pchisq(Q, dfQ),
#              Uscaled = Uscaled, dfU = dfU, pval.asymp.U = 1- pchisq(Uscaled, dfU), 
#              Iscaled = Iscaled, dfI = dfI, pval.asymp.I = 1- pchisq(Iscaled, dfI),
#              Rscaled = Rscaled, dfR = dfR, pval.asymp.R = 1- pchisq(Rscaled, dfR)))

}
