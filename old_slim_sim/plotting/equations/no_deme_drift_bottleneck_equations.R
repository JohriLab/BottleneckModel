#nolint start
no_deme_drift_bottleneck_coalescent <- function(m, n, N, k, u){

    # Defining substitutions
    a <- (1 - m)^2 + (m^2)/(n-1)
    b <- m * (2 - m)/(n)

    # Defining A matrix
    a11 <- 1 - a *(1 - 1/(2*k)) * (1 - 1/(2*N))
    a12 <- -(1 - a) * (1 - 1/(2*k))
    a21 <- -b * (1 - 1/(2*N))
    a22 <- 1 - (1 - b)

    A <- matrix(c(a11, a12, a21, a22), nrow = 2, byrow = TRUE)

    # Defining B vector
    B <- c(1, 1)

    # if (det(A) <= 1e-10) {
    #     stop(sprintf("Error: Matrix A is singular (det(A) <= 1e-10). Parameters: m=%g, n=%g, N=%g, k=%g, u=%g", m, n, N, k, u))
    # }

    solution_vector <- solve(A) %*% B

    T0 <- solution_vector[1]
    T1 <- solution_vector[2]

    TT <- T0 / n + (1 - 1/n) * T1

    piS <- 2 * u * T0
    piT <- 2 * u * TT
    FST <- (TT - T0)/TT

    return(list(no_deme_drift_bottleneck_T0 = T0, no_deme_drift_bottleneck_T1 = T1))

}

#nolint end