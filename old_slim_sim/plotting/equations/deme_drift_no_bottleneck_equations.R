#nolint start

deme_drift_coalescence <- function(m, d, N, u){

    # Defining substitutions
    alpha <- 2*m*(1-m)
    beta <- (1 - m)^2 + 2*m*(1 - m)*(d - 2)/(d - 1) + m^2 * (d - 2)/(d - 1)

    a11 = 1 - (1 - alpha)*(1 - 1/N)
    a12 = -alpha
    a21 = -1/d*(1 - alpha)*(1 - 1/N) - (1-1/d)*(1 - beta)*(1 - 1/N)
    a22 = 1 - 1/d*alpha - (1 - 1/d)*beta

    A <- matrix(c(a11, a12, a21, a22), nrow = 2, byrow = TRUE)

    b1 = 1
    b2 = 1

    B = c(b1, b2)

    solution_vector <- solve(A) %*% B

    T0 <- solution_vector[1]
    T1 <- solution_vector[2]

    return(list(deme_drift_T0 = T0, deme_drift_T1 = T1))
}


#nolint end