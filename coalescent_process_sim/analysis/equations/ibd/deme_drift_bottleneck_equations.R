#nolint start


deme_drift_bottleneck_ibd_Ne <- function(m, d, N, k, u){

    # Defining substitutions
    alpha <- 2*m*(1-m)
    beta <- (1 - m)^2 + 2*m*(1 - m)*(d - 2)/(d - 1) + m^2 * (d - 2)/(d - 1)

    a11 = 1 - (1 - 1/k)*(1 - alpha)*(1 - 1/N)*(1 - u)^2
    a12 = -(1 - 1/k)*alpha*(1 - u)^2
    a21 = -1/d*(1 - alpha)*(1 - 1/N)*(1 - u)^2 - (1-1/d)*(1 - beta)*(1 - 1/N)*(1 - u)^2
    a22 = 1 - 1/d*alpha*(1 - u)^2 - (1 - 1/d)*beta*(1 - u)^2

    A <- matrix(c(a11, a12, a21, a22), nrow = 2, byrow = TRUE)

    b1 = (1 - u)^2 * (1/k) + (1 - u)^2 * (1 - 1/k) * (1 - alpha) * (1/N)
    b2 = (1 - u)^2 * (1/d) * (1 - alpha) * (1/N) + (1 - u)^2 * (1 - 1/d) * (1 - beta) * (1/N)

    B = c(b1, b2)

    solution_vector <- solve(A) %*% B

    f0 <- solution_vector[1]
    f1 <- solution_vector[2]

    f_bar <- 1/d * f0 + (1 - 1/d) * f1

    Ne <- 1/f_bar

    return(Ne)
}

deme_drift_bottleneck_ibd <- function(m, d, N, k, u){

    # Defining substitutions
    alpha <- 2*m*(1-m)
    beta <- (1 - m)^2 + 2*m*(1 - m)*(d - 2)/(d - 1) + m^2 * (d - 2)/(d - 1)

    a11 = 1 - (1 - 1/k)*(1 - alpha)*(1 - 1/N)*(1 - u)^2
    a12 = -(1 - 1/k)*alpha*(1 - u)^2
    a21 = -1/d*(1 - alpha)*(1 - 1/N)*(1 - u)^2 - (1-1/d)*(1 - beta)*(1 - 1/N)*(1 - u)^2
    a22 = 1 - 1/d*alpha*(1 - u)^2 - (1 - 1/d)*beta*(1 - u)^2

    A <- matrix(c(a11, a12, a21, a22), nrow = 2, byrow = TRUE)

    b1 = (1 - u)^2 * (1/k) + (1 - u)^2 * (1 - 1/k) * (1 - alpha) * (1/N)
    b2 = (1 - u)^2 * (1/d) * (1 - alpha) * (1/N) + (1 - u)^2 * (1 - 1/d) * (1 - beta) * (1/N)

    B = c(b1, b2)

    solution_vector <- solve(A) %*% B

    f0 <- solution_vector[1]
    f1 <- solution_vector[2]

    f_bar <- 1/d * f0 + (1 - 1/d) * f1

    Ne <- 1/f_bar

    return(f_bar)
}




#nolint end