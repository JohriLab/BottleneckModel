#nolint start
WF_coalescent_k <- function(m, d, N, u){

    # Substitutions
    a <- (1 - m)^2 + (m^2)/(d-1)
    b <- m * (2 - m)/(d)


    # Define A
    a11 <- 1 - a * (1 - 1/(N))
    a12 <- -(1 - a)
    a21 <- -b * (1 - 1/(N))
    a22 <- 1 - (1 - b)

    A <- matrix(c(a11, a12, a21, a22), nrow = 2, byrow = TRUE)

    # Define B
    b1 <- 1
    b2 <- 1

    B <- c(b1, b2)

    solution_vector <- solve(A) %*% B

    T0 <- solution_vector[1]
    T1 <- solution_vector[2]

    return(list(T0_WF = T0, T1_WF = T1))

}



#nolint end