rank_distance_upgraded_cpp <- function(X, Y, Z, use_avx512 = TRUE, use_parallel = TRUE, use_float = FALSE) {
  n <- nrow(X)
  m <- nrow(Y)
  c <- nrow(Z)
  ZN <- Z

  # Sprawdzenie wsparcia dla AVX-512
  if (use_avx512 && !check_avx512_support()) {
    message("AVX-512 not supported on this machine. Falling back to non-AVX-512 method.")
    use_avx512 <- FALSE
  }

  # Wybór odpowiedniej funkcji na podstawie argumentów use_avx512, use_parallel i use_float
  if (use_float) {
    if (use_avx512 && use_parallel) {
      dist_function <- dist_avx512_parallel_float
      method_name <- "dist_avx512_parallel_float"
    } else if (use_avx512 && !use_parallel) {
      dist_function <- dist_avx512_float
      method_name <- "dist_avx512_float (without parallel)"
    } else if (!use_avx512 && use_parallel) {
      dist_function <- dist_no_avx_parallel_float
      method_name <- "dist_no_avx_parallel_float"
    } else {
      dist_function <- dist_no_avx_no_parallel_float
      method_name <- "dist_no_avx_no_parallel_float"
    }
  } else {
    if (use_avx512 && use_parallel) {
      dist_function <- dist_avx512_parallel
      method_name <- "dist_avx512_parallel"
    } else if (use_avx512 && !use_parallel) {
      dist_function <- dist_avx512
      method_name <- "dist_avx512 (without parallel)"
    } else if (!use_avx512 && use_parallel) {
      dist_function <- dist_no_avx_parallel
      method_name <- "dist_no_avx_parallel"
    } else {
      dist_function <- function(mat) as.matrix(dist(mat, method = "euclidean", diag = TRUE, upper = TRUE))
      method_name <- "dist() from standard R library"
    }
  }

  # Wyświetlenie informacji o wybranej metodzie
  message("Using method: ", method_name)

  DZ <- dist_function(ZN)
  DZX <- dist_function(rbind(X, ZN))
  DZY <- dist_function(rbind(Y, ZN))
  
  DX <- DZX[1:n, (n+1):(c+n)]
  DY <- DZY[1:m, (m+1):(c+m)]
  
  WX <- matrix(0, n, c)
  WY <- matrix(0, m, c)
  Tz <- numeric(c)

  for (i in 1:c) {
    sorted_DZ <- sort(DZ[i, ])
    
    for (k in 1:n) {
      WX[k, i] <- sum(sorted_DZ <= DX[k, i])
    }
    
    for (l in 1:m) {
      WY[l, i] <- sum(sorted_DZ <= DY[l, i])
    }
    
    WX_sum_i <- sum(WX[, i])
    WY_sum_i <- sum(WY[, i])
    Tz[i] <- (n + m)^(-2) * (WX_sum_i / n - WY_sum_i / m)^2
  }
  
  stat <- mean(Tz)
  return(stat)
}
