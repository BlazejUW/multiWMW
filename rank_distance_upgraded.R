rank_distance_uprgaded<- function(X, Y, Z) {
  n <- nrow(X)
  m <- nrow(Y)
  c <- nrow(Z)
  ZN <- Z
  DZ <- as.matrix(dist(ZN, method = "euclidean", diag = TRUE, upper = TRUE))
  DZX <- as.matrix(dist(rbind(X, ZN), method = "euclidean", diag = TRUE, upper = TRUE))
  DX <- DZX[1:n, (n+1):(c+n)]
  DZY <- as.matrix(dist(rbind(Y, ZN), method = "euclidean", diag = TRUE, upper = TRUE))
  DY <- DZY[1:m, (m+1):(c+m)]
  
  WX <- matrix(0, n, c)
  WY <- matrix(0, m, c)
  Tz=NULL

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
  