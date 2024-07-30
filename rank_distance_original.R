rank_distance_original=function(X,Y,Z){
  n=nrow(X)
  m=nrow(Y)
  ZN=Z
  DZ=as.matrix(dist(ZN, method = "euclidean", diag = TRUE, upper = TRUE))
  DZX=as.matrix(dist(rbind(X,ZN), method = "euclidean", diag = TRUE, upper = TRUE))
  DX=DZX[(1:n),((n+1):(m+2*n))]
  DZY=as.matrix(dist(rbind(Y,ZN), method = "euclidean", diag = TRUE, upper = TRUE))
  DY=DZY[(1:m),((m+1):(2*m+n))]
  Tz=NULL
  for (i in (1:(m+n))){
    WX=0
    for (k in (1:n)){
      WX=WX+sum(DZ[i,]<=DX[k,i])
    }
    WY=0
    for (l in (1:m)){
      WY=WY+sum(DZ[i,]<=DY[l,i])
    }
    Tz[i]=(m+n)^(-2)*(WX/n-WY/m)^2  
  }
  stat=mean(Tz)
  return(stat)
}
