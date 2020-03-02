get_polynom <- function(deg, bc_points, other_points){
  
  mat <- matrix(nrow = deg+1, ncol = deg+1)
  y <- c()
  x <- c(bc_points[1], other_points, bc_points[3])
  for(points in 2:(deg)){
    y[points] <- sin(3*x[points])
  }
  y[1] <- bc_points[2]
  y[length(x)] <- bc_points[4]
  
  for(degs in 0:deg){
    sec_pol <- degs-2
    if(sec_pol < 0) {
      mat[,degs+1] <- 4*x^degs
    }else{
      mat[,degs+1] <- -degs*(degs-1)*x^sec_pol+4*x^degs
    }
  }
  # boundaries
  for(bc1 in 1:(deg+1)){
    mat[c(1,length(x)), bc1] <- c(bc_points[1],bc_points[3])^(bc1-1)
  }
  print("assembled matrix")
  print(mat)
  print("assembled b-vector")
  print(y)
  for(k in 1:(ncol(mat)-1)){
    for(i in (k+1):ncol(mat)){
      xmult <- mat[i,k]/mat[k,k]
      for(j in k:ncol(mat)){
        mat[i,j]<-mat[i,j]-xmult*mat[k,j] 
        if(abs(mat[i,j]) < 1e-14){mat[i,j] <- 0}
      }
      y[i] <- y[i]-xmult*y[k]
      
    }
  }
  print("matrix after Gaussian forward elimination")
  print(mat)
  print("b vector after Gaussian forward elimination")
  print(y)
  as <- c()
  as[deg+1] <- y[deg+1]/mat[deg+1,deg+1]
  for(i in deg:1){
    sum <- y[i]
    for(j in (i+1):(deg+1)){
      sum <- sum - mat[i,j]*as[j]
    }
    as[i] <- sum/mat[i,i]
  }
  print("identified polynomial coefficients")
  print(as)
  return(as)
}

get_y <- function(as, x){
  deg <- length(as)-1
  y = as[1]
  for(degs in 1:deg){
    y <- y + as[degs+1]*x^degs
  }
  return(y)
}
