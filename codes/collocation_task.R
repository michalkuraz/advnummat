
# set your working directory to folder with R functions
setwd("./")

# links to second R scripts with functions
source('collocation_functions.R') 

# boundary condition
bc_points <- c(0,1,1,exp(-2)+sin(3)/13) # Order is: X1, Y1, X2, Y2

# x used for plotting
xs <- seq(0,1,0.001)

# degree of polynomials
deg <- 10

#defines collocation points
points <- seq(1/deg,(deg-1)/deg,length.out = deg-1)

# gets polynoam coefficiants
as <- get_polynom(deg, bc_points, points)

# evaluates the polynomial with the estimated coefficients and saves it to y
ys10 <- c()
for(i in 1:length(xs)){
  ys10[i] <- get_y(as, xs[i])
}

# See above
deg <- 5
points <- seq(1/deg,(deg-1)/deg,length.out = deg-1)
as <- get_polynom(deg, bc_points, points)
ys5 <- c()
for(i in 1:length(xs)){
  ys5[i] <- get_y(as, xs[i])
}

deg <- 4
points <- seq(1/deg,(deg-1)/deg,length.out = deg-1)
as <- get_polynom(deg, bc_points, points)
ys4 <- c()
for(i in 1:length(xs)){
  ys4[i] <- get_y(as, xs[i])
}

points <- c(1/3,2/3)
deg <- 3
as <- get_polynom(deg, bc_points, points)
ys3 <- c()
for(i in 1:length(xs)){
  ys3[i] <- get_y(as, xs[i])
}

deg <- 2
points <- seq(1/deg,(deg-1)/deg,length.out = deg-1)
as <- get_polynom(deg, bc_points, points)
ys2 <- c()
for(i in 1:length(xs)){
  ys2[i] <- get_y(as, xs[i])
}

pdf("collocation_adv_numerics.pdf")
plot(xs, ys10, col = "darkblue", type = "l", ylim = c(0,1), ylab = "y", xlab = "x", main = "Collocation method with different polynomial degrees")
par(new = T)
plot(xs, ys5, col = "royalblue4", type = "l", ylim = c(0,1), ylab = "y", xlab = "x")
par(new = T)
plot(xs, ys4, col = "royalblue3", type = "l", ylim = c(0,1), xaxt = "n", yaxt = "n", ylab = "", xlab = "")
par(new = T)
plot(xs, ys3, col = "royalblue2", type = "l", ylim = c(0,1), xaxt = "n", yaxt = "n", ylab = "", xlab = "")
par(new = T)
plot(xs, ys2, col = "lightblue", type = "l", ylim = c(0,1), xaxt = "n", yaxt = "n", ylab = "", xlab = "")
par(new = T)
plot(xs, exp(-2*xs)+sin(3*xs)/13, type = "l", col = "black", ylim = c(0,1), xaxt = "n", yaxt = "n", ylab = "", xlab = "")
leg.txt <- c("n = 10", paste("n =",5:2), "exp(-2x)+sin(3x)/13")
legend("bottom",leg.txt, col = c("darkblue","royalblue4","royalblue3","royalblue2", "lightblue", "black", "black"), lty = c(rep(1, 6),2), bty = "n", ncol =2)
dev.off()
