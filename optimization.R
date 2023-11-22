newton_solve <- function(f, df, x0, eps=1e-10, max_iter=10) {
  k <- 1
  
  while (k < max_iter & abs(f(x0)) > eps) {
    # Check if the derivative at the point is non-zero
    if (df(x0) != 0) {
      
      # Update estimate
      x1 <- x0 - f(x0) / df(x0)
      
      x0 <- x1
    }
    k <- k + 1
  }
  
  print(sprintf("Solution: %.10f, Error: %e, Iterations: %d", x1, abs(f(x1)), k))
  return(x1)
}

secant_solve <- function(f, x0, x1, eps=1e-10, max_iter=10) {
  k <- 1
  
  while (k < max_iter & abs(f(x1)) > eps) {
    # Update estimate
    x2 <- x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0))
    
    x0 <- x1
    x1 <- x2
    k <- k + 1
  }
  
  print(sprintf("Solution: %.10f, Error: %e, Iterations: %d", x1, abs(f(x1)), k))
  return(x1)
}

# CH4
T <- 180
P <- 4.5
R <- 8.314413
Tc <- 191.15
Pc <- 4.641
w <- 0.0115
k <- 0.37464+1.54226*w-0.26992*w^2
alpha <- (1 + k*(1 - sqrt(T/Tc)))^2
a <- (0.45724/Pc)*(R*Tc)^2
b <- (0.07780*R*Tc)/Pc

c1 <- b - R*T/P
c2 <- -3*b^2 + (a*alpha - 2*b*R*T)/P
c3 <- b^3 + (b^2 * R*T - a*b*alpha)/P

f <- function(Vm) Vm^3 + c1*Vm^2 + c2*Vm + c3

df <- function(Vm) 3*Vm^2 + 2*c1*Vm + c2

x0 <- R*T/P # Ideal gas solution

result <- newton_solve(f, df, x0)

result <- secant_solve(f, x0 - 10, x0 + 10)

# H2
T <- 273.15
P <- 101.325
R <- 8.314413
Tc <- 33.15
Pc <- 1298
w <- -0.22
k <- 0.37464+1.54226*w-0.26992*w^2
alpha <- (1 + k*(1 - sqrt(T/Tc)))^2
a <- (0.45724/Pc)*(R*Tc)^2
b <- (0.07780*R*Tc)/Pc

c1 <- b - R*T/P
c2 <- -3*b^2 + (a*alpha - 2*b*R*T)/P
c3 <- b^3 + (b^2 * R*T - a*b*alpha)/P

f1 <- function(Vm) Vm^3 + c1*Vm^2 + c2*Vm + c3

df <- function(Vm) 3*Vm^2 + 2*c1*Vm + c2

x0 <- R*T/P # Ideal gas solution

result <- newton_solve(f1, df, x0)

result <- secant_solve(f1, x0 - 10, x0 + 10)

# O2
T <- 273.15
P <- 101.325
R <- 8.314413
Tc <- 154.6
Pc <- 5050
w <- 0.022
k <- 0.37464+1.54226*w-0.26992*w^2
alpha <- (1 + k*(1 - sqrt(T/Tc)))^2
a <- (0.45724/Pc)*(R*Tc)^2
b <- (0.07780*R*Tc)/Pc

c1 <- b - R*T/P
c2 <- -3*b^2 + (a*alpha - 2*b*R*T)/P
c3 <- b^3 + (b^2 * R*T - a*b*alpha)/P

f2 <- function(Vm) Vm^3 + c1*Vm^2 + c2*Vm + c3

df <- function(Vm) 3*Vm^2 + 2*c1*Vm + c2

x0 <- R*T/P # Ideal gas solution

result <- newton_solve(f2, df, x0)

result <- secant_solve(f2, x0 - 10, x0 + 10)

# N2
T <- 273.15
P <- 101.325
R <- 8.314413
Tc <- 126.2
Pc <- 3390
w <- 0.04
k <- 0.37464+1.54226*w-0.26992*w^2
alpha <- (1 + k*(1 - sqrt(T/Tc)))^2
a <- (0.45724/Pc)*(R*Tc)^2
b <- (0.07780*R*Tc)/Pc

c1 <- b - R*T/P
c2 <- -3*b^2 + (a*alpha - 2*b*R*T)/P
c3 <- b^3 + (b^2 * R*T - a*b*alpha)/P

f3 <- function(Vm) Vm^3 + c1*Vm^2 + c2*Vm + c3

df <- function(Vm) 3*Vm^2 + 2*c1*Vm + c2

x0 <- R*T/P # Ideal gas solution

result <- newton_solve(f3, df, x0)

result <- secant_solve(f3, x0 - 10, x0 + 10)



# O2
T <- 273.15 + 550
P <- 5000
R <- 8.314413
Tc <- 154.6
Pc <- 5050
w <- 0.022
k <- 0.37464+1.54226*w-0.26992*w^2
alpha <- (1 + k*(1 - sqrt(T/Tc)))^2
a <- (0.45724/Pc)*(R*Tc)^2
b <- (0.07780*R*Tc)/Pc

c14 <- b - R*T/P
c24 <- -3*b^2 + (a*alpha - 2*b*R*T)/P
c34 <- b^3 + (b^2 * R*T - a*b*alpha)/P

f4 <- function(Vm) Vm^3 + c14*Vm^2 + c24*Vm + c34

df <- function(Vm) 3*Vm^2 + 2*c1*Vm + c2

x0 <- R*T/P # Ideal gas solution

result <- newton_solve(f4, df, x0)

result <- secant_solve(f4, x0 - 10, x0 + 10)

curve(f1, xlim=c(-10,10), col='black', lwd=2, lty=1, ylim=c(-500, 500), xlab="Vm", ylab="Resíduo")
curve(f4, add=TRUE, col='black', lty = 2, lwd=2)
abline(h=0, lty=2, col="lightgrey")
abline(v=0, lty=2, col="lightgrey")
legend("topleft", legend = c("O2 em condições normais", "O2 em condições extremas"), lty = c(1,2))
