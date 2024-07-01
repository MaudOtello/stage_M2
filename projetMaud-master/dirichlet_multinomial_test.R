library(MGLM) # for Dirichlet-Multinomial
library(rstan) # for parameters distributions estimation

# softmax function:
softmax <- function(x) exp(x) / sum(exp(x))

ourfun  <- function(x) x / sum(x)

############### Data simulation
# Let say we have N trees of three species A, B and C
# and an environmental descriptor x (here sampled in a uniform distribution for the simulations)
# species A likes high values of x, B likes medium values of x and C likes low values of x

# trees
N = 50

# environmental variable
x <- runif(N, -10, 10)

# parameters of the model for species A
a_A <- 10
b_A <- 3
# linear combination of parameter and environmental variable for species A
t_A <- b_A * x + a_A

# same for B
a_B <- 10
b_B <- 6
t_B <- b_B * x + a_B 

# same for C
a_C <- 10
b_C <- -3
t_C <- b_C * x + a_C

# aggragated in a table
t <- rbind(t_A, t_B, t_C)
# apply softmax to each column
s <- apply(t, 2, softmax)

## Generate Dirichlet-Multinomial data
# example : sample one time, one tree with parameter c(1,1,50) : rdirmn(n = 1, size = 1, alpha = c(1,1,50))
Y <- rdirmn(n = N, size = rep(1,N), alpha = t(s))

# test without softmax
# i <- apply(t, 2, ourfun)
# Y <- rdirmn(n = N, size = rep(1,N), alpha = t(i))

####################### Parameters estimation
# now let's try to estimate the parameters
fit = stan(file = "DM1.stan", 
           data=list(N = N, S = 3, Y = Y, X = x), 
           iter = 1000)
plot(fit)
traceplot(fit)



