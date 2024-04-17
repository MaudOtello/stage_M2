data {
  int<lower=0> N;
  int y[N] ;
}

// Nous disposons d'un parametre p.
parameters {
  real<lower = 0, upper = 1 > p;
}

// The model to be estimated. Notre model 'y' suit une loi binomiale.
model {
  y ~ binomial(N, p); //likehood
  p ~ uniform (0,1); // priors
}


