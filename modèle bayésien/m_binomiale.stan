// Title : m_binomial
// Auteur : Maud OTELLO
// Date : 17.04.2024

// Donnees
data {
  int <lower = 0> N ; //nombre d'observation
  int <lower = 0, upper = 1> y[N]; //vecteur de présence abscence
}

// Nous disposons d'un parametre p.
parameters {
  real<lower = 0, upper = 1 > p;
}

// The model to be estimated. Notre model 'y' suit une loi binomiale.
model {
  y ~ binomial (1, p); //likehood
  p ~ uniform (0,1); // priors
}

