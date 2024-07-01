// Title : m_binomial
// Auteur : Maud OTELLO
// Date : 17.04.2024

// Donnees
data {
  int<lower=0> n;         // nombre total d'observations
  int<lower=1,upper=2> y[n];  // observations (1 pour A, 2 pour B)
}
parameters {
  real<lower=0,upper=1> theta;  // probabilité pour l'espèce A
}
model {
  theta ~ beta(1, 1);  // a priori non informatif
  y ~ bernoulli(theta);  // modèle de Bernoulli
}

