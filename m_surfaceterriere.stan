// Title : m_surfaceterriere
// Autor : Maud OTELLO
// Date : 13/06/2024


functions {
  real dirichlet_multinomial_lpmf(int[] y, vector alpha) {
    real alpha_plus = sum(alpha);
    return lgamma(alpha_plus) + sum(lgamma(alpha + to_vector(y)))
    - lgamma(alpha_plus+sum(y)) - sum(lgamma(alpha));
  }
}

data {
  int <lower = 0> N ; // Nombre d'observation (ligne)
  int <lower = 1> S ; // Nombre d'espèce (colonne)
  int <lower = 0 , upper = 1 > Y[N,S]; // Matrice présence absence par espèce et observation
  real <lower = 0 > Z[N] ; // Variable environnementale (Surface terrière)
  vector[S] alpha_s ; // moyenne des priors sur les alpha du modèle 1
  int <lower = 1> N_pred ; // longueur de la prédiction
  real <lower = 0 > Z_pred[N_pred] ; // prédiction sur des séquences
}

transformed data {
  real Z_cr[N] ; 
  real Z_pred_cr[N_pred];
  for (n in 1:N) {
    Z_cr[n] = (Z[n]-mean(Z))/sd(Z) ; // centrée & réduit
  }
  for (n in 1:N_pred) {
        Z_pred_cr[n] = (Z_pred[n]-mean(Z))/sd(Z); // prédiction centré & réduit
  }
}

parameters {
  vector[S] alpha ; // origine à l'ordonnée (par espèce)
  vector[S] delta ; // coefficient de G
  vector[S] epsilon ; // coefficient de G²
}

model {
  alpha ~ normal(alpha_s,1); // prior coefficient de G
  delta ~ normal(0,1) ; // prior coefficient de G
  epsilon ~ normal(0,1); // prior coefficient G²
  
  for (n in 1:N)
  Y[n] ~ dirichlet_multinomial(softmax(alpha + delta*Z_cr[n]+ epsilon*(Z_cr[n])^2)); // likehood
}

generated quantities {
  matrix [N_pred, S] Y_pred ;
  for (n in 1:N_pred)
  Y_pred[n]= to_row_vector(softmax(alpha + delta*Z_pred_cr[n]+ epsilon*(Z_pred_cr[n])^2)); // prédiction
}
