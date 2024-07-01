// 
// Title : m_DBH_pred
// Autor : Maud OTELLO
// Date : 14/06/2024
// 


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
  real <lower = 0 > X[N] ; // Variable environnementale (Surface terrière)
  vector[S] alpha_s ; // moyenne des priors sur les alpha du modèle 1
  int <lower = 1> N_pred ; // longueur de la prédiction
  real <lower = 0 > X_pred[N_pred] ; // prédiction sur des séquences
}

transformed data {
  real X_cr[N] ; 
  real X_pred_cr[N_pred];
  for (n in 1:N) {
    X_cr[n] = (log(X[n])- mean(log(X)))/sd(log(X)); // centrée & réduit
  }
  for (n in 1:N_pred) {
        X_pred_cr[n] = (log(X_pred[n])-mean(log(X)))/sd(log(X)); // prédiction centré & réduit
  }
}

parameters {
  vector[S] alpha ; // origine à l'ordonnée (par espèce)
  vector[S] beta ; // coefficient de DBH
  vector[S] gamma ; // coefficient de DBH²
}

model {
  alpha ~ normal(alpha_s,1); // Intercept
  beta ~ normal (0,1) ; // prior coefficient de DBH
  gamma ~ normal(0,1) ; // prior coefficient de DBH²
  
  for (n in 1:N)
  Y[n] ~ dirichlet_multinomial(softmax(alpha + beta*X_cr[n] + gamma*(X_cr[n])^2)); // likehood
}

generated quantities {
  matrix [N_pred, S] Y_pred ;
  for (n in 1:N_pred)
  Y_pred[n]= to_row_vector(softmax(alpha + beta*X_pred_cr[n] + gamma*(X_pred_cr[n])^2)); // prédiction
}
