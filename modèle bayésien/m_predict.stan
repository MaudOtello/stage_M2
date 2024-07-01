// Title : m_predict
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
  real <lower = 0 > X[N] ; // Variable explicative (DBH)
  real <lower = 0 > Z[N] ; // Variable environnementale (Surface terrière)
  vector[S] alpha_s ; // moyenne des priors sur les alpha du modèle 1
  int <lower = 1> N_pred ; // longueur de la prédiction
  real <lower = 0 > Z_pred[N_pred] ; // prédiction sur des séquences
}

transformed data {
  real X_cr[N] ;
  real Z_cr[N] ; 
  real Z_pred_cr[N_pred];
  for (n in 1:N) {
    X_cr[n] = (log(X[n])- mean(log(X)))/sd(X); // log sur les DBH, centrée et réduit
    Z_cr[n] = (Z[n]-mean(Z))/sd(Z) ; // centrée & réduit
  }
  for (n in 1:N_pred) {
        Z_pred_cr[n] = (Z[n]-mean(Z))/sd(Z); // prédiction centré & réduit
  }
}

parameters {
  vector[S] alpha ; // origine à l'ordonnée (par espèce)
  vector[S] beta ; // coeficien du DBH (par espèce)
  vector[S] gamma ; // coefficien du DBH²
  vector[S] delta ; // coefficien de G
}

model {
  alpha ~ normal(alpha_s,1); // prior
  beta ~ normal(0,1); // prior
  gamma ~ normal(0,1); // prior 
  delta ~ normal(0,1) ; // prior
  
  for (n in 1:N)
  Y[n] ~ dirichlet_multinomial(softmax(alpha + beta*X_cr[n] + gamma *(X_cr[n])^2 + delta*Z_cr[n])); // likehood
}

generated quantities {
  matrix [N_pred, S] Y_pred ;
  for (n in 1:N_pred)
  Y_pred[n]= to_row_vector(softmax(alpha + beta*X_cr[n] + gamma *(X_cr[n])^2 + delta*Z_pred_cr[n])); // prédiction
}
