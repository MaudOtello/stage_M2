// Title : m_centre_reduit
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
  ////////////////// Variables prédites à calculer ///////////
  int <lower = 1> N_pred ; // longueur de la prédiction
  real <lower = 0 > X_pred[N_pred] ; // prédiction sur des séquences
  int <lower = 1> K_pred ; // longueur de la prédiction
  real <lower = 0 > Z_pred[K_pred] ; // prédiction sur des séquences
  
}

transformed data {
  real X_cr[N] ;
  real Z_cr[N] ; 
  real Z_pred_cr[K_pred] ;
  real X_pred_cr[N_pred] ;
  
  // Variables explicatives centrées réduites
  for (n in 1:N) {
    X_cr[n] = (log(X[n])- mean(log(X)))/sd(log(X)); // log sur les DBH, centrée et réduit
    Z_cr[n] = (Z[n]-mean(Z))/sd(Z) ; // centrée & réduit
  }
  
  // Variables explicatives prédites centrées réduites aussi
  for (n in 1:N_pred) {
    X_pred_cr[n] = (log(X_pred[n])-mean(log(X)))/sd(log(X)); // prédiction centré & réduit
  }
  
  for (k in 1:K_pred){
    Z_pred_cr[k] = (Z_pred[k]-mean(Z))/sd(Z);
  }
}
parameters {
  vector[S] alpha ; // origine à l'ordonnée (par espèce)
  vector[S] beta ; // coeficient du DBH (par espèce)
  vector[S] gamma ; // coefficient du DBH² (par espèce)
  vector[S] delta ; // coefficient de G (par espèce)
  vector[S] epsilon ; // coefficient de G² (par espèce)
}

model {
  
  // priors
  alpha ~ normal(alpha_s,1); 
  beta ~ normal(0,1); 
  gamma ~ normal(0,1);
  delta ~ normal(0,1); 
  epsilon ~ normal(0,1);
  
  // modèle
  for (n in 1:N)
  Y[n] ~ dirichlet_multinomial(softmax(alpha + beta*X_cr[n] + gamma*(X_cr[n])^2 + delta*Z_cr[n]+ epsilon*(Z_cr[n])^2)); // likehood
}

generated quantities {
  // Prédiction pour le DBH
  matrix [N_pred, S] Y_pred_x ; // format de la sortie
  for (n in 1:N_pred){
    Y_pred_x[n]= to_row_vector(softmax(alpha + beta*X_pred_cr[n] + gamma*(X_pred_cr[n])^2 )); // prédiction de Y quand X varie (DBH)   
    }
  // Prédiction pour la ST
  matrix [K_pred, S] Y_pred_z ;
  for (k in 1:K_pred){
      Y_pred_z[k]= to_row_vector(softmax(alpha + delta*Z_pred_cr[k]+ epsilon*(Z_pred_cr[k])^2)); // prédiction de Y quand Z varie (ST)
    }
}
