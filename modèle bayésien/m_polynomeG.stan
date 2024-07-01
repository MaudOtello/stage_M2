// Title : m_polynomeG
// Autor : Maud OTELLO
// Date : 

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
  vector[S] alpha_s ; // moyenne des priors sur les alpha modèle 1
}

transformed data {
  real X_std[N] ;
  for (n in 1:N) {
    X_std[n] = log(X[n]); // log sur les DBH
  }
}
parameters {
  vector[S] alpha ; // origine à l'ordonnée (par espèce)
  vector[S] beta ; // coeficien du DBH (par espèce)
  vector[S] gamma ; // coefficien du DBH²
  vector[S] epsilon ; // coefficien de G
}

model {
  alpha ~ normal(alpha_s,1); // prior
  beta ~ normal(0,1); // prior
  gamma ~ normal(0,1); // prior 
  epsilon ~ normal(0,1) ; // prior
  
  for (n in 1:N)
   Y[n] ~ dirichlet_multinomial(softmax(alpha + beta*X_std[n] + gamma *(X_std[n])^2 + epsilon*Z[n])); // likehood
  
}

