// Title : m_varcontinue
// Autor : Maud OTELLO
// Date : 13/05/24

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
  real <lower = 0 > X[N] ; // Variable environnementale
}

parameters {
  vector[S] alpha ; // origine à l'ordonnée (par espèce)
  vector[S] beta ; // coeficien de la courbe (par espèce)
}

model {
  alpha ~ normal(0,1); // prior
  beta ~ normal(0,1); // prior
  
  for (n in 1:N)
   Y[n] ~ dirichlet_multinomial(softmax(alpha + beta*X[n])); // likehood
  
}
