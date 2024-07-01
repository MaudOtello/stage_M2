functions {
  real dirichlet_multinomial_lpmf(int[] y, vector alpha) {
	real alpha_plus = sum(alpha);
	return lgamma(alpha_plus) + sum(lgamma(alpha + to_vector(y)))
            	- lgamma(alpha_plus+sum(y)) - sum(lgamma(alpha));
  }
}


// The input data is a vector 'y' of length 'N'.
data {
  int<lower = 1> N; // number of individuals
  int<lower = 1> S; // number of species
  int<lower = 0, upper=1> Y[N,S]; // individuals presence or absence for each species
  real X[N]; // environmental descriptor
}

// The parameters accepted by the model
parameters {
  vector[S] alpha ; // intercept (one per species)
  vector[S] beta ; // slope (one per species)
}


// The model to be estimated. We model the output
// 'y' in a Binomial distribution 
model {
  alpha ~ normal(0,10);  // prior
  beta ~ normal(0,10);   // prior
  
  for (n in 1:N)
	Y[n] ~ dirichlet_multinomial(softmax(alpha + beta*X[n])); // likelihood
}




