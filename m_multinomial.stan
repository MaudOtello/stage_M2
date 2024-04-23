//Titre  : m_multinomial
//Autrice : Maud OTELLO
//Date : 23/04/2024
data {
  int <lower = 0> N ; // Nombre d'observation (ligne)
  int <lower = 1> S ; // Nombre d'espèce (colonne)
  int <lower = 0 , upper = 1 > y[N,S]; //Contenue de la matrice
}

parameters {
  simplex[S] theta ;   // Paramettre de la multinomiale
}

model {
  for (n in 1:N){
    y[n] ~ multinomial(theta) ; //likehood
  }
  theta ~ uniform(0,1)    ;    // prior
}

