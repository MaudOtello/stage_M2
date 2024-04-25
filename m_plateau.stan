//Titre  : m_plateau
//Autrice : Maud OTELLO
//Date : 24/04/2024
data {
  int <lower = 0> N ; // Nombre d'observation (ligne)
  int <lower = 1> S ; // Nombre d'espèce (colonne)
  int <lower = 0 , upper = 1 > y[N,S]; //Matrice présence absence par espèce et observation
  int <lower = 0 , upper = 1 > x[N] ;  //Topographie, présence en plateau
}

parameters {
  simplex[S] alpha ;              // Paramettre de la multinomiale_diritchelt
  }

model {
  for (n in 1:N){
    y[n] ~ dirichlet_multinomial()(alpha) ; //likehood
  }
  alpha ~ uniform(0,1) ;    // prior
}
