//Titre  : m_plateau
//Autrice : Maud OTELLO
//Date : 24/04/2024
data {
  int <lower = 0> N ; // Nombre d'observation (ligne)
  int <lower = 1> S ; // Nombre d'espèce (colonne)
  int <lower = 0 , upper = 1 > y[N,S]; //Matrice présence absence par espèce et observation
  int <lower = 0 , upper = 1 > z[N] ;  //Topographie, présence en plateau
}

parameters {
  simplex[S] alpha ;              // Paramettre de la multinomiale_diritchelt
  }

model {
  for (n in 1:N){
    if (z[n]== 1){ // précise en présence de plateau
    y[n] ~ multinomial (alpha) ; //likehood
  } else{}
  }
  alpha ~ dirichlet(rep_vector(0,S)) ;    // prior
}
