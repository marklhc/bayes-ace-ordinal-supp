//
// This Stan program estimates the polychoric correlation between two
// ordinal variables, and decompose based on an ACE model. This model
// assumes symmetry in the contingency table.
//
functions {
  #include "ace_fun.stan"
}
data {
  int<lower=2> k;  // number of categories
  array[k * (k + 1) %/% 2] int<lower=0> ltri_ns1;  // mz
  array[k * (k + 1) %/% 2] int<lower=0> ltri_ns2;  // dz
}
transformed data {
  int kprime = k * (k + 1) %/% 2;
}
parameters {
  real<lower=-1, upper=1> rho1;
  real<lower=-1, upper=1> rho2;
  // Match blavaan's prior on the first threshold and log spacings.
  vector[k - 1] thres_un;
}
transformed parameters {
  ordered[k - 1] thres;
  vector<lower=0>[kprime] theta1;
  vector<lower=0>[kprime] theta2;
  thres[1] = thres_un[1];
  if (k > 2) {
    for (j in 2:(k - 1)) {
      thres[j] = thres[j - 1] + exp(thres_un[j]);
    }
  }
  theta1 = vec_symm(polycor_to_prob(thres, rho1));
  theta2 = vec_symm(polycor_to_prob(thres, rho2));
}
// The model to be estimated.
model {
  thres_un ~ normal(0, 1.5);
  ltri_ns1 ~ multinomial(theta1);
  ltri_ns2 ~ multinomial(theta2);
}
