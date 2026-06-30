//
// This Stan program estimates the polychoric correlation between two
// ordinal variables, and decompose based on an ACE model. This model
// assumes symmetry in the contingency table.
//
functions {
  #include "ace_fun.stan"
}
data {
  int<lower=0> k;  // number of categories
  array[k * (k + 1) %/% 2] int<lower=0> ltri_ns1;  // mz
  array[k * (k + 1) %/% 2] int<lower=0> ltri_ns2;  // dz
}
transformed data {
  int kprime = k * (k + 1) %/% 2;
}
parameters {
  real<lower=-1, upper=1> rho1;
  real<lower=-1, upper=1> rho2;
  ordered[k - 1] thres;
}
transformed parameters {
  vector<lower=0>[kprime] theta1 =
    vec_symm(polycor_to_prob(thres, rho1));
  vector<lower=0>[kprime] theta2 =
    vec_symm(polycor_to_prob(thres, rho2));
}
// The model to be estimated.
model {
  thres ~ normal(0, 1.5);
  ltri_ns1 ~ multinomial(theta1);
  ltri_ns2 ~ multinomial(theta2);
}
