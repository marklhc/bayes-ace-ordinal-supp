//
// Unrestricted-correlation ACE model with a standard-normal threshold prior,
// matching the threshold prior used by the simplex simulation models.
// Otherwise identical to ace_ordinal_symmetry_no_x2.stan.
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
  int N1 = sum(ltri_ns1);
  real N1real = N1;
  int N2 = sum(ltri_ns2);
  real N2real = N2;
}
parameters {
  real<lower=0, upper=1> scaled_rho1;
  real<lower=0, upper=1> scaled_rho2;
  ordered[k - 1] thres;
}
transformed parameters {
  real rho1 = scaled_rho1 * 2 - 1;
  real rho2 = scaled_rho2 * 2 - 1;
  vector<lower=0>[kprime] theta1 =
    vec_symm(polycor_to_prob(thres, rho1));
  vector<lower=0>[kprime] theta2 =
    vec_symm(polycor_to_prob(thres, rho2));
}
// The model to be estimated.
model {
  scaled_rho1 ~ beta(1.5, 1.5);
  scaled_rho2 ~ beta(1.5, 1.5);
  thres ~ std_normal();
  ltri_ns1 ~ multinomial(theta1);
  ltri_ns2 ~ multinomial(theta2);
}
generated quantities {
  real h2 = 2 * (rho1 - rho2);
  real c2 = 2 * rho2 - rho1;
}
