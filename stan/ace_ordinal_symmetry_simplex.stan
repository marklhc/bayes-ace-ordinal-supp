//
// This Stan program estimates the polychoric correlation between two
// ordinal variables, and decompose based on an ACE model. The model is
// parameterized using a simplex for the ACE components to enforce
// nonzero variances, and assumes symmetry in the contingency table.
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
  simplex[3] ace;
  ordered[k - 1] thres;
}
transformed parameters {
  real rho1 = ace[1] + ace[2];
  real rho2 = ace[1] / 2 + ace[2];
  vector<lower=0>[kprime] theta1 =
    vec_symm(polycor_to_prob(thres, rho1));
  vector<lower=0>[kprime] theta2 =
    vec_symm(polycor_to_prob(thres, rho2));
}
// The model to be estimated.
model {
  thres ~ student_t(4, 0, 2.5);
  ltri_ns1 ~ multinomial(theta1);
  ltri_ns2 ~ multinomial(theta2);
}
generated quantities {
  real h2 = ace[1];
  real c2 = ace[2];
  real e2 = ace[3];
  array[kprime] int ltri_ns1_tilde = multinomial_rng(theta1, N1);
  array[kprime] int ltri_ns2_tilde = multinomial_rng(theta2, N2);
  real X2_1_samp = 0.0;
  real X2_2_samp = 0.0;
  int<lower=0, upper=1> X2_1_gt;
  int<lower=0, upper=1> X2_2_gt;
  {
    vector[kprime] N1theta = N1real .* theta1;
    vector[kprime] N2theta = N2real .* theta2;
    real X2_1_rep = 0.0;
    real X2_2_rep = 0.0;
    for (idx in 1:kprime) {
      real ns1_exp = N1theta[idx];
      real ns2_exp = N2theta[idx];
      X2_1_samp += (ltri_ns1[idx] - ns1_exp)^2 / ns1_exp;
      X2_2_samp += (ltri_ns2[idx] - ns2_exp)^2 / ns2_exp;
      X2_1_rep += (ltri_ns1_tilde[idx] - ns1_exp)^2 / ns1_exp;
      X2_2_rep += (ltri_ns2_tilde[idx] - ns2_exp)^2 / ns2_exp;
    }
    X2_1_gt = X2_1_rep > X2_1_samp;
    X2_2_gt = X2_2_rep > X2_2_samp;
  }
}
