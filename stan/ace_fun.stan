// Adapted from https://mc-stan.org/docs/stan-users-guide/custom-probability.html
real bi_normal_cdf(real z1, real z2, real rho) {
  if (z1 == 0 && z2 == 0) {
    return 0.25 + asin(rho) / (2 * pi());
  }
  real denom = sqrt((1 + rho) * (1 - rho));
  real term1 = z1 == 0
    ? (z2 > 0 ? 0.25 : -0.25)
    :  owens_t(z1, (z2 / z1 - rho) / denom);
  real term2 = z2 == 0
    ? (z1 > 0 ? 0.25 : -0.25)
    : owens_t(z2, (z1 / z2 - rho) / denom);
  real z1z2 = z1 * z2;
  real delta = z1z2 < 0 || (z1z2 == 0 && (z1 + z2) < 0);
  return 0.5 * (Phi(z1) + Phi(z2) - delta) - term1 - term2;
}
matrix polycor_to_prob(vector b, real rho) {
  int k = size(b) + 1;
  matrix[k, k] prob;
  matrix[k, k] cprob;
  for (j in 1:k) {
    for (i in j:k) {
      if (i == k && j == k) {
        cprob[i, j] = 1.0;
      } else if (i == k) {
        cprob[i, j] = Phi(b[j]);
      } else {
        cprob[i, j] = bi_normal_cdf(b[i] | b[j], rho);
      }
      if (i != j) {
        cprob[j, i] = cprob[i, j];
      }
      if (j == 1 && i == 1) {
        prob[i, j] = cprob[i, j];
      } else if (j == 1) {
        prob[i, j] = cprob[i, j] - cprob[i - 1, j];
      } else {
        prob[i, j] = cprob[i, j] - cprob[i - 1, j] - cprob[i, j - 1] +
                       cprob[i - 1, j - 1];
      }
      if (i != j) {
        prob[j, i] = prob[i, j];
      }
    }
  }
  return prob;
}
vector vec_symm(matrix m) {
  int k = rows(m);
  vector[k * (k + 1) %/% 2] v;
  int idx = 1;
  for (j in 1:k) {
    for (i in j:k) {
      if (i == j) {
        v[idx] = m[i, j];
      } else {
        v[idx] = m[i, j] * 2;
      }
      idx += 1;
    }
  }
  return v;
}
