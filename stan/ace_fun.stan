// Adapted from https://mc-stan.org/docs/stan-users-guide/custom-probability.html
real bi_normal_cdf(real z1, real z2, real rho) {
  real a = 1 / sqrt(1 - square(rho));
  real u = Phi(z1);
  real v = Phi(z2);
  real avg_uv = mean([u, v]);
  real alpha_u = a * (z2 / z1 - rho);
  real alpha_v = a * (z1 / z2 - rho);
  real d = 0;
  
  if (v == 0.5 || u == 0.5 && v != u) {
    real x = u == 0.5 ? z2 : z1;
    if (rho < 0) 
      return 0.5 * bi_normal_cdf(x | x, 1 - 2 * rho ^ 2);
    else 
      return u - 0.5 * bi_normal_cdf(x | x, 1 - 2 * rho ^ 2);
  }
  
  if (v == u) 
    return u - 2 * owens_t(inv_Phi(u), sqrt((1 - rho) / (1 + rho)));
  
  if (u < 0.5 && v >= 0.5 || u >= 0.5 && v < 0.5) 
    d = 0.5;
  
  return avg_uv - owens_t(z1, alpha_u) - owens_t(z2, alpha_v) - d;
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
