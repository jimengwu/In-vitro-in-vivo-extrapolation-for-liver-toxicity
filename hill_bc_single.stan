data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y;

  real<lower=1e-12> a_mean;
  real<lower=1e-12> a_sd;
  real<lower=1e-12> b_mean;
  real<lower=1e-12> b_sd;
  real<lower=1e-12> d_mean;
  real<lower=1e-12> d_sd;
  real<lower=1e-12> c_mean;
  real<lower=1e-12> c_sd;
}

parameters {
  real<lower=0> a;
  real<lower=0> b;
  real<lower=0> d;
  real<lower=0> c;
  real<lower=1e-12> sigma;
}

transformed parameters {
  vector[N] mu;

  for (i in 1:N) {
    real denom = fmax(1e-12, pow(b, d) + pow(x[i], d));  // prevent zero
    //mu[i] = a * pow(c, pow(x[i], d) / (pow(b, d) + pow(x[i], d)));
    mu[i] = a * (1 + (c - 1) * pow(x[i], d) /denom);
  }
}

model {
  
  a ~ lognormal(log(a_mean), a_sd);
  b ~ lognormal(log(b_mean), b_sd);
  d ~ lognormal(log(d_mean), d_sd);
  c ~ lognormal(log(c_mean), c_sd);
  sigma ~ normal(0, 0.3);

  y ~ normal(mu, sigma);
}

generated quantities {
  real ec5;
  if (d > 1e-6)
    ec5 = b * pow(1.0 / 19.0, 1.0 / d);
  else
    ec5 = negative_infinity();  // clearly flagged as invalid
  real ec10;
  if (d > 1e-6)
    ec10 = b * pow(1.0 / 9.0, 1.0 / d);
  else
    ec10 = negative_infinity();  // clearly flagged as invalid
    
}
