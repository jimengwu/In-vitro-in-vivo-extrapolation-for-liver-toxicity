data {
  int<lower=0> N;        // number of observations
  vector[N] x;           // predictor (concentration)
  vector[N] y;           // observed response
  
  real a_mean;
  real<lower=0> a_sd;
  real b_mean;
  real<lower=0> b_sd;
  real d_mean;
  real<lower=0> d_sd;
}

parameters {
  real<lower=0, upper=1> a;
  real log_b;
  real log_d;
  real<lower=0> sigma;
}

transformed parameters {
  real b = exp(log_b);
  real d = exp(log_d);
  vector[N] mu;

  for (i in 1:N) {
    mu[i] = a + (1 - a) * (pow(x[i], d) / (pow(b, d) + pow(x[i], d)));
  }
}

model {
  a ~ normal(a_mean, a_sd);
  log_b ~ normal(log(b_mean), b_sd);
  log_d ~ normal(log(d_mean), d_sd);
  sigma ~ normal(0, 0.1);

  y ~ normal(mu, sigma);
}

generated quantities {
  real ec5;
  {
    // EC5: concentration at 5% effect relative to range from a to 1
    real target_resp = a + 0.05 * (1 - a);
    real numerator = target_resp - a;
    real denominator = 1 - a;

    if (numerator <= 0 || denominator <= 0 || numerator / denominator >= 1) {
      ec5 = 0.00001;
    } else {
      real frac = numerator / denominator;
      ec5 = b * pow(frac / (1 - frac), 1 / d);
    }
  }
  real ec10;
  {
    // EC10: concentration at 10% effect relative to range from a to 1
    real target_resp = a + 0.1 * (1 - a);
    real numerator = target_resp - a;
    real denominator = 1 - a;

    if (numerator <= 0 || denominator <= 0 || numerator / denominator >= 1) {
      ec10 = 0.00001;
    } else {
      real frac = numerator / denominator;
      ec10 = b * pow(frac / (1 - frac), 1 / d);
    }
  }
}
