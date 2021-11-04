// Stick breaking priorr
data {
  int<lower=0> n;
  int<lower=0> n_knots;
  real x[n];
  real y[n];
  real K[n_knots];
  real<lower=0> lambda;
  real<lower=0> alpha;  // Parameter for stickbreaking process, each componentn Beta(1,alpha)
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real <lower=0 , upper = 1> C; 
  real <lower=0.0> sigma2;
  vector<lower=0,upper=1>[n_knots - 1] v;  // stickbreak components
}

transformed parameters {
  simplex[n_knots] a;  // stick breaking process
  vector<lower=0,upper=1>[n_knots - 1] cumprod_one_minus_v;
  cumprod_one_minus_v = exp(cumulative_sum(log1m(v)));
  a[1] = v[1];
  a[2:(n_knots-1)] = v[2:(n_knots-1)] .* cumprod_one_minus_v[1:(n_knots-2)];
  a[n_knots] = cumprod_one_minus_v[n_knots - 1];
}


model {
  C ~ beta(1,1);
  sigma2 ~ chi_square(2);
  v ~ beta(1, alpha);
  for(i in 1:n){
      y[i] ~ normal(
        C + (1-C)*sum( a * (1-normal_cdf(rep_vector(x[i], n_knots), K, rep_vector(sqrt(lambda), n_knots) ) ) )
        ,sqrt(sigma2)
        );
  }
}

