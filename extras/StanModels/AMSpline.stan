
data {
  int<lower=0> n;
  int<lower=0> m;
  int<lower=0> n_knots;
  real x[n];
  matrix [n, m] y;
  real K[n_knots];
  real<lower=0> lambda;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real <lower=0 , upper = 1> C; 
  real <lower=0.0> sigma2;
  simplex[n_knots] a;
  //vector<lower=0.0, upper=1>[n_knots] a;
}

model {
  C ~ beta(1,1);
  sigma2 ~ chi_square(2);
  a ~ dirichlet(rep_vector(1.0,n_knots));
  for(i in 1:n){
    for(j in 1:m){
      y[i,j] ~ normal(
        C + (1-C)*sum( a*(1-normal_cdf(rep_vector(x[i],n_knots) ,K, rep_vector(sqrt(lambda), n_knots) ) ) )
        ,sqrt(sigma2)
        );
    }
  }
}

