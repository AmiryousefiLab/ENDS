#### Nonparametric Bayesian Model *(NPB)*

 Given points $(x_i, y_i)$ for $i=1..n$, with knots $K=\{k_i\}_{i=1}^J$ and parameter $\lambda>0$ we will model the relation between *x* and *y* with the Bayesian Hierarchical model

<img src="images/fig5.png" alt="drawing" style="width:400px;"/>

Weights *a*​ follow are non-negative and sum up to one. We can use MCMC to sample from the posterior distribution.​
$$
p( C ,a,\sigma^{2} \mid y,x,K,\lambda)
$$
 An inhouse implementation of the Metropolis-Hastings (MH) sampling algorithm is used for posterior inference. Knots are chosen at the unique doses and $\lambda$​​​​​ is chosen as the squared error minimizer over a grid of values including the mean variance estimate of the samples at each knot. We also chose a uniform dirichlet distribution as prior for *a* this is $a\sim Dir(1,..,1)$, slightly different model choice from ["A Bayesian Monotonic Non-parametric Dose-Response Model"](https://www.tandfonline.com/doi/abs/10.1080/10807039.2021.1956298) where a stick breaking process is used to generate the weights of the dirichlet distribution.

