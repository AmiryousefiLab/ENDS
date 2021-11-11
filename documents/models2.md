#### Parametric Logistic *(pL)*

This is the usual four-parameter log logistic curve fitting to the data to minimize a squared error loss function. The model fit is the logistic function given by
$$
f(x) = c + \frac{d-c}{1+\exp{b(\log x - \log e)}}
$$
where *c* is the asymptotic minimum value, *d* is the asymptotic maximum value and *e* is the $IC_{50}$. Note that we are using   $\log x$ in the formula since the fit is done over the logarithm of the doses. The parameters are estimated by minimizing the squared error function $\|y-f(x)\|^2$. There is not an analytical soluiton to this problem, different numerical minimization algorithms exist for estimating the parameters, our implementation uses the **R** package *drc*, which internally uses the base **R** optimizer *optim*. Note that the function is monotonic decreasing unless we have a degenerate fit, in which the adjusted function is a linear fit with non negative slope. 

Since the function is monotonic and continuous then for any value in the interval $p \in (0,1)$  the $IC_p$ is the unique value which maps $f(IC_{p}) = (d-c)p$.

