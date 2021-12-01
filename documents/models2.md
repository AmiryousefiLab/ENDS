#### Parametric Logistic *(pL)*

This model is the logistic function adjusted to the data with four parameters found by least square estimation. The model assumes a fixed "S" shape decreasing curve for the responses.  

#### Details

This is the usual four parameter logistic fit to the data to minimize a squared error loss function, the model fit is the logistif function given by
$$
f(x) = c + \frac{d-c}{1+\exp{b(\log x - \log e)}},
$$
where *c* is the assympotic minimun value, *d* is the assympotic maximum value and *e* is the $IC_{50}$. Note that we are using  $\log x$ in the formula since the fit is done over the logarithm of the doses. The parameters are estimated by minimizing the squared error function $\|y-f(x)\|^2$. There is not an analytical soluiton to this problem, different numerical minimization algorithms exisit for estimating the parameters, our implementation uses the **R** package *drc*, which internally uses the base **R** optimizer *optim*. Note that the function is monotonic decreasing unless we have a degenerate fit, in which the adjusted function is a linear fit with non negative slope. 

Since the function is monotonic and continuous then for any value in the interval $p \in (0,1)$  the $IC_p$ is the unique value wich maps $f(IC_{p}) = c+p(d-c)$.

