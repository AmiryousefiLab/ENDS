#### Parametric Logistic *(PL)*

This is the usual four parameter logistic fit to the data by mean square errror. Using the equation 
$$
f(x) = c + \frac{d-c}{1+\exp{b(\log x - \log e)}}
$$
*c* is the assympotic minimun value, *d* is the assympotic minimum value and *e* is the $IC_{50}$. Note that we are using   $\log x$ in the formula since the fit is done over the logarithm of the doses. The parameters are estimated by minimizing the squared error function $\|y_i-f(x_i)\|^2$. There is not an analytical soluiton to this problem, different numerical minimization algorithms exisit for estimating the parameters. 

