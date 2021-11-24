#### Nonparametric Monotonic *(npM)*

This nonparametric model is fit over the mean responses at each dose or the medians and it is analogous to *npS*, with the difference of forcing a non-increasing constrain on the connected linear functions.

#### Details

This alternative to the nonparametric spline imposes a monotone (non increasing) constrain. If the spline between two doses does not have a negative slope, then the average of the previous doses is recursively calculated until the next spline has a negative slope, which connects previously calculaed average with the next data point. It is also called isotonic regression in some literature. 

<img src="images/fig4.png" alt="drawing" style="width:550px;"/>



where  *k* is the smallest integer such that the slope is negative. It has been shown that this fit is the one that minimzed the squared errors $\|y-f(x)\|^2$ with the constrain that $f(x_i)\geq f(x_{i+1})$.

This fit will always produce a non increasing fit, even in cases where the logistic model or the nonparametric Bayesian produce an incrasing fit.

For the calculation of the effective doses at level $p\in (0,1)$ we find the unique value in the x-axis such that 
$$
f(IC_{p}) = \min_{x_1\leq x\leq x_n}f(x)+p\left(\max_{x_1\leq x\leq x_n}{f(x)}-\min_{x_1\leq x\leq x_n}{f(x)}\right).
$$

