#### Nonparametric Spline

Given a sample of doses $d_1,..,d_n$ and response matrix $(y_{ij})$ for $i=1,..,n$ and $j=1,..,m$, then for each dose we can obtain the dose mean $y_i=\frac{1}{m}\sum_j y_{ij}$ or dose medians​​. The simple spline that connects each of this points with a linear function is given by
$$
f(x) = \left(\frac{y_{i+1}-y_i}{x_{i+1}-x_i}\right) x + y_i \quad \textrm{for } x_i\leq x\leq x_{i+1}
$$
