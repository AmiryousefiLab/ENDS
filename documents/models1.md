#### Nonparametric Spline *(npS)*

Given a sample of doses $d_1,..,d_n$ and for each i-th dose $d_i$ *m* responses $y_{i1},..,y_{im}$, then for each dose we can obtain the dose mean $y_i=\frac{1}{m}\sum_j y_{ij}$ or dose medians​​. The simple spline that connects each of this points with a linear function is given by the piece wise linear function
$$
f(x) = \left(\frac{y_{i+1}-y_i}{x_{i+1}-x_i}\right) x + y_i \quad \textrm{for } x_i\leq x\leq x_{i+1}.
$$
The function is defined on the interval $[x_1,x_n]$. Note that this function is not monotonic necessarily, if we define the $IC_{50}$ as the value in the x-axis such that $(\max y + \min y)/2 = f(IC_{50})$, then it might not be unique. Thus we calculate all the possible x-values that map to the halfway point of the responses and select as the effective dose 50% the one closest in absolute value to the one found by the monotonic fit. For any value $p\in (0,1)$ the $IC_{p}$ will be the chosen out of the candidates $(\max y + \min y) p = f(IC_{p})$ such that it is closest in distance to the uniquely one obtained by the monotone fit. 

In order to obtain the angles associated to the slope at each dose we have the i-th angle is $ \theta = \arctan \left(\frac{y_{i+1}-y_i}{x_{i+1}-x_i}\right)$.

It is clear that if the spline is obtained using the means the it will minimize the squared error $\frac{1}{nm}\sum_{ij} (y_{ij}-f(x_i))^2$, since the mean at each dose is the one that minimized this function, note that the median minimizes the absolute error, or L1 norm.  
