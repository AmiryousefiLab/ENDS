#### Nonparametric Spline *(npS)*

The model is the collectiion of linear functions that connect the average responses at each dose. 

#### Details

 Given a sample of doses $x_1,..,x_n$ where $x_i\leq x_{i+1}$ and for each i-th dose we have $m$ responses $y_{i1},..,y_{im}$, then for each dose we can obtain the dose mean $y_i=\frac{1}{m}\sum_j y_{ij}$ or alternatively we can calculate the dose medians. The simple spline that connects each of the means with a linear function is given by the piece wise linear function
$$
f(x) = \left(\frac{y_{i+1}-y_i}{x_{i+1}-x_i}\right) x + y_i \quad \textrm{for } x_i\leq x\leq x_{i+1}.
$$
The function is defined on the interval $x_1\leq x \leq x_n$. Note that this function is not monotonic necessarily, if we define the $IC_{50}$ as the value in the x-axis such that $f(IC_{50})=(\max y + \min y)/2$, then it might not be unique as it crosses the above function in multiple points, if the function is constant on an interval $I$ such that $f(I) = (\max y + \min y)/2$, then  IC$_{50}=\inf I$ is the infimum of the values on the interval. Thus we calculate all the possible x-values that map to the halfway point of the responses and select as the effective dose 50% the one closest in absolute value to the one found by the monotonic fit. For any value $p\in (0,1)$ the $IC_{p}$ will be the chosen out of the candidates $f(IC _{p}) = \min y + p(\max y - \min y)$, such that it is closest in distance to the uniquely value obtained by the monotone fit. Note that, if the function is constant we again choose as the IC$_{p}$ the infinum of the interval.

In order to obtain the angles associated to the slope at each dose we have the i-th angle is $\theta_i = arctan((y_{i+1}-y_i)/(x_{i+1}-x_i))$. Note that if the npS is calculated on the means at each dose then the squared error $\frac{1}{nm}\sum_{ij} (y_{ij}-f(x_i))^2$ will be minimized since the mean minimizes the L2 norm, alternatively if we calculate npS using the medians then the L1 norm would be minimzed. 

