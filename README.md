# chanceDPA

## About this project
This program solves a dynamic programming problem with joint probabilistic constraints. It can be used to solve a path planning problem where the chance of hitting an obstacle needs to be kept below a threshold.

This project implements this paper [1] with a few modifications.  <br />

## Results
### Path planning problem definition 
The code was tested on a path planning problem with the following dynamics (see [1]):<br />
![x_{k+1}=f_k(xk,uk,wk)=xk+uk+wk](https://latex.codecogs.com/svg.latex?x_%7Bk&plus;1%7D%20%3D%20f_k%28x_k%2Cu_k%2Cw_k%29%20%3D%20x_k%20&plus;%20u_k%20&plus;%20w_k)<br />
![norm(u_k), wk_distrib](https://latex.codecogs.com/gif.latex?%5Cleft%20%5C%7C%20u_k%5Cright%20%5C%7C_2%20%5Cleq%20d_k%2C%20w_k%5Csim%20N%280%2C%5Csigma%5E2I%29)<br />
![norm(u_k), wk_distrib](https://latex.codecogs.com/gif.latex?%5Cdpi%7B120%7D%20%5Cleft%20%5C%7C%20u_k%5Cright%20%5C%7C_2%20%5Cleq%20d_k%2C%20w_k%5Csim%20N%280%2C%5Csigma%5E2I%29)<br />
Note that ![xk](https://latex.codecogs.com/gif.latex?x_k) consists of 100x100 states, ![uk](https://latex.codecogs.com/gif.latex?u_k) can take up to 81 different values using our discretized grid. <br />
As in [1], we use ![sigma167](https://latex.codecogs.com/gif.latex?%5Csigma%3D1.67). Therefore, to minimize the complexity of the problem, we approximate ![wk](https://latex.codecogs.com/gif.latex?w_k) to discretized values, with <br />
![abs(wk)<5](https://latex.codecogs.com/gif.latex?%5Cleft%20%7C%20w_k%20%5Cright%20%7C%3C5).  <br />
To compute efficiently the expected value of the cost for each state and action, which is equal to:<br />
![Jk_fixeduk](https://latex.codecogs.com/gif.latex?J_k%5E%5Clambda%28x_k%29_%7B%5Cmu_k%3Du_k%7D%3DE_%7Bw_k%7D%5Cleft%20%5C%7B%20L_k%5E%5Clambda%28x_k%2Cu_k%29&plus;J_%7Bk&plus;1%7D%28f%28x_k%2Cu_k%2Cw_k%29%29%20%5Cright%20%5C%7D),<br />
we then convolve this with a probability filter. (APPROACH TO BE CLARIFIED HERE)<br />
This which avoids the need to loop through different ![wk](https://latex.codecogs.com/gif.latex?w_k) values to obtain an estimate of the cost for each action. This approach (tested) is computationnally too expensive. <br />



### Plots
Typical costs results after N=50 steps, lambda = 0 (no chance constraint). <br /> <br />
![alt text](https://github.com/thomasjlew/chanceDPA/blob/master/old_files_201117/cost.png)

### Execution time
...

## Further work
- The ongoing work aims at improving the performance of the dynamic programming algorithm to allow for more iterations.

## References
- [1] M. Ono, M. Pavone, Y. Kuwata and J. Balaram, "Chance-constrained dynamic programming with application to risk-aware robotic space exploration", Autonomous Robots, 2015.
- [2] D.P. Bertsekas, "Dynamic Programming and Optimal Control", Vol. I, 3rd edition, 2005, 558 pages.
- [3] D.P. Bertsekas, "Nonlinear Programming", Second Edition. Athena Scientific (1999)
