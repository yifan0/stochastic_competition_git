# stochastic_competition_git
sim_git.jl uses DynamicGrids package to simulate one numeric value, Ny, and its variation across a 2D grid over time. <br />
The processes affecting the values of $N_y$ in a cell are <br />
* local speciation events, which can increase or decrease $N_y$, with probabilities of success determined by the size and direction of mutation. <br />
* invasions by neighboring cells, where the value of $N_y$ is (if the invasion is successful) replaced by its neighbor’s value of $N_y$. <br />


The outputs we are interested in relate to the distribution of values of $N_y$ as a function of box size, system size, and speciation rate. <br />
We are computing
* the width function (width_git.jl), which is defined in terms of the variance in $\log(N_y)$ for a set of cells, and in other contexts is called roughness.
* the diversity function (width_git.jl), which is the number of unique values of $N_y$ in a sample, corresponding to the number of distinct species we have in that sample.
* the range function (range_git.jl)

In sim_git.jl, for each of the 100 timesteps and each of the 10 simulations, the code simulates $N_y$ with speciation and invasion rules.

In width_git.jl/range_git.jl, for each sample size $i$, the code
1. takes $(size-i)$ number of $i\times i$ samples in each of the 10 independent simulations and measures the width/diversity/range of each of them.
2. averages the width/diversity/range over the $(size-i)\quad i\times i$ samples.
3. averages the width/diversity/range over the 10 independent simulations