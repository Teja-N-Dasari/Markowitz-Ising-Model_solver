# Markowitz-Ising-Model_solver

We convert the Markowitz Problem, which is a minimization problem by solving which we can obtain a portifolio with minimized risk and maximized returns, into a Ising type Hamiltoninan

We are supposed to input the expected returns $r_i$, Correlations or Covarince matrix $C_ij$ and the budget $b$ into this model

Here, it is the Ising type model but here all spins interact with each other not just nearest neighbours.

Then the ck-annealing method is used to solve that Ising problem, by introducing the magnetization $m_i$ = $\langle s_i \rangle$, which is continous $m_i$ [+1,1].

