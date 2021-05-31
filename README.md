# Gas Field Multi-period Optimization Under Harmonic Decline

## Abstract
In operating a gas field, petroleum firms would like to maximize the condensate production. This process is usually done with a single-period linear optimization. However, some production constraints, especially the Harmonic decline pattern of the wells, are nonlinear and nonconvex. Moreover, the solution to a single-period optimization problem can be drastically different from the multi-period one as the former ignores the time domain completely. 

In this project, we formulate a multi-period production optimization problem with Harmonic decline as a nonlinear optimization problem. We then solve it with the augmented Lagrangian technique and use the gradient projection method with trust-region to solve the subproblem. We show that the technique works perfectly under different scenarios by reporting the corresponding production schemes and explaining their semantics.

**Keywords**: Harmonic decline; Gas field optimization; Multi-period; Augmented Lagrangian method

## Objectives

Solve Gas field multi-period optimization problem under harmonic decline with the nonlinear augmented Lagrangian method.

## Dataset

Simulated data

## Software

MATLAB with a basic linear algebra package

## Key file description
* Paper_Chayut_Wongkamthong.pdf : Final report of the project
* AL_main.m : Main MATLAB script to run simulations for different scenarios.
* computeKKT_AL.m : Function to evaluate the optimality condition.
* gen_case_1.m : Function to generate an instance of case I problem.
* gen_case_2.m : Function to generate an instance of case II problem.
* gen_case_3.m : Function to generate an instance of case III problem.
* gen_case_4.m : Function to generate an instance of case IV problem.
* gradient_checking.m : MATLAB script to perform gradient checking for different constraints and objective function using the finite differences technique.
* gradient_checking_AL.m : MATLAB script to perform gradient checking for the augmented Lagrangian function using the finite differences technique.
* plot_q.m : Function to generate the plot of production contribution over time.
* plot_t.m : Function to generate the plot of the open duration over time.
* project.m : Function to perform projection onto box constraints $P(x,l,u)$.
* solution_case_1.mat : A solution of an instance of case I problem.
* solution_case_2.mat : A solution of an instance of case II problem.
* solution_case_3.mat : A solution of an instance of case III problem.
* solution_case_4_long.mat : A solution of an instance of case IV problem with long field life ($T = 15$).
* solution_case_4_short.mat : A solution of an instance of case IV problem with short field life ($T = 3$).
* gradHess\ALagrangian.m : Function to evaluate the augmented Lagrangian function together with its gradient and Hessian matrices.
* gradHess\combineConst.m : Function to get and combine constraints and their gradients and Hessian matrices.
* gradHess\computeHarmonicConstr.m : Function to evaluate the Harmonic decline constraints together with its gradient and Hessian matrices.
* gradHess\computeNomConstr.m : Function to evaluate the nomination constraints together with its gradient and Hessian matrices.
* gradHess\computeObjGradHess.m : Function to evaluate the objective function together with its gradient and Hessian matrices.
* gradHess\computeTimeConstr.m : Function to evaluate the time constraints together with its gradient and Hessian matrices.
* TR_ALsubproblem\calculate_t_bound.m : Function to calculate $t$ until we can move to hit each of the boundaries in the negative gradient direction in the gradient projection technique in order to find the Cauchy point.
* TR_ALsubproblem\CG_subproblem.m : Function to perform subspace minimization to find a better solution than the Cauchy point. This is used in the second step of the gradient projection technique to solve the AL subproblem.
* TR_ALsubproblem\getActiveSet.m : Function to get the active set indices from the current iterate/Cauchy point.
* TR_ALsubproblem\getCauchypoint.m : Function to find the Cauchy point in the gradient projection technique.
* TR_ALsubproblem\solveWithTR.m : Function to solve the AL subproblem using the gradient projection with trust-region.

**Author**

Chayut Wongkamthong - Social Science Research Institute (SSRI), Duke University
