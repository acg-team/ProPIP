[Back](../Index.md) | [Home](https://github.com/acg-team/ProPIP/blob/master/ProPIP.wiki/ProPIP-Progressive-Multiple-Sequence-Alignment-with-Poisson-Indel-Process.md)

---
#  Numerical parameter optimisation options
---


This program allows to (re-)estimate numerical parameters, including
- Branch lengths
- Entries of the substitution matrices, included base frequencies values)
- Parameters of the rate distribution (currently shape parameter of the gamma law, proportion of invariant sites).


    optimization={method}


*The following methods are currently available:*

    None                                                        No optimization is performed, initial values are kept 'as is'.

    FullD(derivatives={Newton|Gradient})                         Full-derivatives method.
                                                                 Branch length derivatives are computed analytically, others numerically.
                                                                 The derivatives arguments specifies if first or second order derivatives
                                                                 should be used. In the first case, the optimization method used is the
                                                                 so-called conjugate gradient method, otherwise the Newton-Raphson method
                                                                 will be used.

    D-Brent(derivatives={Newton|Gradient|BFGS}, nstep={int>0})   Branch lengths parameters are optimized using either the conjugate gradient,
                                                                 the Newton-Raphson method or BFGS, other parameters are estimated using the
                                                                 Brent method in one dimension. The algorithm then loops over all parameters
                                                                 until convergence. The nstep arguments allow to specify a number of
                                                                 progressive steps to perform during optimization. If nstep=3 and
                                                                 precision=E-6, a first optimization with precision=E-2, will be performed,
                                                                 then a round with precision set to E-4 and finally precision will be set to
                                                                 E-6. This approach generally increases convergence time.

    D-BFGS(derivatives={Newton|Gradient|BFGS}, nstep={int>0})    Branch lengths parameters are optimized using either the conjugate gradient,
                                                                 the Newton-Raphson method or BFGS, other parameters are estimated using the
                                                                 BFGS method in one dimension. The algorithm then loops over all parameters
                                                                 until convergence. The nstep arguments allow to specify a number of
                                                                 progressive steps to perform during optimization. If nstep=3 and
                                                                 precision=E-6, a first optimization with precision=E-2, will be performed,
                                                                 then a round with precision set to E-4 and finally precision will be set to
                                                                 E-6. This approach generally increases convergence time.

    ND-Brent(derivatives={Newton|Gradient|BFGS|Brent}, nstep={int>0})
                                                                  Branch lengths parameters are optimized using either the conjugate gradient,
                                                                  the Newton-Raphson method, BFGS, or Brent (monodimesional) other parameters
                                                                  are estimated using the Brent method in one dimension. The algorithm then
                                                                  loops over all parameters until convergence. The nstep arguments allow to
                                                                  specify a number of progressive steps to perform during optimization. If
                                                                  nstep=3 and precision=E-6, a first optimization with precision=E-2, will
                                                                  be performed, then a round with precision set to E-4 and finally precision
                                                                  will be set to E-6. This approach generally increases convergence time.

    ND-BFGS(derivatives={Newton|Gradient|BFGS|Brent}, nstep={int>0})
                                                                  Branch lengths parameters are optimized using either the conjugate gradient,
                                                                  the Newton-Raphson method, BFGS, or Brent (monodimesional) other parameters
                                                                  are estimated using the BFGS method in one dimension. The algorithm then
                                                                  loops over all parameters until convergence. The nstep arguments allow to
                                                                  specify a number of progressive steps to perform during optimization. If
                                                                  nstep=3 and precision=E-6, a first optimization with precision=E-2, will
                                                                  be performed, then a round with precision set to E-4 and finally precision
                                                                  will be set to E-6. This approach generally increases convergence time.


    optimization.reparametrization=<bool>                         Tells if parameters should be transformed in order to remove constraints (for
                                                                  instance positivie-only parameters will be log transformed in order to obtain
                                                                  parameters defined from -inf to +inf). This may improve the optimization,
                                                                  particularly for parameter-rich models, but the likelihood calculations
                                                                  will take a bit more time.
    optimization.final={powell|simplex|bfgs}                      Optional final optimization step, useful if numerical derivatives are
                                                                  to be used. Leave the field empty in order to skip this step.
    optimization.profiler={{path}|std|none}                       A file where to dump optimization steps (a file path or std for standard
                                                                  output or none for no output).
    optimization.message_handler={{path}|std|none}                A file where to dump warning messages.
    optimization.max_number_f_eval=<int>0>                        The maximum number of likelihood evaluations to perform.
    optimization.tolerance=<float>0>                              The precision on the log-likelihood to reach.
    optimization.constrain_parameter={list<chars=interval>}       A list of parameters on which the authorized values are limited to a given
                                                                  interval. For example:
                                                                  optimization.constrain_parameter = YN98.omega = [-inf;1.9[, *theta* = [0.1;0.7[, BrLen*=[0.01;inf]
    optimization.ignore_parameter={list<chars>}                   A list of parameters to ignore during the estimation process. The parameter
                                                                  name should include there 'namespace', that is their model name, for
                                                                  instance K80.kappa,TN93.theta, GTR.a, Gamma.alpha, etc.'BrLen' will
                                                                  ignore all branch lengths and 'Model' will ignore all model parameters.
                                                                  The '*' wildcard can be used, as in *theta* for all the parameters whose
                                                                  name has theta in it.
