[Back](./Index.md) | [Home](../Home.md)

---
#  Evolutionary model options
---

### Substitution models


    model={string}                                          A description of the substitution model to use, using the keyval syntax.

*The following nucleotide models are currently available:*

    JC69
    K80([kappa={real>0}])
    F84([kappa={real>0}, theta={real]0,1[}, theta1={real]0,1[},theta2={real]0,1[} ,"equilibrium frequencies"])
    HKY85([kappa={real>0}, theta={real]0,1[}, theta1={real]0,1[}, theta2={real]0,1[} ,"equilibrium frequencies"])
    T92([kappa={real>0}, theta={real]0,1[} ,"equilibrium frequencies"])
    TN93([kappa1={real>0}, kappa2={real>0}, theta={real]0,1[}, theta1={real]0,1[}, theta2={real]0,1[} ,"equilibrium frequencies"])
    GTR([a={real>0}, b={real>0}, c={real>0}, d={real>0}, e={real>0}, theta={real]0,1[}, theta1={real]0,1[}, theta2={real]0,1[} ,"equilibrium frequencies"])
    L95([beta={real>0}, gamma={real>0}, delta={real>0}, theta={real]0,1[}, theta1={real]0,1[}, theta2={real]0,1[} ,"equilibrium frequencies"])
    SSR([beta={real>0}, gamma={real>0}, delta={real>0}, theta={real]0,1[}])
    RN95([thetaR={real]0,1[}, thetaC={real]0,1[}, thetaG={real]0,1[}, kappaP={real[0,1[}, gammaP={real[0,1[}, sigmaP={real>1}, alphaP={real>1}])
    RN95s([thetaA={real]0,0.5[}, gamma={real]0,0.5[}, alphaP={real>1}])

*The following protein models are currently available:*

    JC69
    DSO78
    JTT92
    WAG01
    LG08
    LLG08_EX2([relrate1={real]0,1[}, relproba1={real]0,1[}])
    LLG08_EX3([relrate1={real]0,1[}, relrate2={real]0,1[}, relproba1={real]0,1[}, relproba2={real]0,1[}])
    LLG08_EHO([relrate1={real]0,1[}, relrate2={real]0,1[}, relproba1={real]0,1[}, relproba2={real]0,1[}])
    LLG08_UL2([relrate1={real]0,1[}, relproba1={real]0,1[}])
    LLG08_UL3([relrate1={real]0,1[}, relrate2={real]0,1[}, relproba1={real]0,1[}, relproba2={real]0,1[}])
    LGL08_CAT(nbCat={[10,20,30,40,50,60]}, [relrate1={real]0,1[}, relrate2={real]0,1[}, ..., relproba1={real]0,1[}, relproba2={real]0,1[}, ...] ))
    LGL08_CAT_C{[1,...,nbCat]}(nbCat={[10,20,30,40,50,60]})
    DSO78+F([theta={real]0,1[}, theta1={real]0,1[}, theta2={real]0,1[}, ... ,"equilibrium frequencies"])
    JTT92+F([theta={real]0,1[}, theta1={real]0,1[}, theta2={real]0,1[}, ..., "equilibrium frequencies"])
    WAG01+F([theta={real]0,1[}, theta1={real]0,1[}, theta2={real]0,1[}, ..., "equilibrium frequencies"])
    LG08+F([theta={real]0,1[}, theta1={real]0,1[}, theta2={real]0,1[}, ..., "equilibrium frequencies"])
    Empirical(name={chars}, file={path})
    Empirical+F(name={chars}, file={path}, [theta={real]0,1[}, theta1={real]0,1[}, theta2={real]0,1[}, ...,  "equilibrium frequencies"])

*The following meta models are currently available:*

    PIP(model={model description}, {lambda={real>0}(default=0.1), mu={real>0}(default=0.2) | estimated={true|false}})
    TS98(model={model description}, s1={real>0}, s2={real>0} [, "equilibrium frequencies"])
    G01(model={model description}, rdist={rate distribution description}, mu={real>0} [, "equilibrium frequencies"])
    RE08(model={model description}, lambda={real>0}, mu={real>0} [, "equilibrium frequencies"])

### Frequencies distribution sets
*The following frequencies distributions are available:*

    Fixed()                                                 All frequencies are fixed to their initial value and are not estimated.

    GC(theta={real]0,1[})                                   For nucleotides only, set the G content equal to the C content.

    Full(theta1={real]0,1[}, theta2={real]0,1[}, ..., thetaN={real]0,1[})
                                                            Full parametrization. Contains N free parameters.


### Rate across site distribution

    rate_distribution={rate distribution description}       Specify the rate across sites distribution

*The following distributions are currently available:*

**[Most used distributions]**

    Constant                                                Uses a constant rate across sites

    Gamma(n={int>=2}, alpha={float>0})                      A discretized gamma distribution of rates, with n classes, and a given shape,
                                                            with mean 1 (scale=shape).

    Invariant(dist={rate distribution description}, p={real[0,1]})
                                                            A composite distribution allowing a special class of invariant site, with a
                                                            probability p.

**[Standard distributions]**

    Beta(n={int>=2}, alpha={float>0}, beta={float>0})       A discretized beta distribution, with n classes, with standard parameters
                                                            alpha and beta.

    Gamma(n={int>=2}, alpha={float>0}, beta={float>0})      A discretized gamma distribution, with n classes, a shape alpha and a rate
                                                            beta, with parameters alpha and beta.

    Gaussian(n={int>=1}, mu={float}, sigma={float>0})
                                                            A discretized gaussian distribution, with n classes, a mean mu and a
                                                            standard deviation sigma, with parameters mu and sigma.

    Exponential(n={int>=2}, lambda={float>0})               A discretized exponential distribution, with n classes and parameter lambda.

    Simple(values={vector<double>}, probas={vector<double>} [, ranges={vector<parametername[min;max]>}])
                                                             A discrete distribution with specific values (in values) and their respective
                                                             non-negative probabibilities (in probas). The parameters are V1, V2, ..., Vn
                                                             for all the values and the relative probabibility parameters are
                                                             theta1, theta2, ..., thetan-1.
                                                             Optional argument {ranges} sets the allowed ranges of values taken by the
                                                             parameters; usage is like 'ranges=(V1[0.2;0.9],V2[1.1;999])'.

    TruncExponential(n={int>=2}, lambda={float>0}, tp={float>0})
                                                             A discretized truncated exponential distribution, with n classes, parameter
                                                             lambda and a truncation point tp. The parameters are lambda and tp.

    Uniform(n={int>=1}, begin={float>0}, end={float>0})      A uniform distribution, with n classes in interval [begin,end].
                                                             There are no parameters.

**[Mixture Distributions]**

    Invariant(dist={distribution description}, p={float>0})  A Mixture of a given discrete distributution and a 0 Dirac. p is the
                                                             probability of this 0 Dirac.

    Mixture(probas={vector<double>}, dist1={distribution description}, ..., distn={distribution description})
                                                             A Mixture of discrete distributions with specific probabilities (in probas)
                                                             and their respective desccriptions (in probas). The parameters are the relative
                                                             probabibility parameters theta1, theta2, ..., thetan-1, and the parameters
                                                             of the included distributions prefixed by Mixture.i_ where i is the order
                                                             of the distribution.
