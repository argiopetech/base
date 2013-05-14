Bayesian Analysis for Stellar Evolution
===============================================

Introduction
------------

Bayesian Analysis for Stellar Evolution with Nine Variables (BASE-9) is a Bayesian software suite that recovers star cluster and stellar parameters from photometry. BASE-9 is useful for analyzing single-age, single-metallicity star clusters, binaries, or single stars, and for simulating such systems. BASE-9 uses a Markov chain Monte Carlo (MCMC) technique along with bruteforce numerical integration to estimate the posterior probability distribution for the age, metallicity, helium abundance, distance modulus, line-of-sight absorption, and parameters of the initial-final mass relation (IFMR) for a cluster, and for the primary mass, secondary mass (if a binary), and cluster probability for every potential cluster member. The MCMC technique is used for the cluster quantities (first 6 items in the previous list) and numerical integration is used for the stellar quantities (last 3 items in the previous list). BASE-9 is freely available source code that you may use as is or modify for your own research and educational purposes.

BASE-9 may be the code for you if

1. you are dissatisfied with deriving cluster-level parameters by over-plotting isochrones on your data and iteratively adjusting parameters,
2. you wish to recover more than just an average and error bar for each parameter, and instead wish to characterize the probability distributions for these parameters, 
3. you wish to take fuller advantage of ancillary data, such as proper motion membership probabilities, spectroscopic mass estimates, or distancesfrom trigonometric parallaxes

Installation
------------

We designed BASE-9 to run under a variety of UNIX and Linux operating systems, though we have not tested it under a wide variety of systems. Currently we have BASE-9 running under Mac OS 10.6 and 10.8 and Linux Ubuntu 10.04. The code is written in the c programming language, though we expect to shortly add some c++ code because the latter handles memory more easily than does c. BASE-9 also takes advantage of parallel computing for the numerical integrations. You will need gcc (the c language compiler), gsl (the gnu science library), cmake (a cross-platform build system), and Open MPI (a high performance message passing library) to compile the code. To install these software packages, you may need help from your system administrator, though we provide some guidance in [the manual] (http://webfac.db.erau.edu/~vonhippt/base9/Manual_files/BASE-9_Manual.pdf).

References
----------
* [Embry-Riddle Aeronautical University] (http://www.erau.edu)
* [The BASE homepage] (http://webfac.db.erau.edu/~vonhippt/base9/)


