Bayesian Analysis for Stellar Evolution
===============================================

**Note**: Running BASE-9 requires [the models](https://github.com/argiopetech/base-models).


Introduction
------------

Bayesian Analysis for Stellar Evolution with Nine Variables (BASE-9) is a Bayesian software suite that recovers star cluster and stellar parameters from photometry. BASE-9 is useful for analyzing single-age, single-metallicity star clusters, binaries, or single stars, and for simulating such systems. BASE-9 uses a Markov chain Monte Carlo (MCMC) technique along with bruteforce numerical integration to estimate the posterior probability distribution for the age, metallicity, helium abundance, distance modulus, line-of-sight absorption, and parameters of the initial-final mass relation (IFMR) for a cluster, and for the primary mass, secondary mass (if a binary), and cluster probability for every potential cluster member. The MCMC technique is used for the cluster quantities (first 6 items in the previous list) and numerical integration is used for the stellar quantities (last 3 items in the previous list). BASE-9 is freely available source code that you may use as is or modify for your own research and educational purposes.

BASE-9 may be the code for you if

1. you are dissatisfied with deriving cluster-level parameters by over-plotting isochrones on your data and iteratively adjusting parameters,
2. you wish to recover more than just an average and error bar for each parameter, and instead wish to characterize the probability distributions for these parameters, 
3. you wish to take fuller advantage of ancillary data, such as proper motion membership probabilities, spectroscopic mass estimates, or distancesfrom trigonometric parallaxes

Installation
------------

We designed BASE-9 to run under a variety of UNIX and Linux operating systems, though we have not tested it under a wide variety of systems. Currently we have BASE-9 running under Mac OS 10.6 and 10.8 and Linux Ubuntu 10.04. The code is written primarily in the C++ programming language. BASE-9 also takes advantage of parallel computing for the numerical integrations using a hybrid scheme utilizing MPI and C++11 threads. You will need g++ (a C++ compiler associated with GCC) version 4.8 or newer or clang++ version 3.2 or newer (available with XCode 4.6), GSL (the Gnu Science Library), cmake (a cross-platform build system) version 2.8.10 or newer, Boost (a peer-reviewed C++ library) 1.54 or newer, and Open MPI (a high performance message passing library) to compile the code. To install these software packages, you may need help from your system administrator, though we provide some guidance in [the manual] (http://webfac.db.erau.edu/~vonhippt/base9/Manual_files/BASE-9_Manual.pdf).

C++11
-----
As of version 9.3.0, BASE-9 is dependent upon a C++11 compliant compiler for proper operation. Compilation is currently supported with gcc/g++ 4.8, Clang 3.3 (preferred) and 3.2 (for XCode 4.6 compatibility). An incomplete list of utilized C++11 features includes:

* C++11 threads
* Thread-safe static initialization (per [N2660](http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2008/n2660.htm))
* Template aliases
* std::array
* C++11 \<random\> (Mersenne Twister and distributions)

References
----------
* [Embry-Riddle Aeronautical University] (http://www.erau.edu)
* [The BASE homepage] (http://webfac.db.erau.edu/~vonhippt/base9/)


