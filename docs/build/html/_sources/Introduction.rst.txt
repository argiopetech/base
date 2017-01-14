###############
I. Introduction
###############

Bayesian Analysis for Stellar Evolution with Nine Parameters (BASE-9) is a Bayesian software suite that recovers star cluster and stellar parameters from photometry. BASE-9 is useful for analyzing single-age, single-metallicity star clusters, binaries, or single stars, and for simulating such systems. This document assumes you are working with base 9.4.3. We will endeavor to update this manual as we update the code or as libraries or operating systems meaningfully change.

BASE-9 uses a Markov chain Monte Carlo (MCMC) technique along with brute-force numerical integration to estimate the posterior probability distribution for up to six cluster and three stellar properties. The cluster properties are age, metallicity, helium abundance, distance modulus, line- of-sight absorption, and parameters of the initial-final mass relation (IFMR). The stellar properties are primary mass, secondary mass (if a binary), and cluster membership probability. The MCMC technique is used for the cluster quantities and numerical integration is used for the stellar quantities. BASE-9 is freely available source code that you may use as is or modify for your own research and educational purposes.

BASE-9 may be the code for you if

1. you are dissatisfied with deriving cluster-level parameters by over-plotting isochrones on your data and iteratively adjusting parameters,
2. you wish to recover more than just an average and error bar for each parameter, and instead wish to characterize the probability distributions for these parameters,
3. you wish to take fuller advantage of ancillary data, such as proper motion membership probabilities, spectroscopic mass estimates, or distances from trigonometric parallaxes.

This manual is designed to help you install and run BASE-9. If you use BASE-9 in your research, please cite

``von Hippel, T., Jefferys, W. H., Scott, J., Stein, N., Winget, D. E., DeGennaro, S., Dam, A., & Jeffery, E. 2006, Inverting Color-Magnitude Diagrams to Access Precise Star Cluster Parameters: A Bayesian Approach, ApJ, 645, 1436``

and if you find the following helpful, please also cite

``DeGennaro, S., von Hippel, T., Jefferys, W. H., Stein, N., van Dyk, D. A., & Jeffery, E. 2009, Inverting Color-Magnitude Diagrams to Access Precise Star Cluster Parameters: A New White Dwarf Age for the Hyades, ApJ, 696, 12``

``van Dyk, D. A., DeGennaro, S., Stein, N., Jefferys, W. H., & von Hippel, T. 2009, Statistical Analysis of Stellar Evolution, Annals of Applied Statistics, 3, 117``

Depending on how you use BASE-9 (this part is under your control), the software also relies on the stellar evolution models of

``Dotter, A., Chaboyer, B., Jevremovic, D., Kostov, V., Baron, E., & Ferguson, J. W. 2008, The Dartmouth Stellar Evolution Database, ApJS, 178, 89``

``Girardi, L., Bressan, A., Bertelli, G., & Chiosi, C. 2000, Evolutionary tracks and isochrones for low- and intermediate-mass stars: From 0.15 to 7 Msun, and from Z=0.0004 to 0.03, A&AS, 141, 371``

``Yi, S., Demarque, P., Kim, Y.-C., Lee, Y.-W., Ree, C. H., Lejeune, T., & Barnes, S. 2001, Toward Better Age Estimates for Stellar Populations: The Y2 Isochrones for Solar Mixture, ApJS, 136, 417``

the white dwarf atmosphere models of

``Bergeron, P., Wesemael, F., & Beauchamp, A. 1995, Photometric Calibration of Hydrogen- and Helium-Rich White Dwarf Models, PASP, 107, 1047
(as updated and made available at http://www.astro.umontreal.ca/~bergeron/CoolingModels/)``

the white dwarf interior models of

``Althaus, L. G. & Benvenuto, O. G. 1998, Evolution of DA white dwarfs in the context of a
new theory of convection, MNRAS, 296, 206``

``Montgomery, M. H., Klumpe, E. W., Winget, D. E., & Wood, M. A. 1999, Evolutionary Calculations of Phase Separation in Crystallizing White Dwarf Stars, ApJ, 525, 482
(The paper describes the stellar evolution code that M. Montgomery used in 2012 to calculate the WD sequences specifically for use with BASE-9.)``

``Renedo, I., Althaus, L. G., Miller Bertolami, M. M., Romero, A. D., Corsico, A. H., Rohrmann, R. D., & Garcia-Berro, E. 2010, New Cooling Sequences for Old White Dwarfs, ApJ, 717, 183``

``Wood, M. A. 1992, Constraints on the age and evolution of the Galaxy from the white dwarf luminosity function, ApJ, 386, 539``

the Initial Mass Function of

``Miller, G. E., & Scalo, J. M. 1979, The initial mass function and stellar birthrate in the solar
neighborhood, ApJS, 41, 513``

and the Initial-Final Mass Relations of

``Salaris, Salaris, Maurizio; Serenelli, Aldo; Weiss, Achim; Miller Bertolami, Marcelo. 2009,
Semi-empirical White Dwarf Initial-Final Mass Relationships: A Thorough Analysis of Systematic Uncertainties Due to Stellar Evolution Models, ApJ, 692, 1013``

``Weidemann, V. 2000, Revision of the initial-to-final mass relation, A&A, 363, 647``

``Williams, K. A., Bolte, M., & Koester, D. 2009, Probing the Lower Mass Limit for Supernova Progenitors and the High-Mass End of the Initial-Final Mass Relation from White Dwarfs in the Open Cluster M35 (NGC 2168), ApJ, 693, 355``

or a fitted IFMR parameterized as lines, broken lines, or low-order polynomials as described by

``Stein, N. M., van Dyk, D. A., von Hippel, T., DeGennaro, S., Jeffery, E. J., & Jefferys, W. H. 2013, Combining Computer Models in a Principled Bayesian Analysis: From Normal Stars to White Dwarf Cinders, Statistical Analysis and Data Mining, 6, 34``

For a further discussion of what BASE-9 and its precursor, BASE-8, has been used for to date and some indications of how it might be useful in your research, see also the following papers:

``Jeffery, E. J., von Hippel, T., Jefferys, W. H., Winget, D. E., Stein, N., & DeGennaro, S., 2007, New Techniques to Determine Ages of Open Clusters Using White Dwarfs, ApJ, 658, 391``

``Jeffery, E. J., von Hippel, T., DeGennaro, S., Stein, N, Jefferys, W.H., & van Dyk, D. 2011, The White Dwarf Age of NGC 2477, ApJ, 730, 35``