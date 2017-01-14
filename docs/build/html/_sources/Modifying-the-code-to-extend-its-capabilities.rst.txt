##################################################
VII. Modifying the code to extend its capabilities
##################################################

We continue to upgrade BASE-9 for our on-going projects. If you wish to add capability to BASE-9, we will be happy to suggest to you how best to go about this and try to estimate the work involved. Here is an example list of how involved a variety of tasks are likely to be.

*Less than 2 hours:* Modifying the IFMR. You can do this by editing or adding a few lines of code in ifmr.cpp.

*Less than 8 hours:* Change the IMF. You will need to create a subroutine where a random mass value can be drawn from your IMF distribution. This currently takes place in drawFromIMF.cpp. Note that you will also have to normalize the IMF for the Bayesian routine to work properly and that this takes place in densities.cpp and is stored in logMassNorm.

*Less than 16 hours:* Incorporating another set of stellar evolution models – see instructions at the top of msRgbEvol.cpp and possibly wdCooling.cpp and/or gBergMag.cpp.

*Less than a week:* Sampling a new variable (e.g. stellar rotation, alpha-element enhancement). This takes place primarily in singlePopMcmc/MpiMcmcApplication.cpp and base9/densities.cpp.