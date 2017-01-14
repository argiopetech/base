#############################
V. Diagnostics of run quality
#############################

The following two plots show examples of poor and good sampling. In the first (extreme) case, the age sampling is highly correlated and one would need to use post-run thinning, probably by a factor of ~100. This means that one would need to run the code for 100x as many interations. The metallicity sampling displays only minor correlation, and if all other parameters looked this uncorrelated then this run would be sufficient. In this particular case, both plots were generated from the same singlePopMcmc run and because no single parameter is reliable until all parameters are essentially uncorrelated, this run did not reliably determine the metallicity (or any other) posterior distribution.

.. image:: charts_1.png
    :width: 1015px
    :align: center
    :height: 831px