##################
IV. Running BASE-9
##################

In the following subsections, we describe how to run stand-alone portions of the BASE-9 modules from the command line. There are various reasons why you might want to run one or another of these, or some, but not all, so we detail how to run each one. As of BASE-9.2.0, all settings have been moved into a YAML-format configuration file. A sample configuration file with reasonable initial settings can be found in the ``base-9.4.2/conf`` directory under the name base9.yaml. A sample cshell script can be found at scripts/hyades.csh. Individual settings can be changed on a run-by-run basis via the command line options. Run any of the BASE-9 applications with the command line flag “--help” to view a description of available settings.

The following examples assume you have installed the BASE-9 executables in a directory that is in your PATH (e.g., ``/usr/local/bin``). If this is not the case, you may need to use absolute pathnames (e.g., ``/home/me/base-9.4.2/BUILD/bin/singlePopMcmc``).

*************
A. simCluster
*************

The first tool that you are likely to use within BASE-9 is simCluster. This module simulates a stellar cluster for a particular set of models (see references in the Introduction) and user-specified values of various cluster parameters that have been set in the base9.yaml file: the base-10 log of the cluster age, metallicity, the helium mass fraction (only for some Dotter et al. models), the distance modulus, absorption in the V-band, the percent of cluster stars that are binaries, the upper mass limit for creating a white dwarf (WD), and the fraction of WDs that have helium atmospheres (DBs). We recommend leaving these last two parameters at 8 solar masses and 0%, as we have not yet fully implemented and tested them. An additional parameter, the seed to the random number generator, is necessary as the mass of each star is determined by randomly drawing from the IMF. This allows you to specify multiple clusters with the same parameters but different random seeds if you wish to test the effects of, for instance, cluster size on the number of WDs or the clarity of the main sequence turn-off (MSTO). This seed can be set via the ``--seed`` option in the command line. To run simCluster, simply type its name:

::

	linux> ./simCluster
	Seed: 1559729633
	Reading models... Done.
	Properties for cluster:
	 logClusAge     =  8.796
	 [Fe/H]         =  0.07
	 Y              =  0.29
	 modulus        =  0.00
	 Av             =  0.01
	 WDMassUp       =  8.0
	 fractionBinary =  0.00
	Totals:
	 nSystems		=  100
	 nStars			=  100
	 nMSRG			=  98
	 nWD			=  2
	 nNSBH			=  0
	 massTotal 		=  62.72
	 MSRGMassTotal  =  55.10
	 wdMassTotal    =  1.66

The above output is diagnostic and reiterates the settings in the base9.yaml file. The stored output of simCluster is placed in a filename specified by the user with the outputFileBase option in the base9.yaml file. The file contents from the output of simCluster should look like the following:

::

	linux>head -2 hyades.sim.out
	  id      U      B      V      R      I      J      H      K   sigU
	sigB   sigV   sigR   sigI   sigJ   sigH   sigK   mass1 massRatio stage
	Cmprior useDBI
	   1  3.061  3.031  2.694  2.501  2.320  2.127  1.983  1.965  0.000
	0.000  0.000  0.000  0.000  0.000  0.000  0.000   1.555     0.000
	1   0.999      1

There are more columns than can be presented cleanly on a page, but hopefully this is clear enough. The first column lists identification numbers for each star system (single star or binary). This is meant to be useful in tracking down particular stars. The next eight columns list the U- through K-band magnitudes (or ugriz through K) of the primary star. Columns 10 through 17 give the photometric uncertainties for each filter entry (for a simulated cluster, these are zero).

The 18th and 19th columns give the mass of the primary star and the mass of its companion (if applicable). The 20th column lists the stage of stellar evolution for that particular star (1 = MS or RG, 3 = WD, >3 for evolved stars above the WD mass limit). The final two columns are the cluster membership prior (which is essentially ~1 for simulated stars) and the flag (0 or 1) whether to use the star during the burn-in stage. With these final columns, the output file is formatted for input into scatterCluster. You should be able to plot reasonable looking CMDs/isochrones from this file for a wide range of cluster parameters, stellar models, and filters.

*****************
B. scatterCluster
*****************

The scatterCluster module adds Gaussian random errors to the photometry output created by simCluster. To specify the appropriate amount of error to add for your particular simulation, adjust the virtual exposure time in the base9.yaml file.

We use exposure times of 1 hour in each filter to generate a scattered cluster with the above file. The algorithm for adding noise to the cluster photometry is rudimentary and only meant for simple purposes such as preparing for an observing proposal or for creating test files for the Markov chain Monte Carlos (singlePopMcmc) routine. The algorithm is an approximation to the results one would obtain in one hour with the KPNO 4m + Mosaic (UBVRI) or Flamingos (JHK), assuming dark time, seeing=1.1 arcsec, airmass=1.2. Signal-to-noise for the Spitzer bands, if included, are naively set to be the same as for the K-band. For departures from a one- hour exposure the S/N is scaled by sqrt(exptime). These exposure times can be set in the base9.yaml file under exposures for each individual filter.

Additional options for scatterCluster are available in the yaml file. The bright and faint end cut-off mags allow you to narrow the portion of the CMD that you wish to retain. The relevantFilt option specifies which band is the reference filter (in this case, 0=U, 1=B, etc.). The base9.yaml options brightLimit and faintLimit refer to the bright and faint end cut-off magnitudes for the reference filter indicated. You can also clip on S/N with limitS2N and decide to cut out field stars, if they were simulated by simCluster. Additionally, scatterCluster will determine which filters you are using based on the header in the simCluster output file. Again, the integer seed may be set at the command line to allow you to start from the same input file, but create multiple simulated observations of that file with different initial seed values.

::

	linux> scatterCluster
	Seed: 1564704505

The output file of scatterCluster looks like

::

	linux> head -2 hyades/hyades.scatter.out
	  id      U      B      V      R      I      J      H      K   sigU
	sigB   sigV   sigR   sigI   sigJ   sigH   sigK   mass1 massRatio stage
	Cmprior useDBI
	1 3.065 3.016 2.691 2.509 2.328 2.128 1.985 1.977 0.010 0.010 0.010 0.010 0.010 0.010 0.010 0.010 1.555 0.000 1 0.999 1

Notice now that the output includes the estimated errors for each band (sig-). The format of the output file is otherwise the same as the input file for scatterCluster.

In this case, only the id, mass1, and stage1 values are kept from the output of simCluster. The photometry values (here UBVRIJHK) are derived from the photometry values in the simCluster output file, but are different in that they are scattered by adding a Gaussian random deviate with sigma = sigU, sigB, etc. This section of the output file is all one needs to plot realistic CMDs for proposals and possibly to prepare for observing projects. The scatterCluster output file contains additional information, however, and is formatted to be ingested by singlePopMcmc, so that it can be used to test singlePopMcmc and so that you can test the precision and accuracy that you would expect to recover from real data based on a given set of cluster parameters, observational errors, and the number of stars available. The massRatio column lists the ratio by mass of the secondary to primary stars, which in these examples are both 0 since there were no secondaries. The CMprior column is set by default in scatterCluster to 0.99, but the file can easily be edited to set a different prior probability that any particular star is a cluster member. The final column is just a 0 or 1 switch (off or on) of whether to use a particular star during the burnin process. (See DeGennaro et al. 2009 and van Dyk et al. 2009 for a discussion of what the burnin entails and why it is used.) To make it easiest for singlePopMcmc to converge, it is helpful to have this parameter set to 1 for stars that are likely to be cluster members and if there are many field stars, it is helpful if the bulk of them can be set to 0 at this point.

****************
C. singlePopMcmc
****************

The singlePopMcmc module is the workhorse of our software suite. This routine, along with its many subroutines, runs a Markov chain Monte Carlo sampler using a variety of standard Bayesian techniques as well as a few techniques newly developed by us. The approach and mathematics are presented by DeGennaro et al. (2009), van Dyk et al. (2009), and Stein et al. (2013). This code was designed to run on photometry formatted in the same manner as the output of scatterCluster. It can also be run just as easily on the simulated photometry from simCluster + scatterCluster.

The singlePopMcmc module has a variety of values and options set in the base9.yaml file. Under the singlePopMcmc group, the stage2IterMax and stage3Iter set the length of the burnin for singlePopMcmc. The runIter option lets you choose the number of iterations of the Markov chain Monte Carlo. The rule-of-thumb is that one typically wants 10,000 well- sampled points from a Markov chain Monte Carlo in order to draw robust inferences on the posterior distribution. At the other extreme, the Central Limit Theorem dictates that approximately 30 uncorrelated samples are sufficient for a normal distribution. Before running a particular dataset against a specific set of models, you do not know if the posterior distributions will be Gaussian shaped or more complex, so we suggest you take the conservative approach and initially assume complex posterior distributions and run BASE-9 for 10,000 uncorrelated iterations. The parameter thin sets the increment between saved iterations. We recommend that this parameter be left equal to 1 to keep the adaptive sampling routine efficient. If the output of singlePopMcmc is correlated (see below), then each new iteration or step is not independent and you need substantially more than 10,000 iterations to draw robust inferences. In situations like this, we recommend that the output file be thinned afterwards, i.e. that the user uses every nth record where n is large enough to keep the output uncorrelated.

Under the cluster options, there are five parameters for which means and standard deviations can be set: the metallicity prior (Fe_H), the distance modulus prior (distMod), the absorption prior (Av), the helium prior (Y), and the carbon fraction prior for a C+O WD (carbonicity). Note that carbonicity only works with the Montgomery models and is not yet supported because we are currently testing it. If you only have weak priors, that is fine. If you do not want to sample on one or more of these parameters, you can set the sigma for that parameter to 0.0 and this will turn of sampling for that parameter. Under starting, the parameter logClusAge is a starting value for the log of the age in years (e.g. 9.0 for a 1 billion year old cluster). *This is not a prior*, but just tells singlePopMcmc where to start searching for a fit. We have found that although convergence may depend on starting with a roughly reasonable age, the actual posterior age distribution does not depend on what that value is, assuming it does converge.

The msRgbModel lets you choose which set of models to use with your data (the filters available in the models must match the filters of your observed or simulated/scattered cluster). This allows you to derive cluster parameters for a range of models as well as to create simulated clusters under one set of models and use singlePopMcmc to derive the cluster and stellar parameters under another set of models. The latter experiments might be useful, for instance, if you wanted to test the sensitivity of basic cluster or stellar parameters to a given model ingredient. With ancillary data for cluster or stellar parameters this might allow you to constrain model ingredients.

Again we mention that the seed can be set inline with ``--seed`` when singlePopMcmc is called. If singlePopMcmc appears to be unable to converge on reasonable cluster values, rerun it with a different initial seed. Changing the seed also allows you to start a new MCMC chain if you ran a prior calculation with too few iterations.

To run singlePopMcmc, using a properly prepared input base9.yaml file, type the following:

::

	linux> singlePopMcmc --verbose
	Bayesian Analysis of Stellar Evolution
	Seed: 1570065938
	Reading models... Done.
	Model boundaries are (7.800, 10.250) log years.
	Binaries are OFF
	Running Stage 1 burnin... Complete (acceptanceRatio = 0.090)
	Running Stage 2 (adaptive) burnin...
	    Acceptance ratio: 0.350. Trying for trend.
	    Acceptance ratio: 0.600. Retrying.
	    Acceptance ratio: 0.380. Trying for trend.
	    Acceptance ratio: 0.520. Retrying.
	    Acceptance ratio: 0.280. Trying for trend.
	    Acceptance ratio: 0.180. Retrying.
	    Acceptance ratio: 0.400. Trying for trend.
	    Acceptance ratio: 0.440. Retrying.
	    Acceptance ratio: 0.320. Trying for trend.
	    Acceptance ratio: 0.500. Retrying.
	    Acceptance ratio: 0.240. Trying for trend.
	  Leaving adaptive burnin early with an acceptance ratio of 0.220
	(iteration 1300)
	Starting adaptive run... Preliminary acceptanceRatio = 0.300

The singlePopMcmc routine creates multiple output files. In this case, it created:

::

	-rw-r--r--  1 comp  staff  57955 Nov  3 17:27 hyades/hyades.res
	-rw-r--r--  1 comp  staff  50317 Nov  3 17:25 hyades/hyades.res.burnin

The .burnin files provide the sampling patterns during the burnin process and may be useful for diagnostic purposes, especially if singlePopMcmc is not sampling well (see below). The .res.burnin files look like:

::

	linux> head -2 hyades/hyades.res.burnin
	    logAge          Y        FeH    modulus absorption     logPost
	  8.821886   0.280626   0.086646  -0.010790   0.011206 -174.520334

And the .res files have the same format:

::

	linux> head -2 hyades/hyades.res
	    logAge          Y        FeH    modulus absorption     logPost
	  8.843553   0.288468   0.016795  -0.180626   0.013647 -198.263339

After the column headers, there is one record for each iteration of each of the cluster parameters of interest.

If everything goes well, all you really need to do is plot histograms for any column of interest. These are the posterior parameter distributions. You can also calculate moments of these columns if you’d like, and look at correlations among the columns, e.g. by plotting logAge vs. modulus.

******************************
D. sampleMass and sampleWDMass
******************************

These modules are useful for anyone interested in the masses of some or all of the stars in their database. Running them is unnecessary if you are only interested in the cluster parameters. The module sampleMass reports the primary mass and secondary mass ratio at all iterations for every star in the database, and sampleWDMass reports the primary mass for the subset of database stars that are being fit as WDs.

Running these programs is quite simple:

::

	linux> sampleWDMass
	Seed: 1690745648
	Warming up generator... Done.
	Generated 10000 values.
	Reading models... Done.

	sampledPars.at(0).age    = 8.78411
	sampledPars.at(last).age = 8.74765
	Part 2 completed successfully

Running sampleMass is effectively identical.

These output files names end with .wdMassSamples, .wdMassSamples.membership, .massSamples, and .massSamples.membership. These correspond to the WD mass outputs from sampleWDMass, the membership likelihood of those masses, the mass and secondary mass ratio outputs from sampleMass, and the membership likelihood of those pairs.

sampleWDMass output files consist of the same number of columns as there are WDs, and the same number of rows as there are in the results (.res) file. Each item in a row corresponds to the mass of a WD (ordered as in the database) given the sampled parameters in the results file. The membership file shares this format, though the values correspond to the likelihood that the given star is a member of a cluster with the given parameters.

sampleMass output files are similar to sampleWDMass output files but have two columns per star in the database. For every 0-indexed star 𝑘 in the database, column 2𝑘 corresponds to that star’s primary mass, and 2𝑘 + 1 to that star’s secondary mass ratio. The membership file is identical to that of sampleWDMass, though the values correspond to the likelihood that the given unresolved binary is a cluster member.

sampleWDMass has no configurable parameters.

sampleMass takes two parameters in the YAML file: deltaMass and deltaMassRatio. These values are used as starting step sizes for the adaptive MCMC process used to obtain mass and mass ratio. We recommended that you change these parameters only if you are manipulating the code for diagnostic purposes.

**********
E. makeCMD
**********

The final module of our software suite, makeCMD, is a small module that calculates a mean fit isochrone. This is helpful for runs that do not converge as well as for situations where the posterior distribution of some key parameter may be multimodal. To run makeCMD

::

	linux> makeCMD
	Seed: 1574116425
	Reading models... Done.
	***Warning: "F435W" is not available in the selected WD Atmosphere
	model
	This is non-fatal if you aren't using the WD models

The output of makeCMD looks like

::

	linux> head -2 hyades/hyades.cmd
	      Mass          U          B          V          R          I
	J          H          K      F435W      F475W      F550M      F555W
	F606W      F625W      F775W      F814W
	  0.150000  16.170454  14.623905  13.035403  11.950929  10.490719
	9.196773   8.640356   8.386555  14.666132  13.907510  12.768529
	13.128038  12.587469  12.243942  10.744691  10.477001


Because makeCMD uses the values of means under cluster in the base9.yaml, one can enter the mean or median values from the singlePopMcmc posterior distributions into the yaml file prior to running makeCMD. The output file from makeCMD can then be used to overplot what essentially amounts to the average fit isochrone from among the posterior parameter distributions. Note that this is not a best-fit isochrone, but rather a representative example drawn from that distribution. In fact, isochrones created from summary statistics such as mean or median parameters may not be truly representative if the distributions are substantially non- Gaussian because that simultaneous combination of parameters may fit the data with low probability.

**************
F. Hyades Test
**************

We have created a script, hyades.csh, which is set up to run on a Hyades data set (Hyades.UBV.testphot). It is a cshell script. If you have problems with this script, you may be using a shell other than the cshell or tshell, e.g. the Bourne shell. You can invoke the cshell as follows:

::

	Bourne shell> csh
	New prompt indicating you are now running csh> hyades.csh

This will allow you to test your code installation and plot results, then compare to the DeGennaro et al. results. Note that you will not obtain an exact correspondence to the results of DeGennaro et al. because we have updated the Hyades data set since that publication. Because of the relative depth of the Hyades, which is significant compared to its distance, we have now corrected the cluster stars to lie at the mean cluster distance using individual proper motions from Hipparcos and the cluster converging point method. Because of the way we have corrected distances, this data set is converted to absolute magnitude space (we otherwise always use apparent magnitudes) and for this one test case, you will find a distance modulus of approximately 0.0.

**********************************
G. How long does all of this take?
**********************************

In our tests, it took 147 minutes to run hyades.csh, which in turn ran singlePopMcmc for 152 Hyades stars in three photometric bands for 10,000 iterations on a early 2011 Macbook Pro (2.3 GHz Intel with 8 GB RAM) laptop computer. Increasing the number of filters or number of stars will increase the computation time linearly. Increasing the number of MCMC iterations will increase the run time, but somewhat less than linearly because some of the time is spent during the burnin. You will see substantial increases in runtime if you have much larger data sets and/or if you have to increase the total number of calculated iterations.