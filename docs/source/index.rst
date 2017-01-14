==========================================================================================
Bayesian Analysis for Stellar Evolution with Nine Parameters (BASE-9) v9.4.3 User’s Manual
==========================================================================================

* Ted von Hippel
	* Embry-Riddle Aeronautical University, Daytona Beach, FL, USA; ted.vonhippel@erau.edu
* Elliot Robinson
	* Argiope Technical Solutions, Ft White, FL, USA; elliot.robinson@argiopetech.com
* Elizabeth Jeffery
	* Brigham Young University, Provo, UT, USA; ejeffery@byu.edu
* Rachel Wagner-Kaiser
	* University of Florida, Gainesville, FL, USA; rawagnerkaiser@gmail.com
* Steven DeGennaro
	* Studio 42, Austin, TX, USA; studiofortytwo@yahoo.com
* Nathan Stein
	* University of Pennsylvania, Philadelphia, PA, USA; nathanmstein@gmail.com
* David Stenning
	* University of California, Irvine, CA, USA; dstennin@uci.edu
* William H Jefferys
	* University of Texas, Austin, TX, USA and University of Vermont, Burlington, VT, USA;  bill@astro.as.utexas.edu
* David van Dyk
	* Imperial College London, London, UK; d.van-dyk@imperial.ac.uk


BASE-9 is a Bayesian software suite that recovers star cluster and stellar parameters from photometry. BASE-9 is useful for analyzing single-age, single-metallicity star clusters, binaries, or single stars, and for simulating such systems. BASE-9 uses Markov chain Monte Carlo and brute-force numerical integration techniques to estimate the posterior probability distributions for the age, metallicity, helium abundance, distance modulus, and line-of-sight absorption for a cluster, and the mass, binary mass ratio, and cluster membership probability for every stellar object. BASE-9 is provided as open source code on a version-controlled web server. The executables are also available as Amazon Elastic Compute Cloud images. This manual provides potential users with an overview of BASE-9, including instructions for installation and use.

.. toctree::
   :maxdepth: 2

   Introduction
   Skip-the-install-and-go-to-the-cloud 
   Installation
   Running-BASE-9
   Diagnostics-of-run-quality
   Example-uses-of-BASE-9
   Modifying-the-code-to-extend-its-capabilities