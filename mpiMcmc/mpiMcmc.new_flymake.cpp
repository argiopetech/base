int main (int argc, char *argv[])
{
  /// Load Settings
  {
    // Load CLI settings (do this first to get custom config files)

    // Load YAML settings

    // Reload CLI settings (override config file)
  }

  /// Initialize things
  {
    // Initialize RNG

    // Construct MS, WD, and BD models (generically)

    // Open output files with line buffering
  }

  /// Run simulation
  {
    /// Run Burnin
    {

    }

    /// Run MCMC
    {

    }
  }

}



/********* MAIN LOOP *********/
for (iteration = 0; iteration < ctrl.burnIter + ctrl.nIter * ctrl.thin; iteration++)
{
  propClust = mc.clust;

  /* propose and broadcast new value */
  propClustMarg (&propClust, &ctrl, iteration);
  logPostProp = logPriorClust (&propClust);

  index = 0;

  logPostProp = logPriorClust (&propClust);
  if (isInf (logPostProp))
  {
    /* don't bother computing, already know this cluster will be rejected */
    for (i = index; i < index + chunksize; i++)
    {
      logPostEachStar[i] = -HUGE_VAL;
    }
  }
  else
  {
    /* loop over assigned stars */
    for (i = index; i < index + chunksize; i++)
    {
      /* loop over all (mass1, mass ratio) pairs */
      if (mc.stars[i].status[0] == WD)
      {
	postClusterStar = 0.0;
	double tmpLogPost;

	for (j = 0; j < N_WD_MASS1; j++)
	{
	  wd[j] = mc.stars[i];
	  wd[j].boundsFlag = 0;
	  wd[j].isFieldStar = 0;
	  wd[j].U = wdMass1Grid[j];
	  wd[j].massRatio = 0.0;
	  evolve (&propClust, wd, j);

	  if (!wd[j].boundsFlag)
	  {
	    tmpLogPost = logPost1Star (&wd[j], &propClust);
	    tmpLogPost += log ((mc.clust.M_wd_up - MIN_MASS1) / (double) N_WD_MASS1);

	    postClusterStar += exp (tmpLogPost);
	  }
	}

	postClusterStar *= mc.stars[i].clustStarPriorDens;
      }
      else
      {
	/* marginalize over isochrone */
	postClusterStar = margEvolveWithBinary (&propClust, &mc.stars[i]);
	postClusterStar *= mc.stars[i].clustStarPriorDens;
      }

      /* marginalize over field star status */
      logPostEachStar[i] = log ((1.0 - mc.stars[i].clustStarPriorDens) * fsLike + postClusterStar);
    }
  }
  index = index + chunksize;
}

for (i = 0; i < mc.clust.nStars; i++)
{
  logPostProp += logPostEachStar[i];
}

/* accept/reject */
doAccept = acceptClustMarg (logPostCurr, logPostProp);
if (doAccept)
{
  mc.clust = propClust;
  logPostCurr = logPostProp;
  accept++;
}
else
{
  reject++;
}

/* save draws to estimate covariance matrix for more efficient Metropolis */
if (iteration >= ctrl.burnIter / 2 && iteration < ctrl.burnIter)
{
  if (iteration % increment == 0)
  {
    /* save draws */
    for (p = 0; p < NPARAMS; p++)
    {
      if (ctrl.priorVar[p] > EPSILON)
      {
	params[p][(iteration - ctrl.burnIter / 2) / increment] = mc.clust.parameter[p];
      }
    }
  }
  if (iteration == ctrl.burnIter - 1)
  {
    /* compute Cholesky decomposition of covariance matrix */
    int h, k;
    gsl_matrix *covMat = gsl_matrix_alloc (nParamsUsed, nParamsUsed);

    h = 0;

    double cholScale = 1000;	/* for numerical stability */

    for (i = 0; i < NPARAMS; i++)
    {
      if (ctrl.priorVar[i] > EPSILON)
      {
	k = 0;
	for (j = 0; j < NPARAMS; j++)
	{
	  if (ctrl.priorVar[j] > EPSILON)
	  {
	    cov = gsl_stats_covariance (params[i], 1, params[j], 1, nSave);
	    gsl_matrix_set (covMat, h, k, cov * cholScale * cholScale);	/* for numerical stability? */

	    if (h != k)
	    {
	      gsl_matrix_set (covMat, k, h, cov * cholScale * cholScale);
	    }

	    k++;
	  }
	}
	h++;
      }
    }

    for (i = 0; i < nParamsUsed; i++)
    {
      for (j = 0; j < nParamsUsed; j++)
      {
	printf ("%g ", gsl_matrix_get (covMat, i, j));
      }
      printf ("\n");
    }
    fflush (stdout);

    /* Cholesky decomposition */
    gsl_linalg_cholesky_decomp (covMat);

    /* compute proposal matrix from Cholesky factor */

    /* Gelman, Roberts, Gilks scale */
    double GRGscale = 0.97;	/* = 2.38 / sqrt(6) */

    h = 0;
    for (i = 0; i < NPARAMS; i++)
    {
      if (ctrl.priorVar[i] > EPSILON)
      {
	k = 0;
	for (j = 0; j < NPARAMS; j++)
	{
	  if (ctrl.priorVar[j] > EPSILON)
	  {
	    if (j <= i)
	    {
	      ctrl.propMatrix[i][j] = GRGscale * gsl_matrix_get (covMat, h, k) / cholScale;
	    }
	    else
	    {
	      ctrl.propMatrix[i][j] = 0.0;
	    }
	    k++;
	  }
	  else
	  {
	    ctrl.propMatrix[i][j] = 0.0;
	  }
	}
	h++;
      }
      else
      {
	for (j = 0; j < NPARAMS; j++)
	{
	  ctrl.propMatrix[i][j] = 0.0;
	}
      }
    }
  }

  /* Write output */
  if (iteration < ctrl.burnIter)
  {
    for (p = 0; p < NPARAMS; p++)	// For all parameters
    {
      if (ctrl.priorVar[p] > EPS || p == FEH || p == MOD || p == ABS)	// If
      {
	fprintf (ctrl.wClusterFile[1], "%10.6f ", mc.clust.parameter[p]);
      }
    }
    fprintf (ctrl.wClusterFile[1], "%10.6f\n", logPostCurr);
    fflush (ctrl.wClusterFile[1]);
  }
  else if (iteration % ctrl.thin == 0)
  {
    for (p = 0; p < NPARAMS; p++)
    {
      if (ctrl.priorVar[p] > EPS || p == FEH || p == MOD || p == ABS)
      {
	fprintf (ctrl.wClusterFile[0], "%10.6f ", mc.clust.parameter[p]);
      }
    }
    fprintf (ctrl.wClusterFile[0], "%10.6f\n", logPostCurr);
  }
}
