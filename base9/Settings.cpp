#include <string>
#include <iostream>
#include <map>
#include <vector>

#include <cstdio>
#include <cstring>
#include <getopt.h>

#include "Base9Config.h"
#include "yaml-cpp/yaml.h"
#include "Settings.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::istringstream;
using std::string;
using std::vector;
using std::map;
using YAML::Node;
using YAML::LoadFile;

// Forward declaration
static void printUsage ();
static void printVersion ();


void Settings::loadSettings(int argc, char** argv, const string& defaultFile)
{
    fromCLI (argc, argv);

    if (!files.config.empty())
    {
        fromYaml (files.config);
    }
    else
    {
        fromYaml (defaultFile);
    }

    fromCLI (argc, argv);
}


void Settings::fromYaml (const string& yamlFile)
{
    Node configNode;

    try
    {
        configNode = LoadFile (yamlFile);
    }
    catch (YAML::BadFile e)
    {
        cerr << "Configuration file '" << yamlFile << "' not found." << endl;
        exit(1);
    }

    Node generalNode = getNode (configNode, "general");
    Node filesNode = getNode (generalNode, "files");
    Node mainSequenceNode = getNode (generalNode, "main_sequence");
    Node whiteDwarfNode = getNode (generalNode, "white_dwarfs");
    Node clusterNode = getNode (generalNode, "cluster");
    Node startingNode = getNode(clusterNode, "starting");
    Node priorNode = getNode(clusterNode, "priors");
    Node meansNode = getNode (priorNode, "means");
    Node sigmasNode = getNode (priorNode, "sigmas");
    Node singlePopConfNode = getNode (configNode, "singlePopMcmc");
    Node mpiAdaptiveNode = getNode(singlePopConfNode, "adaptive");
    Node mpiStepNode = getNode (singlePopConfNode, "stepSizes");
    Node multiPopConfNode = getNode (configNode, "multiPopMcmc");
    Node cmdConfNode = getNode (configNode, "makeCMD");
    Node simConfNode = getNode (configNode, "simCluster");
    Node isoConfNode = getNode (configNode, "makeIsochrone");
    Node scatterConfNode = getNode (configNode, "scatterCluster");
    Node sampleMassNode = getNode (configNode, "sampleMass");

    mainSequence.msRgbModel = static_cast<MsModel>(getOrRequest <int>(mainSequenceNode, "msRgbModel"));

    whiteDwarf.ifmr = getOrRequest <int>(whiteDwarfNode, "ifmr");
    whiteDwarf.wdModel = static_cast<WdModel>(getOrRequest <int>(whiteDwarfNode, "wdModel"));
    whiteDwarf.M_wd_up = getOrRequest <double>(whiteDwarfNode, "M_wd_up");

    cluster.starting.Fe_H = getOrRequest <double>(startingNode, "Fe_H");
    cluster.priorMeans.Fe_H = getOrRequest <double>(meansNode, "Fe_H");
    cluster.priorSigma.Fe_H = getOrRequest <double>(sigmasNode, "Fe_H");

    cluster.starting.distMod = getOrRequest <double>(startingNode, "distMod");
    cluster.priorMeans.distMod = getOrRequest <double>(meansNode, "distMod");
    cluster.priorSigma.distMod = getOrRequest <double>(sigmasNode, "distMod");

    cluster.starting.Av = getOrRequest <double>(startingNode, "Av");
    cluster.priorMeans.Av = getOrRequest <double>(meansNode, "Av");
    cluster.priorSigma.Av = getOrRequest <double>(sigmasNode, "Av");

    cluster.starting.Y = getOrRequest <double>(startingNode, "Y");
    cluster.priorMeans.Y = getOrRequest <double>(meansNode, "Y");
    cluster.priorSigma.Y = getOrRequest <double>(sigmasNode, "Y");

    cluster.starting.carbonicity = getOrRequest <double>(startingNode, "carbonicity");
    cluster.priorMeans.carbonicity = getOrRequest <double>(meansNode, "carbonicity");
    cluster.priorSigma.carbonicity = getOrRequest <double>(sigmasNode, "carbonicity");

    cluster.starting.logAge = getOrRequest <double>(startingNode, "logAge");
    cluster.priorMeans.logAge = getOrRequest <double>(meansNode, "logAge");
    cluster.priorSigma.logAge = getOrRequest <double>(sigmasNode, "logAge");

    cluster.minMag = getOrRequest <double>(clusterNode, "minMag");
    cluster.maxMag = getOrRequest <double>(clusterNode, "maxMag");
    cluster.index = getOrRequest <int>(clusterNode, "index");

    singlePopMcmc.burnIter = getOrRequest <int>(singlePopConfNode, "stage2IterMax");
    singlePopMcmc.stage3Iter = getOrRequest <int>(singlePopConfNode, "stage3Iter");
    singlePopMcmc.maxIter = getOrRequest <int>(singlePopConfNode, "runIter");
    singlePopMcmc.thin = getOrRequest <int>(singlePopConfNode, "thin");

    singlePopMcmc.adaptiveBigSteps = getOrRequest <int>(mpiAdaptiveNode, "bigStepIter");
    singlePopMcmc.trialIter = getOrRequest <int>(mpiAdaptiveNode, "trialIter");

    if (singlePopMcmc.trialIter <= 0)
        exitWith("singlePopMcmc:adaptive:trialIter must be greater than 0");

    singlePopMcmc.stepSize[AGE] = getOrRequest <double>(mpiStepNode, "age");
    singlePopMcmc.stepSize[FEH] = getOrRequest <double>(mpiStepNode, "Fe_H");
    singlePopMcmc.stepSize[MOD] = getOrRequest <double>(mpiStepNode, "distMod");
    singlePopMcmc.stepSize[ABS] = getOrRequest <double>(mpiStepNode, "Av");
    singlePopMcmc.stepSize[YYY] = getOrRequest <double>(mpiStepNode, "Y");
    singlePopMcmc.stepSize[CARBONICITY] = getOrRequest <double>(mpiStepNode, "carbonicity");
    singlePopMcmc.stepSize[IFMR_INTERCEPT] = getOrRequest <double>(mpiStepNode, "ifmrIntercept");
    singlePopMcmc.stepSize[IFMR_SLOPE] = getOrRequest <double>(mpiStepNode, "ifmrSlope");
    singlePopMcmc.stepSize[IFMR_QUADCOEF] = getOrRequest <double>(mpiStepNode, "ifmrQuadCoef");

    multiPopMcmc.YA_start = getOrRequest <double>(multiPopConfNode, "YA_start");
    multiPopMcmc.YB_start = getOrRequest <double>(multiPopConfNode, "YB_start");
    multiPopMcmc.lambda_start = getOrRequest <double>(multiPopConfNode, "lambda_start");

    multiPopMcmc.YA_lo = getOrRequest <double>(multiPopConfNode, "YA_lo");
    multiPopMcmc.YA_hi = getOrRequest <double>(multiPopConfNode, "YA_hi");
    multiPopMcmc.YB_hi = getOrRequest <double>(multiPopConfNode, "YB_hi");

    multiPopMcmc.lambdaStep = getOrRequest <double>(multiPopConfNode, "lambdaStep");

    simCluster.nStars = getOrRequest <int>(simConfNode, "nStars");
    simCluster.percentBinary = getOrRequest <int>(simConfNode, "percentBinary");
    simCluster.percentDB = getOrRequest <int>(simConfNode, "percentDB");
    simCluster.nFieldStars = getOrRequest <int>(simConfNode, "nFieldStars");
//    simCluster.nBrownDwarfs = getOrRequest <int>(simConfNode, "nBrownDwarfs");

    scatterCluster.brightLimit = getOrRequest <double>(scatterConfNode, "brightLimit");
    scatterCluster.faintLimit = getOrRequest <double>(scatterConfNode, "faintLimit");
    scatterCluster.relevantFilt = getOrRequest <int>(scatterConfNode, "relevantFilt");
    scatterCluster.limitS2N = getOrRequest <double>(scatterConfNode, "limitS2N");
    scatterCluster.crowded  = getOrRequest <bool>(scatterConfNode, "crowded");

    sampleMass.burnIters = getOrRequest <int>(sampleMassNode, "burnIters");
    sampleMass.iters     = getOrRequest <int>(sampleMassNode, "iters");
    sampleMass.deltaMass = getOrRequest <double>(sampleMassNode, "deltaMass");
    sampleMass.deltaMassRatio = getOrRequest <double>(sampleMassNode, "deltaMassRatio");

    {
        auto tNode = getNode(scatterConfNode, "exposures");
        getOrRequest <double>(tNode, "U");
        scatterCluster.exposures = tNode.as<map<string, double>>();
    }

    verbose = getOrRequest <int>(generalNode, "verbose");

    // When we switch to C++11, we can change these to std::string and remove most of the cruft
    files.backend = static_cast<Backend>(getOrRequest <int>(filesNode, "backend"));
    files.phot = getOrRequest <string>(filesNode, "photFile");
    files.output = getOrRequest <string>(filesNode, "outputFileBase");
    files.scatter = getOrRequest <string>(filesNode, "scatterFile");
    files.models = getOrRequest <string>(filesNode, "modelDirectory");
}

void Settings::fromCLI (int argc, char **argv)
{
    char **t_argv = new char*[argc];

    for (int i = 0; i<argc; i++)
    {
        t_argv[i] = new char[strlen (argv[i]) + 1];

        strcpy (t_argv[i], argv[i]);
    }

    static struct option long_options[] = {
        // These all have to be parsed
        {"msRgbModel", required_argument, 0, 0xFE},
        {"ifmr", required_argument, 0, 0xFD},
        {"wdModel", required_argument, 0, 0xFC},
        {"M_wd_up", required_argument, 0, 0xFA},
        {"bdModel", required_argument, 0, 0xF9},
        {"priorFe_H", required_argument, 0, 0xF8},
        {"sigmaFe_H", required_argument, 0, 0xF7},
        {"priorDistMod", required_argument, 0, 0xF6},
        {"sigmaDistMod", required_argument, 0, 0xF5},
        {"priorAv", required_argument, 0, 0xF4},
        {"sigmaAv", required_argument, 0, 0xF3},
        {"priorY", required_argument, 0, 0xF2},
        {"sigmaY", required_argument, 0, 0xF1},
        {"priorLogAge", required_argument, 0, 0xF0},
        {"minMag", required_argument, 0, 0xEF},
        {"maxMag", required_argument, 0, 0xEE},
        {"index", required_argument, 0, 0xED},
        {"burnIter", required_argument, 0, 0xEC},
        {"maxIter", required_argument, 0, 0xEB},
        {"thin", required_argument, 0, 0xEA},
        {"nStars", required_argument, 0, 0xE9},
        {"percentBinary", required_argument, 0, 0xE8},
        {"percentDB", required_argument, 0, 0xE7},
        {"nFieldStars", required_argument, 0, 0xE6},
        {"modelDirectory", required_argument, 0, 0xE5},
        {"brightLimit", required_argument, 0, 0xE4},
        {"faintLimit", required_argument, 0, 0xE3},
        {"relevantFilt", required_argument, 0, 0xE2},
        {"limitS2N", required_argument, 0, 0xE1},
        {"seed", required_argument, 0, 0xE0},
        {"photFile", required_argument, 0, 0xDF},
        {"scatterFile", required_argument, 0, 0xDE},
        {"outputFileBase", required_argument, 0, 0xDD},
        {"config", required_argument, 0, 0xDC},
        {"help", no_argument, 0, 0xDB},
        {"version", no_argument, 0, 0xDA},
        {"bigStepBurnin", no_argument, 0, 0xCF},
        {"threads", required_argument, 0, 0xCE},
        {"priorCarbonicity", required_argument, 0, 0xCD},
        {"sigmaCarbonicity", required_argument, 0, 0xCC},
        {"deltaMass", required_argument, 0, 0xCB},
        {"deltaMassRatio", required_argument, 0, 0xCA},
        {"sigmaLogAge", required_argument, 0, 0xC9},
        {"startingFe_H", required_argument, 0, 0xC8},
        {"startingDistMod", required_argument, 0, 0xC7},
        {"startingAv", required_argument, 0, 0xC6},
        {"startingY", required_argument, 0, 0xC5},
        {"startingLogAge", required_argument, 0, 0xC4},
        {"startingCarbonicity", required_argument, 0, 0xC3},
        {"backend", required_argument, 0, 0xC2},
        {"run", required_argument, 0, 0xC1},


        // Various flags
        // These are now handled the same way as the parameters due to occasional compiler weirdness
        {"verbose", no_argument, 0, 0xAF},
        {"noBinaries", no_argument, 0, 0xAE},
        {"overrideBounds", no_argument, 0, 0xAD},
        {"noWDs", no_argument, 0, 0xAC},
        {"development", no_argument, 0, 0xAB},
        {"details", no_argument, 0, 0xAA},
        {0, 0, 0, 0}
    };

    int c, option_index;

    optind = 0;

    while ((c = getopt_long (argc, t_argv, "", long_options, &option_index)) != (-1))
    {
        int i;

        switch (c)
        {
            case 0:
                /* If this option set a flag, do nothing else now. */
                if (long_options[option_index].flag != 0)
                    break;

                cout << "option " << long_options[option_index].name;

                if (optarg)
                    cout << " with arg " << optarg;

                cout << endl;
                break;

            case 0xFE:
                istringstream (string (optarg)) >> i;
                mainSequence.msRgbModel = static_cast<MsModel>(i);
                break;

            case 0xFD:
                istringstream (string (optarg)) >> whiteDwarf.ifmr;
                break;

            case 0xFC:
                istringstream (string (optarg)) >> i;
                whiteDwarf.wdModel = static_cast<WdModel>(i);
                break;

            case 0xC2:
                istringstream (string (optarg)) >> i;
                files.backend = static_cast<Backend>(i);
                break;

            case 0xFA:
                istringstream (string (optarg)) >> whiteDwarf.M_wd_up;
                break;

            case 0xF8:
                istringstream (string (optarg)) >> cluster.priorMeans.Fe_H;
                break;

            case 0xF7:
                istringstream (string (optarg)) >> cluster.priorSigma.Fe_H;
                break;

            case 0xF6:
                istringstream (string (optarg)) >> cluster.priorMeans.distMod;
                break;

            case 0xF5:
                istringstream (string (optarg)) >> cluster.priorSigma.distMod;
                break;

            case 0xF4:
                istringstream (string (optarg)) >> cluster.priorMeans.Av;
                break;

            case 0xF3:
                istringstream (string (optarg)) >> cluster.priorSigma.Av;
                break;

            case 0xF2:
                istringstream (string (optarg)) >> cluster.priorMeans.Y;
                break;

            case 0xF1:
                istringstream (string (optarg)) >> cluster.priorSigma.Y;
                break;

            case 0xCD:
                istringstream (string (optarg)) >> cluster.priorMeans.carbonicity;
                break;

            case 0xCC:
                istringstream (string (optarg)) >> cluster.priorSigma.carbonicity;
                break;

            case 0xF0:
                istringstream (string (optarg)) >> cluster.priorMeans.logAge;
                break;

            case 0xEF:
                istringstream (string (optarg)) >> cluster.minMag;
                break;

            case 0xEE:
                istringstream (string (optarg)) >> cluster.maxMag;
                break;

            case 0xED:
                istringstream (string (optarg)) >> cluster.index;
                break;

            case 0xEC:
                istringstream (string (optarg)) >> singlePopMcmc.burnIter;
                break;

            case 0xEB:
                istringstream (string (optarg)) >> singlePopMcmc.maxIter;
                break;

            case 0xEA:
                istringstream (string (optarg)) >> singlePopMcmc.thin;
                break;

            case 0xE9:
                istringstream (string (optarg)) >> simCluster.nStars;
                break;

            case 0xE8:
                istringstream (string (optarg)) >> simCluster.percentBinary;
                break;

            case 0xE7:
                istringstream (string (optarg)) >> simCluster.percentDB;
                break;

            case 0xE6:
                istringstream (string (optarg)) >> simCluster.nFieldStars;
                break;

            case 0xE5:
                files.models = optarg;
                break;

            case 0xE4:
                istringstream (string (optarg)) >> scatterCluster.brightLimit;
                break;

            case 0xE3:
                istringstream (string (optarg)) >> scatterCluster.faintLimit;
                break;

            case 0xE2:
                istringstream (string (optarg)) >> scatterCluster.relevantFilt;
                break;

            case 0xE1:
                istringstream (string (optarg)) >> scatterCluster.limitS2N;
                break;

            case 0xE0:
                istringstream (string (optarg)) >> seed;
                break;

            case 0xDF:
                files.phot = optarg;
                break;

            case 0xDE:
                files.scatter = optarg;
                break;

            case 0xDD:
                files.output = optarg;
                break;

            case 0xDC:
                files.config = optarg;
                break;

            case 0xDB:                  // --help
                printUsage ();
                exit (EXIT_SUCCESS);

            case 0xDA:
                printVersion ();
                exit (EXIT_SUCCESS);

            case 0xCF:
                singlePopMcmc.bigStepBurnin = true;
                break;

            case 0xCE:
                istringstream (string (optarg)) >> threads;
                if (threads <= 0)
                {
                    cerr << "You must have at least one thread" << endl;
                    exit(EXIT_FAILURE);
                }
                break;

            case 0xCB:
                istringstream (string (optarg)) >> sampleMass.deltaMass;
                break;

            case 0xCA:
                istringstream (string (optarg)) >> sampleMass.deltaMassRatio;
                break;

            case 0xC9:
                istringstream (string (optarg)) >> cluster.priorSigma.logAge;
                break;

            case 0xC8:
                istringstream (string (optarg)) >> cluster.starting.Fe_H;
                break;

            case 0xC7:
                istringstream (string (optarg)) >> cluster.starting.distMod;
                break;

            case 0xC6:
                istringstream (string (optarg)) >> cluster.starting.Av;
                break;

            case 0xC5:
                istringstream (string (optarg)) >> cluster.starting.Y;
                break;

            case 0xC4:
                istringstream (string (optarg)) >> cluster.starting.logAge;
                break;

            case 0xC3:
                istringstream (string (optarg)) >> cluster.starting.carbonicity;
                break;

            case 0xC1:
                istringstream (string (optarg)) >> run;
                break;

            case 0xAF:
                verbose = true;
                break;

            case 0xAE:
                noBinaries = true;
                break;

            case 0xAD:
                overrideBounds = true;
                break;

            case 0xAC:
                simCluster.noWDs = true;
                break;

            case 0xAB:
                development = true;
                break;

            case 0xAA:
                details = true;
                break;

            case '?':
                // getopt_long already printed an error message.
                printUsage ();
                exit (EXIT_FAILURE);

            default:
                abort ();
        }
    }

    for (int i = 0; i<argc; i++)
    {
        delete[] t_argv[i];
    }
    delete[] t_argv;

    // Print any remaining command line arguments (not options). This is mainly for debugging purposes.
    if (optind < argc)
    {
        cerr << "Unrecognized options: ";
        while (optind<argc)
            cerr << argv[optind++];
        cerr << endl;
    }
}

template<typename T>
T Settings::getDefault (Node & n, string f, T def)
{
    if (n[f])
    {
        return n[f].as<T> ();
    }
    else
    {
        return def;
    }
}

template<typename T>
T Settings::getOrDie (Node & n, string f)
{
    if (n[f])
    {
        return n[f].as<T> ();
    }
    else
    {
        exitWith ("Field '" + f + "' was not set");
    }
}

template<typename T>
T Settings::getOrRequest (Node & n, string f)
{
    // In the event that there is an issue with loading the YAML file, we
    // request more information from the user.
    //
    // There is no sensible way to request information from the user in some
    // circumstances (e.g., supercomputing clusters, where the client may be run
    // many nodes at once and where the user may not have access to stdin). For
    // these cases, we provide a build flag that makes this function call
    // getOrDie instead.
    #ifndef BUILD_NONINTERACTIVE

    if (n[f])
    {
        return n[f].as<T> ();
    }
    else
    {
        cerr << "No value was found in the YAML file for parameter '" + f + "'.\n Please enter your desired value for this parameter: ";

        string input = "";
        getline(std::cin, input);

        T t;
        istringstream(input) >> t;

        return t;
    }

    #else

    return getOrDie<T>(n, f);

    #endif
}

Node Settings::getNode (Node & n, string f)
{
    if (n[f])
    {
        return n[f];
    }
    else
    {
        exitWith ("Node '" + f + "' was not present\nIs your YAML file up to date?\n");
    }
}

[[noreturn]] void Settings::exitWith (string s)
{
    cerr << s << endl;
    abort ();
}

static void printUsage ()
{
    cerr << "\nUsage:" << endl;
    cerr <<   "======" << endl;
    cerr << "\t--help\t\t\tPrints help" << endl;
    cerr << "\t--version\t\tPrints version string\n" << endl;
    cerr << "\t--config\t\tYAML configuration file\n" << endl;
    cerr << "\t--msRgbModel\t\t0 = Girardi\n\t\t\t\t1 = Chaboyer-Dotter w/He sampling\n\t\t\t\t2 = Original Yale-Yonsei (not currently supported)\n\t\t\t\t3 = Old (jc2mass) DSED\n\t\t\t\t4 = New DSED\n\t\t\t\t5 = PARSEC/PanSTARRS\n\t\t\t\t6 = New Yale-Yonsei models (2018)\n" << endl;
    cerr << "\t--ifmr\t\t\t0 = Weidemann\n\t\t\t\t1 = Williams\n\t\t\t\t2 = Salaris lin\n\t\t\t\t3 = Salaris pw lin\n\t\t\t\t4+ = tunable\n" << endl;
    cerr << "\t--wdModel\t\t0 = Wood\n\t\t\t\t1 = Montgomery (original)\n\t\t\t\t2 = Althaus\n\t\t\t\t3 = Renedo\n\t\t\t\t4 = Montgomery (2018)\n" << endl;
    cerr << "\t--M_wd_up\t\tThe maximum mass for a WD-producing star\n" << endl;
    cerr << "\t--bdModel\t\t0 = None\n\t\t\t\t1 = Baraffe\n" << endl;
    cerr << "\t--startingFe_H" << endl;
    cerr << "\t--priorFe_H" << endl;
    cerr << "\t--sigmaFe_H\n" << endl;
    cerr << "\t--startingDistMod" << endl;
    cerr << "\t--priorDistMod" << endl;
    cerr << "\t--sigmaDistMod\n" << endl;
    cerr << "\t--startingAv" << endl;
    cerr << "\t--priorAv" << endl;
    cerr << "\t--sigmaAv\n" << endl;
    cerr << "\t--startingY" << endl;
    cerr << "\t--priorY" << endl;
    cerr << "\t--sigmaY\n" << endl;
    cerr << "\t--startingCarbonicity" << endl;
    cerr << "\t--priorCarbonicity" << endl;
    cerr << "\t--sigmaCarbonicity\n" << endl;
    cerr << "\t--startingLogAge" << endl;
    cerr << "\t--priorLogAge" << endl;
    cerr << "\t--sigmaLogAge" << endl;
    cerr << "\t--minMag" << endl;
    cerr << "\t--maxMag" << endl;
    cerr << "\t--index\t\t\t0 being the first filter in the dataset" << endl;
    cerr << "\t--burnIter" << endl;
    cerr << "\t--maxIter" << endl;
    cerr << "\t--thin" << endl;
    cerr << "\t--nStars" << endl;
    cerr << "\t--percentBinary\t\tpercent binaries (drawn randomly)" << endl;
    cerr << "\t--percentDB\t\tpercent of WDs that have He atmospheres (drawn randomly)" << endl;
    cerr << "\t--nFieldStars" << endl;
    cerr << "\t--brightLimit\t\tapparant mags, can remove bright stars, e.g. RGB" << endl;
    cerr << "\t--faintLimit\t\tapparant mags, can remove faint stars, e.g. faint MS and WDs" << endl;
    cerr << "\t--relevantFilt\t\t0=bluest band available" << endl;
    cerr << "\t--limitS2N\t\tuse to remove objects with overly low signal-to-noise" << endl;
    cerr << "\t--seed\t\t\tinitialize the random number generator" << endl;
    cerr << "\t--photFile" << endl;
    cerr << "\t--scatterFile" << endl;
    cerr << "\t--outputFileBase\tRun information is appended to this name" << endl;
    cerr << "\t--modelDirectory\tThe directory in which models are located\n" << endl;
    cerr << "\t--deltaMass\t\tSet the delta for primary mass in sampleMass" << endl;
    cerr << "\t--deltaMassRatio\tSet the delta for mass ratio in sampleMass" << endl;

    cerr << "\n9.3.0 flags" << endl;
    cerr <<   "===========" << endl;
    cerr << "\t--threads\t\tSpecify the number of local threads to run with" << endl;
    cerr << "\t--bigStepBurnin\t\tRun the burnin only using the \"propClustBigSteps\" algorithm" << endl;

    cerr << "\n9.4.0 flags" << endl;
    cerr <<   "===========" << endl;
    cerr << "\t--noBinaries\t\tTurns off integration over secondary mass" << endl;

    cerr << "\n9.4.4 flags" << endl;
    cerr <<   "===========" << endl;
    cerr << "\t--noWDs\t\t\tKeeps simCluster from generating WDs" << endl;

    cerr << "\n9.5.0 flags" << endl;
    cerr <<   "===========" << endl;
    cerr << "\t--backend <int>" << endl;
    cerr << "\t\tSpecify the desired back end:" << endl;
    cerr << "\t\t\t0 = File" << endl;
    cerr << "\t\t\t1 = SQLite" << endl;
    cerr << "\t--run <runID>" << endl;
    cerr << "\t\tSpecify a previous run ID in the DB on which to base this run." << endl;
    cerr << "\t\tCurrently only supported by sampleWDMass and sampleMass." << endl;
    cerr << "\t\tOnly makes sense if you're using the SQLite backend." << endl;
    cerr << endl;
}

static void printVersion()
{
    cerr << "Base-" << Base9_VERSION << endl;
}
