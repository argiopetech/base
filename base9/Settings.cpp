#include <array>
#include <string>
#include <iostream>
#include <vector>

#include <cstdio>
#include <cstring>
#include <getopt.h>
#include <iostream>

#include "Base9Config.h"
#include "yaml-cpp/yaml.h"
#include "Settings.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::istringstream;
using std::string;
using std::vector;
using std::array;
using YAML::Node;
using YAML::LoadFile;

// Forward declaration
static void printUsage ();
static void printVersion ();

void Settings::fromYaml (const string yamlFile)
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
    Node priorsNode = getNode (clusterNode, "priors");
    Node sigmasNode = getNode (clusterNode, "sigmas");
    Node mpiConfNode = getNode (configNode, "mpiMcmc");
    Node mpiAdaptiveNode = getNode(mpiConfNode, "adaptive");
    Node mpiStepNode = getNode (mpiConfNode, "stepSizes");
    Node cmdConfNode = getNode (configNode, "makeCMD");
    Node simConfNode = getNode (configNode, "simCluster");
    Node isoConfNode = getNode (configNode, "makeIsochrone");
    Node scatterConfNode = getNode (configNode, "scatterCluster");

    mainSequence.filterSet = static_cast<MsFilter>(getOrDie<int>(mainSequenceNode, "filterSet"));
    mainSequence.msRgbModel = static_cast<MsModel>(getOrDie<int>(mainSequenceNode, "msRgbModel"));

    whiteDwarf.ifmr = getOrDie<int>(whiteDwarfNode, "ifmr");
    whiteDwarf.wdModel = static_cast<WdModel>(getOrDie<int>(whiteDwarfNode, "wdModel"));
    whiteDwarf.M_wd_up = getOrDie<double>(whiteDwarfNode, "M_wd_up");

    cluster.Fe_H = getOrDie<double>(priorsNode, "Fe_H");
    cluster.sigma.Fe_H = getOrDie<double>(sigmasNode, "Fe_H");

    cluster.distMod = getOrDie<double>(priorsNode, "distMod");
    cluster.sigma.distMod = getOrDie<double>(sigmasNode, "distMod");

    cluster.Av = getOrDie<double>(priorsNode, "Av");
    cluster.sigma.Av = getOrDie<double>(sigmasNode, "Av");

    cluster.Y = getOrDie<double>(priorsNode, "Y");
    cluster.sigma.Y = getOrDie<double>(sigmasNode, "Y");

    cluster.carbonicity = getOrDie<double>(priorsNode, "carbonicity");
    cluster.sigma.carbonicity = getOrDie<double>(sigmasNode, "carbonicity");

    cluster.logClusAge = getOrDie<double>(clusterNode, "logClusAge");

    cluster.minMag = getOrDie<double>(clusterNode, "minMag");
    cluster.maxMag = getOrDie<double>(clusterNode, "maxMag");
    cluster.index = getOrDie<int>(clusterNode, "index");

    mpiMcmc.burnIter = getOrDie<int>(mpiConfNode, "burnIter");
    mpiMcmc.maxIter = getOrDie<int>(mpiConfNode, "runIter");
    mpiMcmc.thin = getOrDie<int>(mpiConfNode, "thin");

    mpiMcmc.adaptiveBigSteps = getOrDie<int>(mpiAdaptiveNode, "bigStepIter");
    mpiMcmc.trialIter = getOrDie<int>(mpiAdaptiveNode, "trialIter");

    if (mpiMcmc.trialIter <= 0)
        exitWith("mpiMcmc:adaptive:trialIter must be greater than 0");

    if (mpiMcmc.trialIter > mpiMcmc.adaptiveBigSteps)
        cerr << "(bigStepIter > trialIter): Are you sure this is what you want?" << endl;

    if (mpiMcmc.trialIter > mpiMcmc.burnIter)
        exitWith("trialIter must be greater than burnIter (may cause invalide covariance matrix)");

    mpiMcmc.stepSize[AGE] = getOrDie<double>(mpiStepNode, "age");
    mpiMcmc.stepSize[FEH] = getOrDie<double>(mpiStepNode, "Fe_H");
    mpiMcmc.stepSize[MOD] = getOrDie<double>(mpiStepNode, "distMod");
    mpiMcmc.stepSize[ABS] = getOrDie<double>(mpiStepNode, "Av");
    mpiMcmc.stepSize[YYY] = getOrDie<double>(mpiStepNode, "Y");
    mpiMcmc.stepSize[CARBONICITY] = getOrDie<double>(mpiStepNode, "carbonicity");
    mpiMcmc.stepSize[IFMR_INTERCEPT] = getOrDie<double>(mpiStepNode, "ifmrIntercept");
    mpiMcmc.stepSize[IFMR_SLOPE] = getOrDie<double>(mpiStepNode, "ifmrSlope");
    mpiMcmc.stepSize[IFMR_QUADCOEF] = getOrDie<double>(mpiStepNode, "ifmrQuadCoef");

    simCluster.nStars = getOrDie<int>(simConfNode, "nStars");
    simCluster.percentBinary = getOrDie<int>(simConfNode, "percentBinary");
    simCluster.percentDB = getOrDie<int>(simConfNode, "percentDB");
    simCluster.nFieldStars = getOrDie<int>(simConfNode, "nFieldStars");
//    simCluster.nBrownDwarfs = getOrDie<int>(simConfNode, "nBrownDwarfs");

    scatterCluster.brightLimit = getOrDie<double>(scatterConfNode, "brightLimit");
    scatterCluster.faintLimit = getOrDie<double>(scatterConfNode, "faintLimit");
    scatterCluster.relevantFilt = getOrDie<int>(scatterConfNode, "relevantFilt");
    scatterCluster.limitS2N = getOrDie<double>(scatterConfNode, "limitS2N");


    auto temp = getOrDie<vector<double>>(scatterConfNode, "exposures");
    std::copy(temp.begin(), temp.end(), scatterCluster.exposures.begin());

    seed = getOrDie<uint32_t>(generalNode, "seed");
    verbose = getOrDie<int>(generalNode, "verbose");

    // When we switch to C++11, we can change these to std::string and remove most of the cruft
    files.phot = getOrDie<string>(filesNode, "photFile");

    files.output = getOrDie<string>(filesNode, "outputFileBase");

    files.scatter = getOrDie<string>(filesNode, "scatterFile");

    files.models = getOrDie<string>(filesNode, "modelDirectory");
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
        // Thie one just sets a flag outright
        {"verbose", no_argument, &(verbose), 1},

        // These all have to be parsed
        {"filterSet", required_argument, 0, 0xFF},
        {"msRgbModel", required_argument, 0, 0xFE},
        {"ifmr", required_argument, 0, 0xFD},
        {"wdModel", required_argument, 0, 0xFC},
        {"M_wd_up", required_argument, 0, 0xFA},
        {"bdModel", required_argument, 0, 0xF9},
        {"priorFe_H", required_argument, 0, 0xF8},
        {"sigmaFe_H", required_argument, 0, 0xF7},
        {"priordistMod", required_argument, 0, 0xF6},
        {"sigmadistMod", required_argument, 0, 0xF5},
        {"priorAv", required_argument, 0, 0xF4},
        {"sigmaAv", required_argument, 0, 0xF3},
        {"priorY", required_argument, 0, 0xF2},
        {"sigmaY", required_argument, 0, 0xF1},
        {"priorCarbonicity", required_argument, 0, 0xCD},
        {"sigmaCarbonicity", required_argument, 0, 0xCC},
        {"logClusAge", required_argument, 0, 0xF0},
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

            case 0xFF:
                istringstream (string (optarg)) >> i;
                mainSequence.filterSet = static_cast<MsFilter>(i);
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

            case 0xFA:
                istringstream (string (optarg)) >> whiteDwarf.M_wd_up;
                break;

            case 0xF8:
                istringstream (string (optarg)) >> cluster.Fe_H;
                break;

            case 0xF7:
                istringstream (string (optarg)) >> cluster.sigma.Fe_H;
                break;

            case 0xF6:
                istringstream (string (optarg)) >> cluster.distMod;
                break;

            case 0xF5:
                istringstream (string (optarg)) >> cluster.sigma.distMod;
                break;

            case 0xF4:
                istringstream (string (optarg)) >> cluster.Av;
                break;

            case 0xF3:
                istringstream (string (optarg)) >> cluster.sigma.Av;
                break;

            case 0xF2:
                istringstream (string (optarg)) >> cluster.Y;
                break;

            case 0xF1:
                istringstream (string (optarg)) >> cluster.sigma.Y;
                break;

            case 0xCD:
                istringstream (string (optarg)) >> cluster.carbonicity;
                break;

            case 0xCC:
                istringstream (string (optarg)) >> cluster.sigma.carbonicity;
                break;

            case 0xF0:
                istringstream (string (optarg)) >> cluster.logClusAge;
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
                istringstream (string (optarg)) >> mpiMcmc.burnIter;
                break;

            case 0xEB:
                istringstream (string (optarg)) >> mpiMcmc.maxIter;
                break;

            case 0xEA:
                istringstream (string (optarg)) >> mpiMcmc.thin;
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
                mpiMcmc.bigStepBurnin = true;
                break;

            case 0xCE:
                istringstream (string (optarg)) >> threads;
                if (threads <= 0)
                {
                    cerr << "You must have at least one thread" << endl;
                    exit(EXIT_FAILURE);
                }
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

template<typename T> T Settings::getDefault (Node & n, string && f, T def)
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

template<typename T> T Settings::getOrDie (Node & n, string && f)
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

Node Settings::getNode (Node & n, string && f)
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

[[noreturn]] void Settings::exitWith (string && s)
{
    cerr << s << endl;
    abort ();
}

static void printUsage ()
{
    cerr << "\nUsage:" << endl;
    cerr << "=======" << endl;
    cerr << "\t--help\t\t\tPrints help" << endl;
    cerr << "\t--version\t\tPrints version string" << endl << endl;
    cerr << "\t--config\t\tYAML configuration file" << endl << endl;
    cerr << "\t--filterSet\t\t0 = UBVRIJHK\n\t\t\t\t1 = ACS\n\t\t\t\t2 = SDSS + JHK" << endl << endl;
    cerr << "\t--msRgbModel\t\t0 = Girardi\n\t\t\t\t1 = Chaboyer-Dotter w/He sampling\n\t\t\t\t2 = Yale-Yonsei\n\t\t\t\t3 = DSED" << endl << endl;
    cerr << "\t--ifmr\t\t\t0 = Weidemann\n\t\t\t\t1 = Williams\n\t\t\t\t2 = Salaris lin\n\t\t\t\t3 = Salaris pw lin\n\t\t\t\t4+ = tunable" << endl << endl;
    cerr << "\t--wdModel\t\t0 = Wood\n\t\t\t\t1 = Montgomery" << endl << endl;
    cerr << "\t--M_wd_up\t\tThe maximum mass for a WD-producing star" << endl << endl;
    cerr << "\t--bdModel\t\t0 = None\n\t\t\t\t1 = Baraffe" << endl << endl;
    cerr << "\t--priorFe_H" << endl;
    cerr << "\t--sigmaFe_H" << endl << endl;
    cerr << "\t--priordistMod" << endl;
    cerr << "\t--sigmadistMod" << endl << endl;
    cerr << "\t--priorAv" << endl;
    cerr << "\t--sigmaAv" << endl << endl;
    cerr << "\t--priorY" << endl;
    cerr << "\t--sigmaY" << endl << endl;
    cerr << "\t--priorCarbonicity" << endl;
    cerr << "\t--sigmaCarbonicity" << endl << endl;
    cerr << "\t--logClusAge" << endl;
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
    cerr << "\t--modelDirectory\tThe directory in which models are located" << endl;
    cerr << "\n9.3.0 flags" << endl;
    cerr <<   "===========" << endl;
    cerr << "\t--threads\t\tSpecify the number of local threads to run with" << endl;
    cerr << "\t--bigStepBurnin\t\tRun the burnin only using the \"propClustBigSteps\" algorithm" << endl;
    cerr << endl;
}

static void printVersion()
{
    cerr << "Base-" << Base9_VERSION_MAJOR << "." << Base9_VERSION_MINOR << "." << Base9_VERSION_PATCH << endl;
}
