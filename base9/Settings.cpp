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
    // test
    cout << "whiteDwarf.ifmr: " << whiteDwarf.ifmr << "\n";
    cout << "whiteDwarf.M_wd_up: " << whiteDwarf.M_wd_up << "\n";

    cout << "cluster.starting.Fe_H: " << cluster.starting.Fe_H << "\n";
    cout << "cluster.priorMeans.Fe_H: " << cluster.priorMeans.Fe_H << "\n";
    cout << "cluster.priorSigma.Fe_H: " << cluster.priorSigma.Fe_H << "\n";

    cout << "cluster.starting.distMod: " << cluster.starting.distMod << "\n";
    cout << "cluster.priorMeans.distMod: " << cluster.priorMeans.distMod << "\n";
    cout << "cluster.priorSigma.distMod: " << cluster.priorSigma.distMod << "\n";

    cout << "cluster.starting.Av: " << cluster.starting.Av << "\n";
    cout << "cluster.priorMeans.Av: " << cluster.priorMeans.Av << "\n";
    cout << "cluster.priorSigma.Av: " << cluster.priorSigma.Av << "\n";

    cout << "cluster.starting.Y: " << cluster.starting.Y << "\n";
    cout << "cluster.priorMeans.Y: " << cluster.priorMeans.Y << "\n";
    cout << "cluster.priorSigma.Y: " << cluster.priorSigma.Y << "\n";

    cout << "cluster.starting.carbonicity: " << cluster.starting.carbonicity << "\n";
    cout << "cluster.priorMeans.carbonicity: " << cluster.priorMeans.carbonicity << "\n";
    cout << "cluster.priorSigma.carbonicity: " << cluster.priorSigma.carbonicity << "\n";

    cout << "cluster.starting.logAge: " << cluster.starting.logAge << "\n";
    cout << "cluster.priorMeans.logAge: " << cluster.priorMeans.logAge << "\n";
    cout << "cluster.priorSigma.logAge: " << cluster.priorSigma.logAge << "\n";

    cout << "cluster.minMag: " << cluster.minMag << "\n";
    cout << "cluster.maxMag: " << cluster.maxMag << "\n";
    cout << "cluster.index: " << cluster.index << "\n";

    cout << "singlePopMcmc.burnIter: " << singlePopMcmc.burnIter << "\n";
    cout << "singlePopMcmc.stage3Iter: " << singlePopMcmc.stage3Iter << "\n";
    cout << "singlePopMcmc.maxIter: " << singlePopMcmc.maxIter << "\n";
    cout << "singlePopMcmc.thin: " << singlePopMcmc.thin << "\n";

    cout << "singlePopMcmc.adaptiveBigSteps: " << singlePopMcmc.adaptiveBigSteps << "\n";
    cout << "singlePopMcmc.trialIter: " << singlePopMcmc.trialIter << "\n";

    cout << "singlePopMcmc.stepSize[AGE]: " << singlePopMcmc.stepSize[AGE] << "\n";
    cout << "singlePopMcmc.stepSize[FEH]: " << singlePopMcmc.stepSize[FEH] << "\n";
    cout << "singlePopMcmc.stepSize[MOD]: " << singlePopMcmc.stepSize[MOD] << "\n";
    cout << "singlePopMcmc.stepSize[ABS]: " << singlePopMcmc.stepSize[ABS] << "\n";
    cout << "singlePopMcmc.stepSize[YYY]: " << singlePopMcmc.stepSize[YYY] << "\n";
    cout << "singlePopMcmc.stepSize[CARBONICITY]: " << singlePopMcmc.stepSize[CARBONICITY] << "\n";
    cout << "singlePopMcmc.stepSize[IFMR_INTERCEPT]: " << singlePopMcmc.stepSize[IFMR_INTERCEPT] << "\n";
    cout << "singlePopMcmc.stepSize[IFMR_SLOPE]: " << singlePopMcmc.stepSize[IFMR_SLOPE] << "\n";
    cout << "singlePopMcmc.stepSize[IFMR_QUADCOEF]: " << singlePopMcmc.stepSize[IFMR_QUADCOEF] << "\n";

    cout << "multiPopMcmc.YA_start: " << multiPopMcmc.YA_start << "\n";
    cout << "multiPopMcmc.YB_start: " << multiPopMcmc.YB_start << "\n";
    cout << "multiPopMcmc.lambda_start: " << multiPopMcmc.lambda_start << "\n";

    cout << "multiPopMcmc.YA_lo: " << multiPopMcmc.YA_lo << "\n";
    cout << "multiPopMcmc.YA_hi: " << multiPopMcmc.YA_hi << "\n";
    cout << "multiPopMcmc.YB_hi: " << multiPopMcmc.YB_hi << "\n";

    cout << "multiPopMcmc.lambdaStep: " << multiPopMcmc.lambdaStep << "\n";

    cout << "simCluster.nStars: " << simCluster.nStars << "\n";
    cout << "simCluster.percentBinary: " << simCluster.percentBinary << "\n";
    cout << "simCluster.percentDB: " << simCluster.percentDB << "\n";
    cout << "simCluster.nFieldStars: " << simCluster.nFieldStars << "\n";

    cout << "scatterCluster.brightLimit: " << scatterCluster.brightLimit << "\n";
    cout << "scatterCluster.faintLimit: " << scatterCluster.faintLimit << "\n";
    cout << "scatterCluster.relevantFilt: " << scatterCluster.relevantFilt << "\n";
    cout << "scatterCluster.limitS2N: " << scatterCluster.limitS2N << "\n";
    cout << "scatterCluster.crowded: " << scatterCluster.crowded << "\n";

    cout << "sampleMass.deltaMass: " << sampleMass.deltaMass << "\n";
    cout << "sampleMass.deltaMassRatio: " << sampleMass.deltaMassRatio << "\n";

    cout << "verbose: " << verbose << "\n";

    cout << "files.phot: " << files.phot << "\n";

    cout << "files.output: " << files.output << "\n";

    cout << "files.scatter: " << files.scatter << "\n";

    cout << "files.models: " << files.models << "\n";
    //
    whiteDwarf.ifmr = whiteDwarf.ifmr == 0 ? getOrRequest <int>(whiteDwarfNode, "ifmr") : whiteDwarf.ifmr;
    whiteDwarf.wdModel = static_cast<int>(whiteDwarf.wdModel) == 0 ? static_cast<WdModel>(getOrRequest <int>(whiteDwarfNode, "wdModel")) : whiteDwarf.wdModel;
    whiteDwarf.M_wd_up = whiteDwarf.M_wd_up == 0.0 ? getOrRequest <double>(whiteDwarfNode, "M_wd_up") : whiteDwarf.M_wd_up;
    
    cluster.starting.Fe_H = cluster.starting.Fe_H == 0.0 ? getOrRequest <double>(startingNode, "Fe_H") : cluster.starting.Fe_H;
    cluster.priorMeans.Fe_H = cluster.priorMeans.Fe_H == 0.0 ? getOrRequest <double>(meansNode, "Fe_H") : cluster.priorMeans.Fe_H;
    cluster.priorSigma.Fe_H = cluster.priorSigma.Fe_H == 0.0 ? getOrRequest <double>(sigmasNode, "Fe_H") : cluster.priorSigma.Fe_H;
    
    cluster.starting.distMod = cluster.starting.distMod == 0.0 ? getOrRequest <double>(startingNode, "distMod") : cluster.starting.distMod;
    cluster.priorMeans.distMod = cluster.priorMeans.distMod == 0.0 ? getOrRequest <double>(meansNode, "distMod") : cluster.priorMeans.distMod;
    cluster.priorSigma.distMod = cluster.priorSigma.distMod == 0.0 ? getOrRequest <double>(sigmasNode, "distMod") : cluster.priorSigma.distMod;
    
    cluster.starting.Av = cluster.starting.Av == 0.0 ? getOrRequest <double>(startingNode, "Av") : cluster.starting.Av;
    cluster.priorMeans.Av = cluster.priorMeans.Av == 0.0 ? getOrRequest <double>(meansNode, "Av") : cluster.priorMeans.Av;
    cluster.priorSigma.Av = cluster.priorSigma.Av == 0.0 ? getOrRequest <double>(sigmasNode, "Av") : cluster.priorSigma.Av;
    
    cluster.starting.Y = cluster.starting.Y == 0.0 ? getOrRequest <double>(startingNode, "Y") : cluster.starting.Y;
    cluster.priorMeans.Y = cluster.priorMeans.Y == 0.0 ? getOrRequest <double>(meansNode, "Y") : cluster.priorMeans.Y;
    cluster.priorSigma.Y = cluster.priorSigma.Y == 0.0 ? getOrRequest <double>(sigmasNode, "Y") : cluster.priorSigma.Y;
    
    cluster.starting.carbonicity = cluster.starting.carbonicity == 0.38 ? getOrRequest <double>(startingNode, "carbonicity") : cluster.starting.carbonicity;
    cluster.priorMeans.carbonicity = cluster.priorMeans.carbonicity == 0.38 ? getOrRequest <double>(meansNode, "carbonicity") : cluster.priorMeans.carbonicity;
    cluster.priorSigma.carbonicity = cluster.priorSigma.carbonicity == 0.0 ? getOrRequest <double>(sigmasNode, "carbonicity") : cluster.priorSigma.carbonicity;
    
    cluster.starting.logAge = cluster.starting.logAge == 0.0 ? getOrRequest <double>(startingNode, "logAge") : cluster.starting.logAge;
    cluster.priorMeans.logAge = cluster.priorMeans.logAge == 0.0 ? getOrRequest <double>(meansNode, "logAge") : cluster.priorMeans.logAge;
    cluster.priorSigma.logAge = cluster.priorSigma.logAge == 0.0 ? getOrRequest <double>(sigmasNode, "logAge") : cluster.priorSigma.logAge;
    
    cluster.minMag = cluster.minMag == 0.0 ? getOrRequest <double>(clusterNode, "minMag") : cluster.minMag;
    cluster.maxMag = cluster.maxMag == 0.0 ? getOrRequest <double>(clusterNode, "maxMag") : cluster.maxMag;
    cluster.index = cluster.index == 0 ? getOrRequest <int>(clusterNode, "index") : cluster.index;
    
    singlePopMcmc.burnIter = singlePopMcmc.burnIter == 0 ? getOrRequest <int>(singlePopConfNode, "stage2IterMax") : singlePopMcmc.burnIter;
    singlePopMcmc.stage3Iter = singlePopMcmc.stage3Iter == 0 ? getOrRequest <int>(singlePopConfNode, "stage3Iter") : singlePopMcmc.stage3Iter;
    singlePopMcmc.maxIter = singlePopMcmc.maxIter == 0 ? getOrRequest <int>(singlePopConfNode, "runIter") : singlePopMcmc.maxIter;
    singlePopMcmc.thin = singlePopMcmc.thin == 0 ? getOrRequest <int>(singlePopConfNode, "thin") : singlePopMcmc.thin;
    
    singlePopMcmc.adaptiveBigSteps = singlePopMcmc.adaptiveBigSteps == 0 ? getOrRequest <int>(mpiAdaptiveNode, "bigStepIter") : singlePopMcmc.adaptiveBigSteps;
    singlePopMcmc.trialIter = singlePopMcmc.trialIter == 0 ? getOrRequest <int>(mpiAdaptiveNode, "trialIter") : singlePopMcmc.trialIter;

    if (singlePopMcmc.trialIter <= 0)
        exitWith("singlePopMcmc:adaptive:trialIter must be greater than 0");
    
    singlePopMcmc.stepSize[AGE] = singlePopMcmc.stepSize[AGE] == 0.0 ? getOrRequest <double>(mpiStepNode, "age") : singlePopMcmc.stepSize[AGE];
    singlePopMcmc.stepSize[FEH] = singlePopMcmc.stepSize[FEH] == 0.0 ? getOrRequest <double>(mpiStepNode, "Fe_H") : singlePopMcmc.stepSize[FEH];
    singlePopMcmc.stepSize[MOD] = singlePopMcmc.stepSize[MOD] == 0.0 ? getOrRequest <double>(mpiStepNode, "distMod") : singlePopMcmc.stepSize[MOD];
    singlePopMcmc.stepSize[ABS] = singlePopMcmc.stepSize[ABS] == 0.0 ? getOrRequest <double>(mpiStepNode, "Av") : singlePopMcmc.stepSize[ABS];
    singlePopMcmc.stepSize[YYY] = singlePopMcmc.stepSize[YYY] == 0.0 ? getOrRequest <double>(mpiStepNode, "Y") : singlePopMcmc.stepSize[YYY];
    singlePopMcmc.stepSize[CARBONICITY] = singlePopMcmc.stepSize[CARBONICITY] == 0.0 ? getOrRequest <double>(mpiStepNode, "carbonicity") : singlePopMcmc.stepSize[CARBONICITY];
    singlePopMcmc.stepSize[IFMR_INTERCEPT] = singlePopMcmc.stepSize[IFMR_INTERCEPT] == 0.0 ? getOrRequest <double>(mpiStepNode, "ifmrIntercept") : singlePopMcmc.stepSize[IFMR_INTERCEPT];
    singlePopMcmc.stepSize[IFMR_SLOPE] = singlePopMcmc.stepSize[IFMR_SLOPE] == 0.0 ? getOrRequest <double>(mpiStepNode, "ifmrSlope") : singlePopMcmc.stepSize[IFMR_SLOPE];
    singlePopMcmc.stepSize[IFMR_QUADCOEF] = singlePopMcmc.stepSize[IFMR_QUADCOEF] == 0.0 ? getOrRequest <double>(mpiStepNode, "ifmrQuadCoef") : singlePopMcmc.stepSize[IFMR_QUADCOEF];
    
    multiPopMcmc.YA_start = multiPopMcmc.YA_start == 0.0 ? getOrRequest <double>(multiPopConfNode, "YA_start") : multiPopMcmc.YA_start;
    multiPopMcmc.YB_start = multiPopMcmc.YB_start == 0.0 ? getOrRequest <double>(multiPopConfNode, "YB_start") : multiPopMcmc.YB_start;
    multiPopMcmc.lambda_start = multiPopMcmc.lambda_start == 0.0 ? getOrRequest <double>(multiPopConfNode, "lambda_start") : multiPopMcmc.lambda_start;
    
    multiPopMcmc.YA_lo = multiPopMcmc.YA_lo == 0.0 ? getOrRequest <double>(multiPopConfNode, "YA_lo") : multiPopMcmc.YA_lo;
    multiPopMcmc.YA_hi = multiPopMcmc.YA_hi == 0.0 ? getOrRequest <double>(multiPopConfNode, "YA_hi") : multiPopMcmc.YA_hi;
    multiPopMcmc.YB_hi = multiPopMcmc.YB_hi == 0.0 ? getOrRequest <double>(multiPopConfNode, "YB_hi") : multiPopMcmc.YB_hi;
    
    multiPopMcmc.lambdaStep = multiPopMcmc.lambdaStep == 0.0 ? getOrRequest <double>(multiPopConfNode, "lambdaStep") : multiPopMcmc.lambdaStep;
    
    simCluster.nStars = simCluster.nStars == 0 ? getOrRequest <int>(simConfNode, "nStars") : simCluster.nStars;
    simCluster.percentBinary = simCluster.percentBinary == 0 ? getOrRequest <int>(simConfNode, "percentBinary") : simCluster.percentBinary;
    simCluster.percentDB = simCluster.percentDB == 0 ? getOrRequest <int>(simConfNode, "percentDB") : simCluster.percentDB;
    simCluster.nFieldStars = simCluster.nFieldStars == 0 ? getOrRequest <int>(simConfNode, "nFieldStars") : simCluster.nFieldStars;
//    simCluster.nBrownDwarfs = getOrRequest <int>(simConfNode, "nBrownDwarfs");
    
    scatterCluster.brightLimit = scatterCluster.brightLimit == 0.0 ? getOrRequest <double>(scatterConfNode, "brightLimit") : scatterCluster.brightLimit;
    scatterCluster.faintLimit = scatterCluster.faintLimit == 0.0 ? getOrRequest <double>(scatterConfNode, "faintLimit") : scatterCluster.faintLimit;
    scatterCluster.relevantFilt = scatterCluster.relevantFilt == 0 ? getOrRequest <int>(scatterConfNode, "relevantFilt") : scatterCluster.relevantFilt;
    scatterCluster.limitS2N = scatterCluster.limitS2N == 0.0 ? getOrRequest <double>(scatterConfNode, "limitS2N") : scatterCluster.limitS2N;
    scatterCluster.crowded = scatterCluster.crowded == false ? getOrRequest <bool>(scatterConfNode, "crowded") : scatterCluster.crowded;
    
    sampleMass.deltaMass = sampleMass.deltaMass == 0.0 ? getOrRequest <double>(sampleMassNode, "deltaMass") : sampleMass.deltaMass;
    sampleMass.deltaMassRatio = sampleMass.deltaMassRatio == 0.0 ? getOrRequest <double>(sampleMassNode, "deltaMassRatio") : sampleMass.deltaMassRatio;

    {
        auto tNode = getNode(scatterConfNode, "exposures");
        getOrRequest <double>(tNode, "U");
        scatterCluster.exposures = tNode.as<map<string, double>>();
    }
    
    verbose = verbose == 0 ? getOrRequest <int>(generalNode, "verbose") : verbose;
    
    // When we switch to C++11, we can change these to std::string and remove most of the cruft
    files.phot = files.phot == "" ? getOrRequest <string>(filesNode, "photFile") : files.phot;

    files.output = files.output == "" ? getOrRequest <string>(filesNode, "outputFileBase") : files.output;

    files.scatter = files.scatter == "" ? getOrRequest <string>(filesNode, "scatterFile") : files.scatter;

    files.models = files.models == "" ? getOrRequest <string>(filesNode, "modelDirectory") : files.models;
}

void Settings::fromCLI (int argc, char **argv)
{
    int i_noBinaries = 0; // Has to be set false
    int i_overrideBounds = 0; // Has to be set false
    int i_development = 0; // Has to be set false
    int i_noWDs = 0; // Has to be set false
    int i_details = 0; // Has to be set false

    char **t_argv = new char*[argc];

    for (int i = 0; i<argc; i++)
    {
        t_argv[i] = new char[strlen (argv[i]) + 1];

        strcpy (t_argv[i], argv[i]);
    }

    static struct option long_options[] = {
        // Thie one just sets a flag outright
        {"verbose", no_argument, &(verbose), 1},
        {"noBinaries", no_argument, &(i_noBinaries), 1},
        {"overrideBounds", no_argument, &(i_overrideBounds), 1},
        {"noWDs", no_argument, &(i_noWDs), 1},
        {"development", no_argument, &(i_development), 1},
        {"details", no_argument, &(i_details), 1},

        // These all have to be parsed
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
        {"priorlogAge", required_argument, 0, 0xF0},
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
        {"sigmalogAge", required_argument, 0, 0xC9},
        {"startingFe_H", required_argument, 0, 0xC8},
        {"startingdistMod", required_argument, 0, 0xC7},
        {"startingAv", required_argument, 0, 0xC6},
        {"startingY", required_argument, 0, 0xC5},
        {"startinglogAge", required_argument, 0, 0xC4},
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

    noBinaries = i_noBinaries;
    overrideBounds = i_overrideBounds;
    simCluster.noWDs = i_noWDs;
    development = i_development;
    details = i_details;

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

template<typename T> T Settings::getOrRequest (Node & n, string && f)
{
    if (n[f])
    {
        return n[f].as<T> ();
    }
    else
    {
        string input = "";
        cout << "Please enter your desired setting for field: '" + f + "'\n";
        getline(std::cin, input);
        T t;
        istringstream(input) >> t;
        return t;
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
    cerr << "\t--msRgbModel\t\t0 = Girardi\n\t\t\t\t1 = Chaboyer-Dotter w/He sampling\n\t\t\t\t2 = Yale-Yonsei\n\t\t\t\t3 = Old (jc2mass) DSED\n\t\t\t\t4 = New DSED" << endl << endl;
    cerr << "\t--ifmr\t\t\t0 = Weidemann\n\t\t\t\t1 = Williams\n\t\t\t\t2 = Salaris lin\n\t\t\t\t3 = Salaris pw lin\n\t\t\t\t4+ = tunable" << endl << endl;
    cerr << "\t--wdModel\t\t0 = Wood\n\t\t\t\t1 = Montgomery" << endl << endl;
    cerr << "\t--M_wd_up\t\tThe maximum mass for a WD-producing star" << endl << endl;
    cerr << "\t--bdModel\t\t0 = None\n\t\t\t\t1 = Baraffe" << endl << endl;
    cerr << "\t--startingFe_H" << endl;
    cerr << "\t--priorFe_H" << endl;
    cerr << "\t--sigmaFe_H" << endl << endl;
    cerr << "\t--startingdistMod" << endl;
    cerr << "\t--priordistMod" << endl;
    cerr << "\t--sigmadistMod" << endl << endl;
    cerr << "\t--startingAv" << endl;
    cerr << "\t--priorAv" << endl;
    cerr << "\t--sigmaAv" << endl << endl;
    cerr << "\t--startingY" << endl;
    cerr << "\t--priorY" << endl;
    cerr << "\t--sigmaY" << endl << endl;
    cerr << "\t--startingCarbonicity" << endl;
    cerr << "\t--priorCarbonicity" << endl;
    cerr << "\t--sigmaCarbonicity" << endl << endl;
    cerr << "\t--startinglogAge" << endl;
    cerr << "\t--priorlogAge" << endl;
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
    cerr << "=============" << endl;
    cerr << "\t--noWDs\t\t\tKeeps simCluster from generating WDs";
    cerr << endl;
}

static void printVersion()
{
    cerr << "Base-" << Base9_VERSION << endl;
}
