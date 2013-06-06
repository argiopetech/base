#include <cstdio>
#include <cstring>
#include <getopt.h>
#include <iostream>
#include <memory>
#include <mpi.h>
#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>

#include "Settings.hpp"

using std::cerr;
using std::endl;
using std::istringstream;
using std::string;
using std::unique_ptr;
using std::vector;
using YAML::Node;
using YAML::LoadFile;

// Forward declaration
void printUsage();

void makeSettings(char *yamlFile, struct Settings *settings)
{
    Node config;

    try
    {
        config = LoadFile(yamlFile);
    }
    catch (YAML::BadFile e)
    {
        cerr << "Configuration file '" << yamlFile << "' not found." << endl;
        abort();
    }

    Node general = getNode(config, "general");
    Node files = getNode(general, "files");
    Node mainSequence = getNode(general, "main_sequence");
    Node whiteDwarfs = getNode(general, "white_dwarfs");
    Node brownDwarfs = getNode(general, "brown_dwarfs");
    Node cluster = getNode(general, "cluster");
    Node priors = getNode(cluster, "priors");
    Node sigmas = getNode(cluster, "sigmas");
    Node mpiConf = getNode(config, "mpiMcmc");
    Node cmdConf = getNode(config, "makeCMD");
    Node simConf = getNode(config, "simCluster");
    Node isoConf = getNode(config, "makeIsochrone");
    Node scatterConf = getNode(config, "scatterCluster");

    settings->mainSequence.filterSet = getOrDie<int>(mainSequence, "filterSet");
    settings->mainSequence.rgbModel = getOrDie<int>(mainSequence, "rgbModel");

    settings->whiteDwarf.ifmr = getOrDie<int>(whiteDwarfs, "ifmr");
    settings->whiteDwarf.wdModel = getOrDie<int>(whiteDwarfs, "wdModel");
    settings->whiteDwarf.carbonicity = getOrDie<double>(whiteDwarfs, "carbonicity");
    settings->whiteDwarf.M_wd_up = getOrDie<double>(whiteDwarfs, "M_wd_up");

    settings->brownDwarf.bdModel = getOrDie<int>(brownDwarfs, "bdModel");

    settings->cluster.Fe_H = getOrDie<double>(priors, "Fe_H");
    settings->cluster.sigma.Fe_H = getOrDie<double>(sigmas, "Fe_H");

    settings->cluster.distMod = getOrDie<double>(priors, "distMod");
    settings->cluster.sigma.distMod = getOrDie<double>(sigmas, "distMod");

    settings->cluster.Av = getOrDie<double>(priors, "Av");
    settings->cluster.sigma.Av = getOrDie<double>(sigmas, "Av");

    settings->cluster.Y = getOrDie<double>(priors, "Y");
    settings->cluster.sigma.Y = getOrDie<double>(sigmas, "Y");

    settings->cluster.logClusAge = getOrDie<double>(cluster, "logClusAge");

    settings->cluster.minMag = getOrDie<double>(cluster, "minMag");
    settings->cluster.maxMag = getOrDie<double>(cluster, "maxMag");
    settings->cluster.index = getOrDie<int>(cluster, "index");

    settings->mpiMcmc.burnIter = getOrDie<int>(mpiConf, "burnIter");
    settings->mpiMcmc.maxIter = getOrDie<int>(mpiConf, "runIter");
    settings->mpiMcmc.thin = getOrDie<int>(mpiConf, "thin");

    settings->simCluster.nStars = getOrDie<int>(simConf, "nStars");
    settings->simCluster.percentBinary = getOrDie<int>(simConf, "percentBinary");
    settings->simCluster.percentDB = getOrDie<int>(simConf, "percentDB");
    settings->simCluster.nFieldStars = getOrDie<int>(simConf, "nFieldStars");
    settings->simCluster.nBrownDwarfs = getOrDie<int>(simConf, "nBrownDwarfs");

    settings->scatterCluster.brightLimit = getOrDie<double>(scatterConf, "brightLimit");
    settings->scatterCluster.faintLimit = getOrDie<double>(scatterConf, "faintLimit");
    settings->scatterCluster.relevantFilt = getOrDie<int>(scatterConf, "relevantFilt");
    settings->scatterCluster.limitS2N = getOrDie<double>(scatterConf, "limitS2N");

    memcpy(settings->scatterCluster.exposures, static_cast<const void*>(getOrDie<vector<double>>(scatterConf, "exposures").data()), 14 * sizeof(double));

    settings->seed = getOrDie<int>(general, "seed");
    settings->verbose = getOrDie<int>(general, "verbose");

    // When we switch to C++11, we can change these to std::string and remove most of the cruft
    settings->files.phot = new char[100];
    strcpy(settings->files.phot, const_cast<char*>(getOrDie<string>(files, "photFile").c_str()));

    settings->files.output = new char[100];
    sprintf(settings->files.output, "%s.a%0.3f.s%d.m%d.maxMag%0.1f.modSigma%0.2f"
            , getOrDie<string>(files, "outputFileBase").c_str()
            , settings->cluster.logClusAge
            , settings->seed
            , settings->mainSequence.rgbModel
            , settings->cluster.maxMag
            , settings->cluster.sigma.distMod);

    settings->files.scatter = new char[100];
    strcpy(settings->files.scatter, const_cast<char*>(getOrDie<string>(files, "scatterFile").c_str()));
}

void settingsFromCLI(int argc, char **argv, struct Settings *settings)
{
    char *t_argv[argc];

    for (int i = 0; i < argc; i++)
    {
        t_argv[i] = new char[strlen(argv[i]) + 1];
        strcpy(t_argv[i], argv[i]);
    }

    static struct option long_options[] =
    {
        // Thie one just sets a flag outright
        {"verbose", no_argument, &(settings->verbose), 1},

        // These all have to be parsed
        {"filterSet",      required_argument, 0, 0xFF},
        {"rgbModel",       required_argument, 0, 0xFE},
        {"ifmr",           required_argument, 0, 0xFD},
        {"wdModel",        required_argument, 0, 0xFC},
        {"carbonicity",    required_argument, 0, 0xFB},
        {"M_wd_up",        required_argument, 0, 0xFA},
        {"bdModel",        required_argument, 0, 0xF9},
        {"priorFe_H",      required_argument, 0, 0xF8},
        {"sigmaFe_H",      required_argument, 0, 0xF7},
        {"priordistMod",   required_argument, 0, 0xF6},
        {"sigmadistMod",   required_argument, 0, 0xF5},
        {"priorAv",        required_argument, 0, 0xF4},
        {"sigmaAv",        required_argument, 0, 0xF3},
        {"priorY",         required_argument, 0, 0xF2},
        {"sigmaY",         required_argument, 0, 0xF1},
        {"logClusAge",     required_argument, 0, 0xF0},
        {"minMag",         required_argument, 0, 0xEF},
        {"maxMag",         required_argument, 0, 0xEE},
        {"index",          required_argument, 0, 0xED},
        {"burnIter",       required_argument, 0, 0xEC},
        {"maxIter",        required_argument, 0, 0xEB},
        {"thin",           required_argument, 0, 0xEA},
        {"nStars",         required_argument, 0, 0xE9},
        {"percentBinary",  required_argument, 0, 0xE8},
        {"percentDB",      required_argument, 0, 0xE7},
        {"nFieldStars",    required_argument, 0, 0xE6},
        {"nBrownDwarfs",   required_argument, 0, 0xE5},
        {"brightLimit",    required_argument, 0, 0xE4},
        {"faintLimit",     required_argument, 0, 0xE3},
        {"relevantFilt",   required_argument, 0, 0xE2},
        {"limitS2N",       required_argument, 0, 0xE1},
        {"seed",           required_argument, 0, 0xE0},
        {"photFile",       required_argument, 0, 0xDF},
        {"scatterFile",    required_argument, 0, 0xDE},
        {"outputFileBase", required_argument, 0, 0xDD},
        {"config",         required_argument, 0, 0xDC},
        {0, 0, 0, 0}
    };

    int c, option_index;

    optind = 0;

    while ((c = getopt_long (argc, t_argv, "", long_options, &option_index)) != (-1))
    {
        switch (c)
        {
            case 0:
                /* If this option set a flag, do nothing else now. */
                if (long_options[option_index].flag != 0)
                    break;
                printf ("option %s", long_options[option_index].name);
                if (optarg)
                    printf (" with arg %s", optarg);
                printf ("\n");
                break;

            case 0xFF:
                istringstream(string(optarg)) >> settings->mainSequence.filterSet;
                break;

            case 0xFE:
                istringstream(string(optarg)) >> settings->mainSequence.rgbModel;
                break;

            case 0xFD:
                istringstream(string(optarg)) >> settings->whiteDwarf.ifmr;
                break;

            case 0xFC:
                istringstream(string(optarg)) >> settings->whiteDwarf.wdModel;
                break;

            case 0xFB:
                istringstream(string(optarg)) >> settings->whiteDwarf.carbonicity;
                break;

            case 0xFA:
                istringstream(string(optarg)) >> settings->whiteDwarf.M_wd_up;
                break;

            case 0xF9:
                istringstream(string(optarg)) >> settings->brownDwarf.bdModel;
                break;

            case 0xF8:
                istringstream(string(optarg)) >> settings->cluster.Fe_H;
                break;

            case 0xF7:
                istringstream(string(optarg)) >> settings->cluster.sigma.Fe_H;
                break;

            case 0xF6:
                istringstream(string(optarg)) >> settings->cluster.distMod;
                break;

            case 0xF5:
                istringstream(string(optarg)) >> settings->cluster.sigma.distMod;
                break;

            case 0xF4:
                istringstream(string(optarg)) >> settings->cluster.Av;
                break;

            case 0xF3:
                istringstream(string(optarg)) >> settings->cluster.sigma.Av;
                break;

            case 0xF2:
                istringstream(string(optarg)) >> settings->cluster.Y;
                break;

            case 0xF1:
                istringstream(string(optarg)) >> settings->cluster.sigma.Y;
                break;

            case 0xF0:
                istringstream(string(optarg)) >> settings->cluster.logClusAge;
                break;

            case 0xEF:
                istringstream(string(optarg)) >> settings->cluster.minMag;
                break;

            case 0xEE:
                istringstream(string(optarg)) >> settings->cluster.maxMag;
                break;

            case 0xED:
                istringstream(string(optarg)) >> settings->cluster.index;
                break;

            case 0xEC:
                istringstream(string(optarg)) >> settings->mpiMcmc.burnIter;
                break;

            case 0xEB:
                istringstream(string(optarg)) >> settings->mpiMcmc.maxIter;
                break;

            case 0xEA:
                istringstream(string(optarg)) >> settings->mpiMcmc.thin;
                break;

            case 0xE9:
                istringstream(string(optarg)) >> settings->simCluster.nStars;
                break;

            case 0xE8:
                istringstream(string(optarg)) >> settings->simCluster.percentBinary;
                break;

            case 0xE7:
                istringstream(string(optarg)) >> settings->simCluster.percentDB;
                break;

            case 0xE6:
                istringstream(string(optarg)) >> settings->simCluster.nFieldStars;
                break;

            case 0xE5:
                istringstream(string(optarg)) >> settings->simCluster.nBrownDwarfs;
                break;

            case 0xE4:
                istringstream(string(optarg)) >> settings->scatterCluster.brightLimit;
                break;

            case 0xE3:
                istringstream(string(optarg)) >> settings->scatterCluster.faintLimit;
                break;

            case 0xE2:
                istringstream(string(optarg)) >> settings->scatterCluster.relevantFilt;
                break;

            case 0xE1:
                istringstream(string(optarg)) >> settings->scatterCluster.limitS2N;
                break;

            case 0xE0:
                istringstream(string(optarg)) >> settings->seed;
                break;

            case 0xDF:
                if (settings->files.phot) // Already something here...
                {
                    delete settings->files.phot;
                }

                assert(optarg); // This is a required parameter, so it should never be null

                settings->files.phot = new char[strlen(optarg) + 1];
                strcpy(settings->files.phot, optarg);
                break;

            case 0xDE:
                if (settings->files.scatter) // Already something here...
                {
                    delete settings->files.scatter;
                }

                assert(optarg); // This is a required parameter, so it should never be null

                settings->files.scatter = new char[strlen(optarg) + 1];
                strcpy(settings->files.scatter, optarg);
                break;

            case 0xDD:
                if (settings->files.output) // Already something here...
                {
                    delete settings->files.output;
                }

                assert(optarg); // This is a required parameter, so it should never be null

                settings->files.output = new char[100];
                sprintf(settings->files.output, "%s.a%0.3f.s%d.m%d.maxMag%0.1f.modSigma%0.2f"
                        , optarg
                        , settings->cluster.logClusAge
                        , settings->seed
                        , settings->mainSequence.rgbModel
                        , settings->cluster.maxMag
                        , settings->cluster.sigma.distMod);
                break;

            case 0xDC:
                assert(optarg); // This is a required parameter, so it should never be null

                settings->files.config = new char[100];
                strcpy(settings->files.config, optarg);
                break;

            case '?':
                // getopt_long already printed an error message.
                printUsage();
                MPI_Finalize();
                exit(EXIT_FAILURE);
                break;

            default:
                abort();
        }
    }

    // Print any remaining command line arguments (not options). This is mainly for debugging purposes.
    if (optind < argc)
    {
        cerr << "Unrecognized options: ";
        while (optind < argc)
            cerr << argv[optind++];
        cerr << endl;
    }
}

template <typename T> T getDefault(Node &n, string &&f, T def)
{
    if (n[f])
    {
        return n[f].as<T>();
    }
    else
    {
        return def;
    }
}

template <typename T> T getOrDie(Node &n, string &&f)
{
    if (n[f])
    {
        return n[f].as<T>();
    }
    else
    {
        exitWith("Field '" + f + "' was not set");
    }
}

Node getNode(Node &n, string &&f)
{
    if (n[f])
    {
        return n[f];
    }
    else
    {
        exitWith("Node '" + f + "' was not present\nIs your YAML file up to date?\n");
    }
}

[[noreturn]] void exitWith (string &&s)
{
    cerr << s << endl;
    abort();
}

void printUsage()
{
    int taskId;

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskId);

    if (taskId == 0)
    {
        cerr << "\nUsage:" << endl;
        cerr << "=======" << endl;
        cerr << "\t--filterSet\t\t0 = UBVRIJHK\n\t\t\t\t1 = ACS\n\t\t\t\t2 = SDSS + JHK" << endl << endl;
        cerr << "\t--rgbModel\t\t0 = Girardi\n\t\t\t\t1 = Chaboyer-Dotter w/He sampling\n\t\t\t\t2 = Yale-Yonsei\n\t\t\t\t3 = DSED" << endl << endl;
        cerr << "\t--ifmr\t\t\t0 = Weidemann\n\t\t\t\t1 = Williams\n\t\t\t\t2 = Salaris lin\n\t\t\t\t3 = Salaris pw lin\n\t\t\t\t4+ = tunable" << endl << endl;;
        cerr << "\t--wdModel\t\t0 = Wood\n\t\t\t\t1 = Montgomery" << endl << endl;
        cerr << "\t--carbonicity\t\t" << endl;
        cerr << "\t--M_wd_up\t\tThe maximum mass for a WD-producing star" << endl << endl;
        cerr << "\t--bdModel\t\t0 = None\n\t\t\t\t1 = Baraffe" << endl << endl;
        cerr << "\t--priorFe_H" << endl;
        cerr << "\t--sigmaFe_H" << endl << endl;;
        cerr << "\t--priordistMod" << endl;
        cerr << "\t--sigmadistMod" << endl << endl;;
        cerr << "\t--priorAv" << endl;
        cerr << "\t--sigmaAv" << endl << endl;;
        cerr << "\t--priorY" << endl;
        cerr << "\t--sigmaY" << endl << endl;;
        cerr << "\t--logClusAge" << endl;
        cerr << "\t--minMag" << endl;
        cerr << "\t--maxMag" << endl;
        cerr << "\t--index\t\t0 being the first filter in the dataset" << endl;
        cerr << "\t--burnIter" << endl;
        cerr << "\t--maxIter" << endl;
        cerr << "\t--thin" << endl;
        cerr << "\t--nStars" << endl;
        cerr << "\t--percentBinary\t\tpercent binaries (drawn randomly)" << endl;
        cerr << "\t--percentDB\t\tpercent of WDs that have He atmospheres (drawn randomly)" << endl;
        cerr << "\t--nFieldStars" << endl;
        cerr << "\t--nBrownDwarfs" << endl;
        cerr << "\t--brightLimit\t\tapparant mags, can remove bright stars, e.g. RGB" << endl;
        cerr << "\t--faintLimit\t\tapparant mags, can remove faint stars, e.g. faint MS and WDs" << endl;
        cerr << "\t--relevantFilt\t\t0=bluest band available" << endl;
        cerr << "\t--limitS2N\t\tuse to remove objects with overly low signal-to-noise" << endl;
        cerr << "\t--seed\t\t\tinitialize the random number generator" << endl;
        cerr << "\t--photFile" << endl;
        cerr << "\t--scatterFile" << endl;
        cerr << "\t--outputFileBase\t\tRun information is appended to this name" << endl;
        cerr << endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);
}
