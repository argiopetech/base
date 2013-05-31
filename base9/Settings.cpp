#include <cstdio>
#include <cstring>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>

#include "Settings.hpp"

using std::string;
using std::unique_ptr;
using std::vector;
using YAML::Node;
using YAML::LoadFile;

struct Settings* makeSettings(char *yamlFile)
{
    unique_ptr<Settings> settings(new Settings);

    std::ostringstream oss;

    Node config = LoadFile(yamlFile);
    Node general = getNode(config, "general");
    Node files = getNode(config, "files");
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

    settings->mainSequence.filterSet = getDefault<int>(mainSequence, "filterSet", 0);
    settings->mainSequence.rgbModel = getDefault<int>(mainSequence, "rgbModel", 2);

    settings->whiteDwarf.ifmr = getDefault<int>(whiteDwarfs, "ifmr", 1);
    settings->whiteDwarf.wdModel = getDefault<int>(whiteDwarfs, "wdModel", 1);
    settings->whiteDwarf.carbonicity = getDefault<double>(whiteDwarfs, "carbonicity", 0.7);
    settings->whiteDwarf.M_wd_up = getOrDie<double>(whiteDwarfs, "M_wd_up");

    settings->brownDwarf.bdModel = getDefault<int>(brownDwarfs, "bdModel", 1);

    settings->cluster.Fe_H = getOrDie<double>(priors, "Fe_H");
    settings->cluster.sigma.Fe_H = getOrDie<double>(sigmas, "Fe_H");

    settings->cluster.distMod = getOrDie<double>(priors, "distMod");
    settings->cluster.sigma.distMod = getOrDie<double>(sigmas, "distMod");

    settings->cluster.Av = getOrDie<double>(priors, "Av");
    settings->cluster.sigma.Av = getOrDie<double>(sigmas, "Av");

    settings->cluster.Y = getOrDie<double>(priors, "Y");
    settings->cluster.sigma.Y = getOrDie<double>(sigmas, "Y");

    settings->cluster.logClusAge = getOrDie<double>(cluster, "logClusAge");

    settings->cluster.minMag = getDefault<double>(cluster, "minMag", 0.0);
    settings->cluster.maxMag = getDefault<double>(cluster, "maxMag", 25.0);
    settings->cluster.index = getDefault<int>(cluster, "index", 2);

    settings->mpiMcmc.burnIter = getDefault<int>(mpiConf, "burnIter", 2000);
    settings->mpiMcmc.maxIter = getDefault<int>(mpiConf, "maxIter", 10000);
    settings->mpiMcmc.thin = getDefault<int>(mpiConf, "thin", 1);

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

    settings->seed = getDefault<int>(general, "seed", 73);
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


    return settings.release();
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
        exitWith("Field " + f + " was not set");
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
        exitWith("Node " + f + " was not present");
    }
}

//[[noreturn]]
void exitWith (string &&s)
{
    std::cout << s << std::endl;
    exit(EXIT_FAILURE);
}
