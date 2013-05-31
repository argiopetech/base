#include <cstdio>
#include <cstring>
#include <iostream>
#include <memory>
#include <string>
#include <yaml-cpp/yaml.h>

#include "Settings.hpp"

using std::string;
using std::unique_ptr;
using YAML::Node;
using YAML::LoadFile;

struct Settings* makeSettings(char *yamlFile)
{
    unique_ptr<Settings> settings(new Settings);

    std::ostringstream oss;

    Node config = LoadFile(yamlFile);
    Node general = getNode(config, "general");
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

    settings->mainSequence.filterSet = getDefault<int>(mainSequence, "filterSet", 0);
    settings->mainSequence.rgbModel = getDefault<int>(mainSequence, "rgbModel", 2);

    settings->whiteDwarf.ifmr = getDefault<int>(whiteDwarfs, "ifmr", 1);
    settings->whiteDwarf.wdModel = getDefault<int>(whiteDwarfs, "wdModel", 1);
    settings->whiteDwarf.carbonicity = getDefault<double>(whiteDwarfs, "carbonicity", 0.7);

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

    settings->makeCMD.M_wd_up = getOrDie<double>(cmdConf, "M_wd_up");
    settings->makeCMD.verbose = getOrDie<int>(cmdConf, "verbose");
    settings->makeCMD.scatterFile = new char[100];
    strcpy(settings->makeCMD.scatterFile, const_cast<char*>(getOrDie<string>(cmdConf, "scatterFile").c_str()));

    settings->simCluster.nStars = getOrDie<int>(simConf, "nStars");
    settings->simCluster.verbose = getOrDie<int>(simConf, "verbose");
    settings->simCluster.M_wd_up = getOrDie<double>(simConf, "M_wd_up");
    settings->simCluster.fractionBinary = getOrDie<int>(simConf, "fractionBinary");
    settings->simCluster.fractionDB = getOrDie<int>(simConf, "fractionDB");
    settings->simCluster.nFieldStars = getOrDie<int>(simConf, "nFieldStars");
    settings->simCluster.nBrownDwarfs = getOrDie<int>(simConf, "nBrownDwarfs");

    settings->makeIso.M_wd_up = getOrDie<double>(isoConf, "M_wd_up");
    settings->makeIso.verbose = getOrDie<int>(isoConf, "verbose");

    settings->seed = getDefault<int>(general, "seed", 73);
    // When we switch to C++11, we can change these to std::string and remove most of the cruft
    settings->photFile = new char[100];
    strcpy(settings->photFile, const_cast<char*>(getOrDie<string>(general, "photFile").c_str()));
    settings->outputFile = new char[100];
    sprintf(settings->outputFile, "%s.a%0.3f.s%d.m%d.maxMag%0.1f.modSigma%0.2f", getOrDie<string>(general, "outputFileBase").c_str(), settings->cluster.logClusAge, settings->seed, settings->mainSequence.rgbModel, settings->cluster.maxMag, settings->cluster.sigma.distMod);

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
