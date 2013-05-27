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
    Node mpiConf = getNode(config, "mpiMcmc");

    settings->mainSequence.ifmr = getDefault<int>(mainSequence, "ifmr", 0);
    settings->mainSequence.filterSet = getDefault<int>(mainSequence, "filterSet", 0);
    settings->mainSequence.rgbModel = getDefault<int>(mainSequence, "rgbModel", 0);

    settings->whiteDwarf.filterSet = getDefault<int>(whiteDwarfs, "filterSet", 0);
    settings->whiteDwarf.carbonicity = getDefault<double>(whiteDwarfs, "carbonicity", 0.6);

    settings->brownDwarf.filterSet = getDefault<int>(brownDwarfs, "filterSet", 1);

    settings->cluster.Fe_H = getOrDie<double>(cluster, "Fe_H");
    settings->cluster.Fe_HSigma = getOrDie<double>(cluster, "Fe_HSigma");

    settings->cluster.distMod = getOrDie<double>(cluster, "distMod");
    settings->cluster.distModSigma = getOrDie<double>(cluster, "distModSigma");

    settings->cluster.Av = getOrDie<double>(cluster, "Av");
    settings->cluster.AvSigma = getOrDie<double>(cluster, "AvSigma");

    settings->cluster.Y = getOrDie<double>(cluster, "Y");
    settings->cluster.YSigma = getOrDie<double>(cluster, "YSigma");

    settings->cluster.logClusAge = getOrDie<double>(cluster, "logClusAge");

    settings->cluster.minMag = getDefault<double>(cluster, "minMag", 0.0);
    settings->cluster.maxMag = getDefault<double>(cluster, "maxMag", 25.0);
    settings->cluster.index = getOrDie<int>(cluster, "index");

    settings->mpiMcmc.burnIter = getDefault<int>(mpiConf, "burnIter", 2000);
    settings->mpiMcmc.maxIter = getDefault<int>(mpiConf, "maxIter", 10000);
    settings->mpiMcmc.thin = getDefault<int>(mpiConf, "thin", 1);

    settings->seed = getDefault<int>(general, "seed", 73);
    // When we switch to C++11, we can change these to std::string and remove most of the cruft
    settings->photFile = new char[100];
    strcpy(settings->photFile, const_cast<char*>(getOrDie<string>(general, "photFile").c_str()));
    settings->outputFile = new char[100];
    sprintf(settings->outputFile, "%s.a%0.3f.s%d.m%d.maxMag%0.1f.modSigma%0.2f.res", getOrDie<string>(general, "outputFileBase").c_str(), settings->cluster.logClusAge, settings->seed, settings->mainSequence.rgbModel, settings->cluster.maxMag, settings->cluster.distModSigma);

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
