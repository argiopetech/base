#ifndef SETTINGS_H
#define SETTINGS_H

#ifdef __cplusplus
#include <string>

extern "C"
{
#endif
    struct MainSequenceSettings
    {
        int filterSet;
        int rgbModel;
    };

    struct WhiteDwarfSettings
    {
        int ifmr;
        int wdModel;
        double carbonicity;
    };

    struct BrownDwarfSettings
    {
        int bdModel;
    };

    struct MpiMcmcSettings
    {
        int burnIter;
        int maxIter;
        int thin;
    };

    struct MakeCMDSettings
    {
        double M_wd_up;
        int verbose;
        char *scatterFile;
    };

    struct ClusterSigmas
    {
        double Fe_H;
        double distMod;
        double Av;
        double Y;
    };

    struct ClusterSettings
    {
        double Fe_H;
        double distMod;
        double Av;
        double Y;

        struct ClusterSigmas sigma;

        double logClusAge;

        double minMag;
        double maxMag;
        int index;
    };

    struct Settings
    {
        int seed;
        char *photFile;
        char *outputFile;
        struct MainSequenceSettings mainSequence;
        struct WhiteDwarfSettings whiteDwarf;
        struct BrownDwarfSettings brownDwarf;
        struct MpiMcmcSettings mpiMcmc;
        struct MakeCMDSettings makeCMD;
        struct ClusterSettings cluster;
    };

    struct Settings* makeSettings(char*);

#ifdef __cplusplus
}

template <typename T> T getDefault(YAML::Node&, std::string&&, T);
template <typename T> T getOrDie(YAML::Node&, std::string&&);
YAML::Node getNode(YAML::Node &n, std::string &&f);
//[[noreturn]]
void exitWith (std::string&&);
#endif
#endif
