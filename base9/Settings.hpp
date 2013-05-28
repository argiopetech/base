#ifndef SETTINGS_H
#define SETTINGS_H

#ifdef __cplusplus
#include <string>

extern "C"
{
#endif
    struct MainSequenceSettings
    {
        int ifmr;
        int filterSet;
        int rgbModel;
    };

    struct WhiteDwarfSettings
    {
        int filterSet;
        double carbonicity;
    };

    struct BrownDwarfSettings
    {
        int filterSet;
    };

    struct MpiMcmcSettings
    {
        int burnIter;
        int maxIter;
        int thin;
    };

    struct ClusterSettings
    {
        double Fe_H;
        double Fe_HSigma;

        double distMod;
        double distModSigma;

        double Av;
        double AvSigma;

        double Y;
        double YSigma;

        double logClusAge;

        double minMag;
        double maxMag;
        int  index;
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
