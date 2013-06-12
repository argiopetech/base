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
        double M_wd_up;
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

    struct SimClusterSettings
    {
        int nStars;
        int nFieldStars;
        int nBrownDwarfs;
        int percentBinary; // Fraction * 100
        int percentDB; // Fraction * 100
    };

    struct ScatterClusterSettings
    {
        int relevantFilt;
        double exposures[14];
        double brightLimit;
        double faintLimit;
        double limitS2N;
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

    struct Files
    {
        char *phot;
        char *output;
        char *scatter;
        char *config;
    };

    struct Settings
    {
        int seed;
        int verbose;

        struct Files files;
        struct MainSequenceSettings mainSequence;
        struct WhiteDwarfSettings whiteDwarf;
        struct BrownDwarfSettings brownDwarf;
        struct MpiMcmcSettings mpiMcmc;
        struct ClusterSettings cluster;
        struct SimClusterSettings simCluster;
        struct ScatterClusterSettings scatterCluster;
    };

    void makeSettings(char*, struct Settings*);
    void zeroSettingPointers(struct Settings*);
    void settingsFromCLI(int argc, char **argv, struct Settings *settings);

#ifdef __cplusplus
}

template <typename T> T getDefault(YAML::Node&, std::string&&, T);
template <typename T> T getOrDie(YAML::Node&, std::string&&);
YAML::Node getNode(YAML::Node &n, std::string &&f);
[[noreturn]] void exitWith (std::string&&);
#endif
#endif
