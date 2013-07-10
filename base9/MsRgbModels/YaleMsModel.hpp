#ifndef GYALEMAG_H
#define GYALEMAG_H

#include <string>

#include "../MsRgbModel.hpp"

const int N_YY_PARAMS    = 17;
const int N_YY_Z         = 11;
const int N_YY_AGES      = 41;
const int N_YY_FILTS     = 8;
const int MAX_YY_ENTRIES = 140;

const double dydz  =  2.0;
const double yp    =  0.23;
const double zp    =  0.0;
const double Zsun  =  0.0181;
const double Xsun  =  0.7148705;
const double FeHa2 = -0.217;
const double FeHa4 = -0.470;

constexpr double QUAD(double x1, double y1, double x2, double y2, double x3, double y3, double x)
{
    return  y1 * (x2 - x) * (x3 - x) / ((x2 - x1) * (x3 - x1))       \
        + y2 * (x1 - x) * (x3 - x) / ((x1 - x2) * (x3 - x2))         \
        + y3 * (x1 - x) * (x2 - x) / ((x1 - x3) * (x2 - x3));
}

constexpr double CUBEINT(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, double x)
{
    return (x-x2)*(x-x3)*(x-x4)*y1/((x1-x2)*(x1-x3)*(x1-x4))            \
        +(x-x1)*(x-x3)*(x-x4)*y2/((x2-x1)*(x2-x3)*(x2-x4))              \
        +(x-x1)*(x-x2)*(x-x4)*y3/((x3-x1)*(x3-x2)*(x3-x4))              \
        +(x-x1)*(x-x2)*(x-x3)*y4/((x4-x1)*(x4-x2)*(x4-x3));
}

constexpr double POLLIN(double x1, double y1, double x2, double y2, double x)
{
    return (x - x2) * y1 / (x1-x2) + (x-x1) * y2 / (x2 - x1);
}

constexpr double SQR(double x)
{
    return x * x;
}

struct yyIsochrone
{
    double FeH;
    double age;
    double logAge;
    double z;
    double mass[MAX_YY_ENTRIES];
    int nEntries;
    double mag[MAX_YY_ENTRIES][N_YY_FILTS];
    double AgbTurnoffMass;
};

class YaleMsModel : public MsRgbModel
{
  public:
    YaleMsModel() {;}
    virtual ~YaleMsModel() {;}

    virtual double deriveAgbTipMass(const std::vector<int>&, double, double, double);
    virtual double msRgbEvol(const std::vector<int>&, std::array<double, FILTS>&, double);
    virtual double wdPrecLogAge(double, double, double);
    virtual void loadModel(std::string, MsFilter);

  private:
    void intpolAge(int, double);
};

#endif
