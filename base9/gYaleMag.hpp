#if defined( GYALEMAG_H )

#else
#define GYALEMAG_H

#define N_YY_PARAMS        17
#define N_YY_Z             11
#define N_YY_AGES          41
#define N_YY_FILTS         8
#define MAX_YY_ENTRIES     140

#define  dydz       2.0
#define  yp         0.23
#define  zp         0.0
#define  Zsun       0.0181
#define  Xsun       0.7148705
#define  FeHa2     -0.217
#define  FeHa4     -0.470

#define  QUAD(x1,y1,x2,y2,x3,y3,x) (y1)*((x2)-(x))*((x3)-(x))/(((x2)-(x1))*((x3)-(x1))) \
    +(y2)*((x1)-(x))*((x3)-(x))/(((x1)-(x2))*((x3)-(x2)))               \
    +(y3)*((x1)-(x))*((x2)-(x))/(((x1)-(x3))*((x2)-(x3)))


#define CUBEINT(x1,y1,x2,y2,x3,y3,x4,y4,x) (x-x2)*(x-x3)*(x-x4)*y1/((x1-x2)*(x1-x3)*(x1-x4)) \
    +(x-x1)*(x-x3)*(x-x4)*y2/((x2-x1)*(x2-x3)*(x2-x4))                  \
    +(x-x1)*(x-x2)*(x-x4)*y3/((x3-x1)*(x3-x2)*(x3-x4))                  \
    +(x-x1)*(x-x2)*(x-x3)*y4/((x4-x1)*(x4-x2)*(x4-x3))

#define POLLIN(x1,y1,x2,y2,x) (x-x2)*y1/(x1-x2) +(x-x1)*y2/(x2-x1)

#define SQR(x) (x)*(x)

struct yyIsochrone
{
    double FeH;
    double age;
    double logAge;
    double z;
    //double y;
    double mass[MAX_YY_ENTRIES];
    int nEntries;
    double mag[MAX_YY_ENTRIES][N_YY_FILTS];
    double AgbTurnoffMass;
};

void loadYale (char *path, int filterSet);
double deriveYYAgbTip (double newFeH, double newY, double newLogAge);
double wdPrecLogAgeYY (double thisFeH, double thisY, double zamsMass);
double getYaleMags (double zamsMass);


#endif
