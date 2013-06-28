#ifndef MODEL_HPP
#define MODEL_HPP

/*** Define a structure model that houses information about the evolution model ***/
struct model
{
    int evoModel;
    int brownDwarfEvol;
    int mainSequenceEvol;
    int IFMR;
    int WDcooling;
    int WDatm;
    int filterSet;
    int numFilts;
};

void initModels (struct model *models);

#endif
