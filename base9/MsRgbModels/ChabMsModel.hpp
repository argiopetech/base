#ifndef GCHABMAG_H
#define GCHABMAG_H

#include <string>
#include <vector>

#include "GenericMsModel.hpp"

const int N_CHAB_FILTS = 8;

class ChabMsModel : public GenericMsModel
{
  public:
    ChabMsModel() {;}
    virtual ~ChabMsModel() {;}

    virtual bool isSupported(FilterSetName filterSet) const
        { return filterSet == FilterSetName::UBVRIJHK; }

  protected:
    virtual int numFilts() const
        { return N_CHAB_FILTS; }

    virtual std::string getFileName (std::string path) const
        { return path + "chaboyer/chaboyer.model"; }
};

#endif
