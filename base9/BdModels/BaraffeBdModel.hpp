#ifndef BARAFFEBDMODEL_HPP
#define BARAFFEBDMODEL_HPP

#include <array>
#include <string>
#include <vector>

#include "BdModel.hpp"

class BaraffeBdModel : public BdModel
{
    const int N_BAR_AGES   = 30;
    const int N_BAR_MASSES = 24;
    const int N_BAR_FILTS  = 6;

  public:
    void loadModel (std::string path);
    void getMags (const std::vector<int>&, std::array<double, FILTS>&, double logAge, double mass);
};
#endif
