#ifndef GBARAFFEMAG_H
#define GBARAFFEMAG_H

#include <string>
#include <vector>

const int N_BAR_AGES   = 30;
const int N_BAR_MASSES = 24;
const int N_BAR_FILTS  = 6;

void loadBaraffe (std::string path);
void getBaraffeMags (const std::vector<int>&, double logAge, double mass);

#endif
