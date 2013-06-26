#ifndef GBARAFFEMAG_H
#define GBARAFFEMAG_H

const int N_BAR_AGES   = 30;
const int N_BAR_MASSES = 24;
const int N_BAR_FILTS  = 6;

void loadBaraffe (char *path);
void getBaraffeMags (double logAge, double mass);

#endif
