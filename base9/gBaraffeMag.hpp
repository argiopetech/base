#if defined( GBARAFFEMAG_H )

#else
#define GBARAFFEMAG_H

#define N_BAR_AGES        30
#define N_BAR_MASSES      24
#define N_BAR_FILTS        6

void loadBaraffe (char *path);
void getBaraffeMags (double logAge, double mass);

#endif
