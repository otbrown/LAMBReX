#ifndef LAMBREX_H
#define LAMBREX_H

#include "Simulation.h"

Simulation * lambrexInit(int, char **, int, int, int, double, double, int (&)[NDIMS]);
void lambrexFinalize(Simulation *);

#endif
