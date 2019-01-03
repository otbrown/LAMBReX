#ifndef LAMBREX_H
#define LAMBREX_H

#include "Lamb.h"

Lamb * lambrexInit(int, char **, int, int, int, double, double, int *);
void lambrexFinalize(Lamb *);

#endif
