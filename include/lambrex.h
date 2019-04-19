#ifndef LAMBREX_H
#define LAMBREX_H

#include "AmrSim.h"

void lambrexInit(const std::array<int,NDIMS>&);
void lambrexFinalise(void);
void lambrexSetAmr(const int, const int, const int, const int);
int binGCD(int, int);

#endif
