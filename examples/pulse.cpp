#include <cstdio>
#include "lambrex.h"

int main (int argc, char * argv[])
{
  lambrexInit();

  printf("Hello world.\n");

  lambrexFinalize();

  return 0;
}
