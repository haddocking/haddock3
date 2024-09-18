#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

void outbuf_(void)
{
  (void) setvbuf(stdout, NULL, _IONBF, 0);
}

int csatty_(fildes)
     int *fildes;
{
  return isatty(*fildes);
}

